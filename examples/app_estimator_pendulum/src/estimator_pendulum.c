/**
 * ,---------,       ____  _ __
 * |  ,-^-,  |      / __ )(_) /_______________ _____  ___
 * | (  O  ) |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * | / ,--´  |    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *    +------`   /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2019 Bitcraze AB
 *
 * estimator_pendulum.c - Estimate pendulum angle when attached
 *
 * Estimator runs concurrently with the existing Kalman filter.
 * This OOT estimator (stabilizer.estimator = 3) implements OOT functions
 * that call the Kalman filter functions (estimator 2). This estimator
 * exposes theta and thetaDot via the pendEKF log group.
 */

#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "app.h"

#include "FreeRTOS.h"
#include "task.h"
#include "semphr.h"
#include "static_mem.h"

#define DEBUG_MODULE "ESTIMATORPENDULUM"
#include "debug.h"
#include "supervisor.h" // get quadIsFlying bool

// Reference
// https://www.bitcraze.io/2023/02/adding-an-estimator-or-controller/

#include "estimator_pendulum.h"   // typedefs and function headers and such
#include "estimator.h"            // this is our header file for the three Out Of Tree functions
#include "stabilizer_types.h"     // defines roll pitch yaw attitude struct
#include "estimator_kalman.h"     // use to run existing estimator underneath
// also use estimator_kalman.h to grab gyroLatest and accelLatest with custom accessors
#include "motors.h"               // motor IDs
#include "log.h"                  // logging
#include "system.h"               // systemWaitStart
#include "cf_math.h"              // for matrix math, includes arm_math.h
#include "pm.h"                   // power management, battery voltage compensation

// Public-facing pendulum state (wrapped + offset-corrected)
static pendulum_state_t pendulumEstimatorState; 

// Task declarations (EP = estimator pendulum)
static SemaphoreHandle_t runTaskSemaphoreEP;
static SemaphoreHandle_t dataMutexEP;
static StaticSemaphore_t dataMutexBufferEP;
static bool isInit = false;
static void pendulumTask(void* parameters);
STATIC_MEM_TASK_ALLOC_STACK_NO_DMA_CCM_SAFE(pendulumTask, PENDULUM_TASK_STACKSIZE);

// Faster than estimator_kalman
// Predict rate was originally RATE_100_HZ. I upped this to 250 Hz after Simulink testing
// 250 Hz was also next available option
#define PREDICT_RATE_PEND RATE_250_HZ
static const uint32_t PREDICTION_UPDATE_INTERVAL_MS_PEND = 1000 / PREDICT_RATE_PEND;

// Store after grabbing from estimator_kalman
static Axis3f accLatest;
static Axis3f gyroLatest;

// Logging
static float Fl_latest;
static float Fr_latest;

// Accelerometer filter
#define ACC_ALPHA 0.2f
static float acc_x_filtered;
static float acc_y_filtered;
static float acc_z_filtered;

// Force filter (unused)
#define FORCE_ALPHA 0.5f

// Phidot filter
#define PHI_D_ALPHA 0.1f
static float phi_d_filtered;

// Debug logs 
static volatile float flex1;
static volatile float flex2;
static volatile float flex3;

// debugging
//static int dbg = 0;

// Helper: wrap angle to [-pi, pi]
/*
#ifndef M_PI_F
#define M_PI_F 3.14159265f
#endif
static float wrapPi(float angle) {
  while (angle > M_PI_F) {
    angle -= 2.0f * M_PI_F;
  }
  while (angle < -M_PI_F) {
    angle += 2.0f * M_PI_F;
  }
  return angle;
}
*/

// Estimator data and parameter declarations
NO_DMA_CCM_SAFE_ZERO_INIT static pendulumCoreData_t pendulumCoreData;
static pendulumCoreParams_t pendulumCoreParams = {
  PENDULUM_DEFAULT_PARAMS_INIT
};
static helperConstants_t helperConstants;


/* ================================
 *   appMain – entry point
 * ================================ */

void appMain() {
  DEBUG_PRINT("Pendulum estimator app starting...\n");

  // Initialize pendulum task (background observer)
  //estimatorPendulumTaskInit(); // moved because main never seemed to run this

  while(1) {
    vTaskDelay(M2T(2000));
    //DEBUG_PRINT("Pendulum estimator alive\n");
  }
}


/* ========================================
 *   Out-of-tree estimator hook functions
 * ======================================== */

/**
 * (1) Initialize estimator - stacks on existing kalman filter outputs
 * 
 * Called by the estimator framework when CONFIG_ESTIMATOR_OOT is enabled.
 * We initialize the underlying Kalman filter and the pendulum core data here,
 * and we start the pendulum background task once.
 */
void estimatorOutOfTreeInit() {
  estimatorKalmanInit();

  uint32_t nowMs = T2M(xTaskGetTickCount());
  memset(&pendulumCoreData, 0, sizeof(pendulumCoreData_t));

  // Initialize states
  pendulumCoreData.S[THETA]     = pendulumCoreParams.initialTheta;
  pendulumCoreData.S[THETA_DOT] = pendulumCoreParams.initialThetaDot;

  // Initialize initial state variances
  pendulumCoreData.P[THETA][THETA] = pendulumCoreParams.stdDevInitialTheta*pendulumCoreParams.stdDevInitialTheta;
  pendulumCoreData.P[THETA_DOT][THETA_DOT] = pendulumCoreParams.stdDevInitialThetaDot*pendulumCoreParams.stdDevInitialThetaDot;

  pendulumCoreData.Pm.numRows = STATE_DIM;
  pendulumCoreData.Pm.numCols = STATE_DIM;
  pendulumCoreData.Pm.pData = (float*)pendulumCoreData.P; // don't forget this - avoid hard fault

  // Initialize Cm pointer
  pendulumCoreData.Cm.numRows = MEAS_DIM;
  pendulumCoreData.Cm.numCols = STATE_DIM;
  pendulumCoreData.Cm.pData = (float*)pendulumCoreData.C; // avoid hard fault

  // Initialize filter vars
  acc_x_filtered = 0.0f;
  acc_y_filtered = 0.0f;
  acc_z_filtered = 0.0f;
  phi_d_filtered = 0.0f;
  
  // Set update flag and time vals
  pendulumCoreData.isUpdated = false;
  pendulumCoreData.lastPredictionMs = nowMs;
  pendulumCoreData.lastProcessNoiseUpdateMs = nowMs;

  // Runtime initialization of helper constants
  float mq = pendulumCoreParams.mq;
  float mb = pendulumCoreParams.mb;
  float ms = pendulumCoreParams.ms;
  float Ixx = pendulumCoreParams.Ixx;
  float g = pendulumCoreParams.g;
  float r = pendulumCoreParams.r;
  float L = pendulumCoreParams.L;
  // For helper constants: 
  // A(1,0) = a1*cosf(theta)*(Fl + Fr)
  // In predict step: dt multiplied to A(0,1) and A(1,0)
  helperConstants.a1 = -(6.0f*(2.0f*mb + ms)) /
      (L*(12.0f*mb*mq + 4.0f*mb*ms + 4.0f*mq*ms + ms*ms));
  // B(1,0) = b1 - b2*sinf(theta) 
  // B(1,1) = -b1 - b2*sinf(theta)
  // In predict step: dt multiplied to B(1,0) and B(1,1)
  helperConstants.b1 = (L*ms*ms*r + 12.0f*L*mb*mq*r + 4.0f*L*mb*ms*r + 4.0f*L*mq*ms*r) /
      (Ixx*L*(12.0f*mb*mq + 4.0f*mb*ms + 4.0f*mq*ms + ms*ms));
  helperConstants.b2 = (12.0f*Ixx*mb + 6.0f*Ixx*ms) /
      (Ixx*L*(12.0f*mb*mq + 4.0f*mb*ms + 4.0f*mq*ms + ms*ms));
  // C(0,0) = c1 * [phi_d^2 + theta_d^2 + 2*phi_d*theta_d]*cosf(phi + theta) + c2 * [(Fl + Fr)*cosf(phi + 2.0*theta)]  
  // C(0,1) = c1 * 2*sinf(phi + theta)*(phi_d + theta_d) 
  // C(1,0) = c1 * [phi_d^2 + theta_d^2 + 2*phi_d*theta_d]*sinf(phi + theta) + c2 * [(Fl + Fr)*sinf(phi + 2.0*theta)]  
  // C(1,1) = -c1 * 2*cosf(phi + theta)*(phi_d + theta_d)
  helperConstants.c1 = L*(2.0f*mb + ms)/(2.0f*(mb + mq + ms));
  helperConstants.c2 = 3.0f*(2.0f*mb + ms)*(2.0f*mb + ms) /
      ((mb + mq + ms)*(12.0f*mb*mq + 4.0f*mb*ms + 4.0f*mq*ms + ms*ms));

  // Prediction step helper constants
  helperConstants.pred1 = L*r*(ms*ms - 12.0f*mb*mq - 4.0f*mb*ms - 4.0f*mq*ms);
  helperConstants.pred2 = L*(12.0f*mb*mq + 4.0f*mb*ms + 4.0f*mq*ms + ms*ms);
  
  // Correction step helper constants for nonlinear function
  helperConstants.cphi2 = 12.0f*mb*mb + 3.0f*ms*ms + 12.0f*mb*ms;
  helperConstants.cphi  = 12.0f*mb*mb + 5.0f*ms*ms + 24.0f*mb*mq +
                          20.0f*mb*ms + 8.0f*mq*ms;
  helperConstants.vphi2 = ms*ms*ms
        + 24.0f*mb*mb*mq
        + 6.0f*mb*ms*ms
        + 8.0f*mb*mb*ms
        + 4.0f*mq*ms*ms
        + 20.0f*mb*mq*ms; // V_phi2, velocity-related coefficient (multiply by L later)
  helperConstants.vphitheta = 2.0f*ms*ms*ms
              + 48.0f*mb*mb*mq
              + 12.0f*mb*ms*ms
              + 16.0f*mb*mb*ms
              + 8.0f*mq*ms*ms
              + 40.0f*mb*mq*ms; // V_phi_theta, velocity-related coefficient (multiply by L later)
  helperConstants.gblock =
      g*( 2.0f*ms*ms*ms
          + 24.0f*mb*mq*mq
          + 24.0f*mb*mb*mq
          + 10.0f*mb*ms*ms
          + 8.0f*mb*mb*ms
          + 10.0f*mq*ms*ms
          + 8.0f*mq*mq*ms
          + 40.0f*mb*mq*ms ); // G_block used only in yexp(2)
  helperConstants.expdenom =
      2.0f*(mb + mq + ms)*(12.0f*mb*mq + 4.0f*mb*ms + 4.0f*mq*ms + ms*ms);
  
  #if 1
  DEBUG_PRINT(
    "HELPER CONSTANTS "
    "a1=%.6f b1=%.6f b2=%.6f c1=%.6f c2=%.6f p1=%.9f p2=%.9f cp2=%.6f cp=%.6f vp2=%.6f vpt=%.6f gb=%.6f ed=%.6f\n",
    (double)helperConstants.a1,
    (double)helperConstants.b1,
    (double)helperConstants.b2,
    (double)helperConstants.c1,
    (double)helperConstants.c2,
    (double)helperConstants.pred1,
    (double)helperConstants.pred2,
    (double)helperConstants.cphi2,
    (double)helperConstants.cphi,
    (double)helperConstants.vphi2,
    (double)helperConstants.vphitheta,
    (double)helperConstants.gblock,
    (double)helperConstants.expdenom
  );
  #endif

  // Start the pendulum background task once
  if (!isInit) {
    estimatorPendulumTaskInit();
  }
}

/* (2) Required test function for estimator */

bool estimatorOutOfTreeTest() {
  return estimatorKalmanTest();
}

/**
 * (3) Estimator update function
 * 
 * When the OOT estimator is selected as stabilizer.estimator, this function is 
 * called. We call the Kalman estimator so it still externalizes its state
 * flight behavior is identical.
 */
void estimatorOutOfTree(state_t *state, const stabilizerStep_t stabilizerStep) {
  estimatorKalman(state, stabilizerStep); // base estimator

  if (runTaskSemaphoreEP != NULL) {
    xSemaphoreGive(runTaskSemaphoreEP);
  }
}


/* ==========================
 *   Pendulum task control
 * ========================== */

/* Starts the pendulum task */

void estimatorPendulumTaskInit() {
  // Binary semaphore used to wake the pendulum task from estimatorOutOfTree()
  runTaskSemaphoreEP = xSemaphoreCreateBinary();
  ASSERT(runTaskSemaphoreEP);

  // Mutex for shared state
  dataMutexEP = xSemaphoreCreateMutexStatic(&dataMutexBufferEP);

  // Create the pendulum RTOS task
  STATIC_MEM_TASK_CREATE(pendulumTask, pendulumTask, PENDULUM_TASK_NAME, NULL, PENDULUM_TASK_PRI); // changed

  isInit = true;
  DEBUG_PRINT("ESTIMATORPENDULUM: Task init complete\n");
}

/* Included for completion sake */

bool estimatorPendulumTaskTest() {
  return isInit;
}


/* ==========================
 *   Core predict/correct
 * ========================== */

/* Main prediction step - where half the magic happens */

void pendulumCorePredict(pendulumCoreData_t* this,
                         const pendulumCoreParams_t *params,
                         const Axis3f *gyro,
                         float Fl, float Fr,
                         const uint32_t nowMs) {
  float dt = (nowMs - this->lastPredictionMs) / 1000.0f;

  // The linearized update matrices
  NO_DMA_CCM_SAFE_ZERO_INIT static float A[STATE_DIM][STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Am = { STATE_DIM, STATE_DIM, (float *)A};

  NO_DMA_CCM_SAFE_ZERO_INIT static float B[STATE_DIM][INPUT_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Bm = { STATE_DIM, INPUT_DIM, (float *)B};

  NO_DMA_CCM_SAFE_ZERO_INIT static float C[STATE_DIM][INPUT_DIM];
  // C calculated here, but only used in correct step

  // Temp matrices for process noise Q calculation
  NO_DMA_CCM_SAFE_ZERO_INIT static float BT[INPUT_DIM][STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 BTm = { INPUT_DIM, STATE_DIM, (float *)BT};

  NO_DMA_CCM_SAFE_ZERO_INIT static float Q[STATE_DIM][STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Qm = { STATE_DIM, STATE_DIM, (float *)Q};

  // Temp matrices for the covariance updates
  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpPred1[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpPred1m = { STATE_DIM, STATE_DIM, tmpPred1};

  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpPred2[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpPred2m = { STATE_DIM, STATE_DIM, tmpPred2};

  // ====== DYNAMICS LINEARIZATION ======

  // Get local copy of state variables
  float theta = this->S[THETA];
  float theta_d = this->S[THETA_DOT];

  // Take phi and phi_d as parameters from underlying Kalman filter
  float R[3][3];
  estimatorKalmanGetEstimatedRot((float*)R);
  float phi = atan2f(R[2][1], R[2][2]); // roll [rad]
  float temp_phi_d = gyro->x * DEG_TO_RAD; // deg/s to rad/s!!!
  // filter phi_d
  #if 1
  phi_d_filtered = (PHI_D_ALPHA * temp_phi_d) + ((1.0f - PHI_D_ALPHA) * phi_d_filtered);
  temp_phi_d = phi_d_filtered;
  #endif
  float phi_d = temp_phi_d;

  // Store for correction step
  this->phi_hold  = phi;
  this->phi_d_hold = phi_d;

  A[0][0] = 1.0f;
  A[0][1] = dt;
  A[1][0] = dt * helperConstants.a1 * cosf(theta) * (Fl + Fr);
  A[1][1] = 1.0f;

  B[0][0] = 0.0f;
  B[0][1] = 0.0f;
  B[1][0] = dt * (helperConstants.b1 + helperConstants.b2 * sinf(theta));
  B[1][1] = dt * (helperConstants.b2 * sinf(theta) - helperConstants.b1);

  C[0][0] = helperConstants.c1 * (phi_d*phi_d + theta_d*theta_d + 2.0f*phi_d*theta_d)
              * cosf(phi + theta)
            + helperConstants.c2 * ((Fl + Fr) * cosf(phi + 2.0f*theta));
  C[0][1] = helperConstants.c1 * 2.0f * sinf(phi + theta) * (phi_d + theta_d);
  C[1][0] = helperConstants.c1 * (phi_d*phi_d + theta_d*theta_d + 2.0f*phi_d*theta_d)
              * sinf(phi + theta)
            + helperConstants.c2 * ((Fl + Fr) * sinf(phi + 2.0f*theta));
  C[1][1] = -helperConstants.c1 * 2.0f * cosf(phi + theta) * (phi_d + theta_d);

  // Store for correct step
  this->C[0][0] = C[0][0];
  this->C[0][1] = C[0][1];
  this->C[1][0] = C[1][0];
  this->C[1][1] = C[1][1];

  // ====== COVARIANCE UPDATE ======

  mat_mult(&Am, &this->Pm, &tmpPred1m); // A P
  mat_trans(&Am, &tmpPred2m); // A'
  mat_mult(&tmpPred1m, &tmpPred2m, &this->Pm); // A P A'

  // Calculate Q = B * B' * 0.0001
  mat_trans(&Bm, &BTm);                 // B'
  mat_mult(&Bm, &BTm, &Qm);             // B * B'
  mat_scale(&Qm, pendulumCoreParams.gamma, &Qm); // Scale by gamma = 0.00001
  #if 0
  if (++dbg % 200 == 0) {
    DEBUG_PRINT(
      "Q ENTRIES "
      "Q1=%.6f Q2=%.6f Q3=%.6f Q4=%.6f\n",
      (double)(Q[0][0]),
      (double)(Q[0][1]),
      (double)(Q[1][0]),
      (double)(Q[1][1]) // is 625*gamma with some spikes to 927*gamma
    );
  }
  #endif
  // Debug test: Overwrite Q 
  #if 1
  Q[0][0] = 0.0f * pendulumCoreParams.gamma;
  Q[0][1] = 0.0f;
  Q[1][0] = 0.0f;
  Q[1][1] = 625.0f * pendulumCoreParams.gamma;
  #endif

  // Add process noise: P = P + Q
  arm_mat_add_f32(&this->Pm, &Qm, &this->Pm);

  // ====== PREDICTION STEP ======

  this->S[THETA] += dt * this->S[THETA_DOT];
  //DEBUG_PRINT("DEBUG theta =%4.4f\n", (double)this->S[THETA]);

  #if 0
  this->S[THETA_DOT] += dt *
    - (Fl - Fr) *
      (-12.0f*pendulumCoreParams.Ixx*pendulumCoreParams.mb*sinf(theta)
       - 6.0f*pendulumCoreParams.Ixx*pendulumCoreParams.ms*sinf(theta)
       + helperConstants.pred1) /
      helperConstants.pred2;
  #endif

  this->S[THETA_DOT] += dt *
  - (Fl - Fr) * (-12.0f*pendulumCoreParams.Ixx*pendulumCoreParams.mb*sinf(theta)
  - 6.0f*pendulumCoreParams.Ixx*pendulumCoreParams.ms*sinf(theta) + helperConstants.pred1)
  / helperConstants.pred2 / pendulumCoreParams.Ixx;

  #if 0
  if (++dbg % 200 == 0) {
    float chunk1 = (-12.0f*pendulumCoreParams.Ixx*pendulumCoreParams.mb*sinf(theta) - 6.0f*pendulumCoreParams.Ixx*pendulumCoreParams.ms*sinf(theta) + helperConstants.pred1);
    float chunk2 = dt * -(Fl - Fr) * chunk1;
    float chunk3 = chunk2/helperConstants.pred2/pendulumCoreParams.Ixx;
    DEBUG_PRINT(
      "HELPER CONSTANTS "
      "Fl-Fr=%.9f dt=%.9f chunk1=%.9f chunk2=%.9f chunk3=%.9f thetaDot=%.9f\n",
      (double)(-(Fl-Fr)),
      (double)dt,
      (double)chunk1,
      (double)chunk2,
      (double)chunk3,
      (double)(this->S[THETA_DOT])
    );
  }
  #endif

  this->isUpdated = true;
  this->lastPredictionMs = nowMs;
}

/* Main correction step - where the other half of the magic happens */

void pendulumCoreCorrect(pendulumCoreData_t* this,
                         const pendulumCoreParams_t *params,
                         const Axis3f *acc,
                         float Fl, float Fr,
                         const uint32_t nowMs) {
  // The Kalman gain matrix
  NO_DMA_CCM_SAFE_ZERO_INIT static float L[STATE_DIM][MEAS_DIM];
  static arm_matrix_instance_f32 Lm = {STATE_DIM, MEAS_DIM, (float *)L};
  
  // The linearized update matrices
  NO_DMA_CCM_SAFE_ZERO_INIT static float C[MEAS_DIM][STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Cm = { MEAS_DIM, STATE_DIM, (float *)C};

  // Grab local C from data
  memcpy(C, this->C, sizeof(C));

  // Temporary matrices for math
  NO_DMA_CCM_SAFE_ZERO_INIT static float CT[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 CTm = { STATE_DIM, STATE_DIM, CT};

  NO_DMA_CCM_SAFE_ZERO_INIT static float PCT[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 PCTm = { STATE_DIM, STATE_DIM, PCT};

  NO_DMA_CCM_SAFE_ZERO_INIT static float CPCT[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 CPCTm = { STATE_DIM, STATE_DIM, CPCT};

  NO_DMA_CCM_SAFE_ZERO_INIT static float SmData[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Sm = { STATE_DIM, STATE_DIM, SmData};

  NO_DMA_CCM_SAFE_ZERO_INIT static float SinvmData[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Sinvm = { STATE_DIM, STATE_DIM, SinvmData};

  NO_DMA_CCM_SAFE_ZERO_INIT static float CP[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 CPm = { STATE_DIM, STATE_DIM, CP};

  NO_DMA_CCM_SAFE_ZERO_INIT static float LCP[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 LCPm = { STATE_DIM, STATE_DIM, LCP};

  // R covariance matrix
  NO_DMA_CCM_SAFE_ZERO_INIT static float RmData[MEAS_DIM][MEAS_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Rm = { MEAS_DIM, MEAS_DIM, (float *)RmData};
  
  RmData[0][0] = 1.0f * pendulumCoreParams.beta;
  RmData[0][1] = 0.0f;
  RmData[1][0] = 0.0f;
  RmData[1][1] = 1.0f * pendulumCoreParams.beta;

  // ====== GAIN UPDATE ======

  // Calculate the Kalman gain and perform the state update
  mat_trans(&Cm, &CTm);       // C'
  mat_mult(&this->Pm, &CTm, &PCTm);   // P C'
  mat_mult(&Cm, &PCTm, &CPCTm);         // C P C'
  arm_mat_add_f32(&Rm, &CPCTm, &Sm);  // R + C P C'
  mat_inv(&Sm, &Sinvm);               // (R + C P C')^-1
  mat_mult(&PCTm, &Sinvm, &Lm);         // L = P C' (R + C P C')^-1

  // ====== COVARIANCE UPDATE ======

  mat_mult(&Cm, &this->Pm, &CPm);        // C P
  mat_mult(&Lm, &CPm, &LCPm);         // L C P
  arm_mat_sub_f32(&this->Pm, &LCPm, &this->Pm); // P - L C P

  // ====== CALCULATE EXPECTED MEASUREMENT ====== 

  float theta   = this->S[THETA];
  float theta_d = this->S[THETA_DOT];
  float phi     = this->phi_hold;
  float phi_d   = this->phi_d_hold;

  float numerator_yexp1 =
      (Fl + Fr) * (helperConstants.cphi2 * sinf(phi + 2.0f*theta) - helperConstants.cphi  * sinf(phi))
    + pendulumCoreParams.L * sinf(phi + theta) *
      (helperConstants.vphi2 * (phi_d*phi_d + theta_d*theta_d) +
       helperConstants.vphitheta * (phi_d*theta_d));

  float numerator_yexp2 =
    - helperConstants.gblock
    + (Fl + Fr) * (helperConstants.cphi * cosf(phi) - helperConstants.cphi2 * cosf(phi + 2.0f*theta))
    - pendulumCoreParams.L * cosf(phi + theta) *
      (helperConstants.vphi2 * (phi_d*phi_d + theta_d*theta_d) +
       helperConstants.vphitheta * (phi_d*theta_d));

  float yexp1 = numerator_yexp1 / helperConstants.expdenom;
  float yexp2 = numerator_yexp2 / helperConstants.expdenom;

  // ====== CORRECT THE STATE ====== 

  float ymeas1 = acc->y;
  float ymeas2 = acc->z;

  // xhat = xbar + Lgain*(y_meas - yexp); 
  flex1 = ymeas1 - yexp1;
  flex2 = ymeas2 - yexp2;
  this->S[THETA]     += L[0][0] * (ymeas1 - yexp1) + L[0][1] * (ymeas2 - yexp2);
  this->S[THETA_DOT] += L[1][0] * (ymeas1 - yexp1) + L[1][1] * (ymeas2 - yexp2);

  #if 0
  if (++dbg % 200 == 0) {
    DEBUG_PRINT(
      "HELPER CONSTANTS "
      "theta=%.6f theta_dot=%.6f ymeas1=%.6f ymeas2=%.6f yexp1=%.6f yexp2=%.6f Fl=%.6f Fr=%.6f\n",
      (double)this->S[THETA],
      (double)this->S[THETA_DOT],
      (double)ymeas1,
      (double)ymeas2,
      (double)yexp1,
      (double)yexp2,
      (double)Fl,
      (double)Fr
    );
  }
  #endif

  this->isUpdated = true;
}

/**
 * Main task that runs the pendulum observer at PREDICT_RATE_PEND.
 * Calculates forces from motor PWM data, then calls predict and correct
 * steps. Finally exposes theta and theta_dot via pendulumEstimatorState.
 * 
 * Excluded rate supervisor, estimator supervisor (to check if in bounds), 
 * and statistics loggers for simplicity
 */
static void pendulumTask(void* parameters) {
  systemWaitStart();
  DEBUG_PRINT("ESTIMATORPENDULUM: pendulumTask started\n");

  uint32_t nowMs = T2M(xTaskGetTickCount());
  uint32_t nextPredictionMs = nowMs;

  // Debug variable
  //static uint32_t dbg = 0;

  while (true) {
    xSemaphoreTake(runTaskSemaphoreEP, portMAX_DELAY);
    uint32_t nowMs = T2M(xTaskGetTickCount());

    // Everything is inside this if-statement to run at PREDICT_RATE_PEND
    if (nowMs >= nextPredictionMs) {

      // ---- 1) Read ALL 4 motor PWM values as 8bit (from 16bit)----
      uint8_t pwm1_raw = (uint8_t) (motorsGetRatio(MOTOR_M1)>>8);
      uint8_t pwm2_raw = (uint8_t) (motorsGetRatio(MOTOR_M2)>>8);
      uint8_t pwm3_raw = (uint8_t) (motorsGetRatio(MOTOR_M3)>>8);
      uint8_t pwm4_raw = (uint8_t) (motorsGetRatio(MOTOR_M4)>>8);
      float pwm1 = (float)pwm1_raw; // avoid overflow errors with a uint8_t in next calculation
      float pwm2 = (float)pwm2_raw;
      float pwm3 = (float)pwm3_raw;
      float pwm4 = (float)pwm4_raw;

      // ---- 2) Convert PWM to Thrust in N from gram force ----
      // https://www.bitcraze.io/documentation/repository/crazyflie-firmware/master/functional-areas/pwm-to-thrust/
      // their test setup is for total thrust - we need to divide by 4
      float f1 = (0.000409f * pwm1 * pwm1 + 0.1405f * pwm1 - 0.099f) * 0.00980665f / 4;
      float f2 = (0.000409f * pwm2 * pwm2 + 0.1405f * pwm2 - 0.099f) * 0.00980665f / 4;
      float f3 = (0.000409f * pwm3 * pwm3 + 0.1405f * pwm3 - 0.099f) * 0.00980665f / 4;
      float f4 = (0.000409f * pwm4 * pwm4 + 0.1405f * pwm4 - 0.099f) * 0.00980665f / 4;

      // each f SHOULD be around or less than 0.20 N
      float exp = 0.78; // experimentally determined to compensate for battery voltage
      // 0.86 too high, 0.50 too low, 0.75 bit low, 0.80 too high,

      /* too noisy
      float bat = pmGetBatteryVoltage();
      float batMin = 3.00; // see charge curve - this is 0%
      float batMax = 4.10;
      exp = (bat - batMin)/(batMax - batMin);
      flex3 = exp;
      */

      float Fl = (f1 + f2)*exp;
      float Fr = (f3 + f4)*exp;
      Fl_latest = Fl; // N 
      Fr_latest = Fr; // N

      #if 0
      if (++dbg % 200 == 0) {
        DEBUG_PRINT(
          "Forces "
          "f1=%.6f f2=%.6f f3=%.6f f4=%.6f Fl=%.6f Fr=%.6f\n",
          (double)f1,
          (double)f2,
          (double)f3,
          (double)f4,
          (double)Fl_latest,
          (double)Fr_latest
        );
      }
      #endif

      // ---- 3) Grab gyro from Kalman filter via custom accessor function ----
      estimatorKalmanGetGyroLatest(&gyroLatest); // added to estimator_kalman.c

      // ---- 4) PREDICT STEP ----
      pendulumCorePredict(&pendulumCoreData, &pendulumCoreParams, &gyroLatest, Fl, Fr, nowMs);
      
      // ---- 5) Grab body acceleration from Kalman filter via custom accessor function ----
      Axis3f tempAccel;
      estimatorKalmanGetAccLatest(&tempAccel); // added to estimator_kalman.c
      // rotate to global frame - measurements need to correct in this frame
      // see kalman_core.c line 734
      float R[3][3];
      estimatorKalmanGetEstimatedRot((float*)R);
      float ax = tempAccel.x;
      float ay = tempAccel.y;
      float az = tempAccel.z;
      tempAccel.x = R[0][0]*ax + R[0][1]*ay + R[0][2]*az;
      tempAccel.y = R[1][0]*ax + R[1][1]*ay + R[1][2]*az;
      tempAccel.z = R[2][0]*ax + R[2][1]*ay + R[2][2]*az - 1; // sub g for coordinate acceleration
      // filter w/ LPF? (moving average) removed for now, acc is already filtered
      #if 0
      acc_x_filtered = (ACC_ALPHA * tempAccel.x) + ((1.0f - ACC_ALPHA) * acc_x_filtered);
      acc_y_filtered = (ACC_ALPHA * tempAccel.y) + ((1.0f - ACC_ALPHA) * acc_y_filtered);
      acc_z_filtered = (ACC_ALPHA * tempAccel.z) + ((1.0f - ACC_ALPHA) * acc_z_filtered);
      tempAccel.x = acc_x_filtered;
      tempAccel.y = acc_y_filtered;
      tempAccel.z = acc_z_filtered;
      #endif
      // convert from reading in g to m/s^2
      accLatest.x = tempAccel.x * pendulumCoreParams.g;
      accLatest.y = tempAccel.y * pendulumCoreParams.g;
      accLatest.z = tempAccel.z * pendulumCoreParams.g;
      
      // ---- 6) CORRECT STEP ----
      pendulumCoreCorrect(&pendulumCoreData, &pendulumCoreParams, &accLatest, Fl, Fr, nowMs);

      // --- Previous debugging: Simple theta update just so it moves ----
      // xSemaphoreTake(dataMutexEP, portMAX_DELAY);
      // pendulumEstimatorState.theta = pendulumCoreData.S[THETA];
      // pendulumEstimatorState.theta_dot = pendulumCoreData.S[THETA_DOT];
      // pendulumEstimatorState.timestamp = nowMs;
      // xSemaphoreGive(dataMutexEP);

      // ---- 7) Export theta in "ball down = 0" frame, wrapped to [-pi, pi] ----
      /*
      float theta_internal = pendulumCoreData.S[THETA];
      // Shift by the initialTheta (≈ pi) so that ball-down is 0 in the exported frame
      float theta_rel = theta_internal - pendulumCoreParams.initialTheta; // formerly pi
      theta_rel = wrapPi(theta_rel);
      */

      // ---- 8) Update public/log-facing copies ----
      xSemaphoreTake(dataMutexEP, portMAX_DELAY);
      pendulumEstimatorState.theta     = pendulumCoreData.S[THETA];//theta_rel;
      pendulumEstimatorState.theta_dot = pendulumCoreData.S[THETA_DOT];
      pendulumEstimatorState.timestamp = nowMs;
      xSemaphoreGive(dataMutexEP);

      // Update time to make next prediction
      nextPredictionMs = nowMs + PREDICTION_UPDATE_INTERVAL_MS_PEND;

      // ---- Debug print ----
      #if 0
      if (++dbg % 200 == 0) {
          DEBUG_PRINT(
            "ESTIMATORPENDULUM: ESTPEND V5 "
            "| g=[%.6f %.6f %.6f] "
            "| a=[%.6f %.6f %.6f] "
            "| pwm=[%f %f %f %f] "
            "| Fl=%.6f Fr=%.6f "
            "| theta=%.8f\n"
            "| thetaDot=%.8f\n",
            (double)gyroLatest.x, (double)gyroLatest.y, (double)gyroLatest.z,
            (double)accLatest.x,  (double)accLatest.y,  (double)accLatest.z,
            (double)f1, (double)f2, (double)f3, (double)f4,
            (double)Fl_latest, (double)Fr_latest,
            (double)pendulumEstimatorState.theta,
            (double)pendulumEstimatorState.theta_dot
          );
        }
      #endif
    }
  }
}


/* ==========================
 *   Getter + logging
 * ========================== */

void pendulumEstimatorGetState(pendulum_state_t *out) {
  xSemaphoreTake(dataMutexEP, portMAX_DELAY);
  memcpy(out, &pendulumEstimatorState, sizeof(pendulum_state_t));
  xSemaphoreGive(dataMutexEP);
}

LOG_GROUP_START(pendEKF)

  // States

  /** @brief State: theta (angle) */
  LOG_ADD(LOG_FLOAT, theta, &pendulumEstimatorState.theta)

  /** @brief State: theta_dot (angular velocity) */
  LOG_ADD(LOG_FLOAT, thetaDot, &pendulumEstimatorState.theta_dot)

  // State covariances
  
  /** @brief Var(theta) */
  LOG_ADD(LOG_FLOAT, varTheta, &pendulumCoreData.P[THETA][THETA])

  /** @brief Var(theta_dot) */
  LOG_ADD(LOG_FLOAT, varThetaDot, &pendulumCoreData.P[THETA_DOT][THETA_DOT])

  // Measurement accels in global frame
  
  /** @brief Measurement: x accel global (m/s^2) */
  LOG_ADD(LOG_FLOAT, xAccelGlobal , &accLatest.x)

  /** @brief Measurement: y accel global (m/s^2) */
  LOG_ADD(LOG_FLOAT, yAccelGlobal, &accLatest.y)

  /** @brief Measurement: z accel global (m/s^2) */
  LOG_ADD(LOG_FLOAT, zAccelGlobal, &accLatest.z)

  // Forces in body frame
  
  /** @brief Force: Left (N) */
  LOG_ADD(LOG_FLOAT, Fl, &Fl_latest)

  /** @brief Force: Right (N) */
  LOG_ADD(LOG_FLOAT, Fr, &Fr_latest)

  // Drone roll
  
  /** @brief Phi (rad) */
  LOG_ADD(LOG_FLOAT, phi, &pendulumCoreData.phi_hold)

  /** @brief Phi_dot (rad) */
  LOG_ADD(LOG_FLOAT, phiDot, &pendulumCoreData.phi_d_hold)

  // Flex debug
  
  /** @brief flex1 */
  LOG_ADD(LOG_FLOAT, flex_1, &flex1)

  /** @brief flex2 */
  LOG_ADD(LOG_FLOAT, flex_2, &flex2)

  /** @brief flex3 */
  LOG_ADD(LOG_FLOAT, flex_3, &flex3)


LOG_GROUP_STOP(pendEKF)
