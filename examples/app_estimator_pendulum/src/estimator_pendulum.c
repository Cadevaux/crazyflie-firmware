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
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * estimator_pendulum.c - Estimate pendulum angle when attached
 * 
 * Estimator can run concurrently with existing kalman filter by enabling
 * ESTIMATOR_OOT in menuconfig (additional selection). This defines
 * CONFIG_ESTIMATOR_OOT which we have =y in the app-config file
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
#include "motors.h"               // motor IDs
#include "sensors.h"              // sample accel and gyro
#include "log.h"                  // logging
#include "system.h"               // systemWaitStart

// Task declarations (EP = estimator pendulum)
static SemaphoreHandle_t runTaskSemaphoreEP;
static SemaphoreHandle_t dataMutexEP;
static StaticSemaphore_t dataMutexBufferEP;
static pendulum_state_t pendulumEstimatorState; // Name of the local state variable. Getter enables others to read
//static state_t taskEstimatorStateEP;  // state_t is for the entire state (pos, vel, att) and we don't wanna mess with this
static bool isInit = false;
static void pendulumTask(void* parameters);
STATIC_MEM_TASK_ALLOC_STACK_NO_DMA_CCM_SAFE(pendulumTask, PENDULUM_TASK_STACKSIZE);

// Declare these up top before they're used
static void pendulumCorePredict(pendulumCoreData_t* this,
                                const pendulumCoreParams_t *params,
                                const Axis3f *gyro,
                                float Fl, float Fr,
                                const uint32_t nowMs,
                                bool quadIsFlying);
static void pendulumCoreCorrect(pendulumCoreData_t* this,
                                const pendulumCoreParams_t *params,
                                const Axis3f *gyro,
                                const Axis3f *acc,
                                float Fl, float Fr,
                                const uint32_t nowMs,
                                bool quadIsFlying);

// Same as estimator_kalman.c
#define PREDICT_RATE RATE_250_HZ // this is slower than the IMU update rate of 1000Hz
// Predict rate was originally RATE_100_HZ. I upped this to 250 Hz after Simulink testing
// 250 Hz was also next available option
const uint32_t PREDICTION_UPDATE_INTERVAL_MS = 1000 / PREDICT_RATE;
static Axis3f accLatest;
static Axis3f gyroLatest;

// Estimator data and parameter declarations
NO_DMA_CCM_SAFE_ZERO_INIT static pendulumCoreData_t pendulumCoreData;
static pendulumCoreParams_t pendulumCoreParams = {
  PENDULUM_DEFAULT_PARAMS_INIT
};
static helperConstants_t helperConstants;

/* Main function - Initialize tasks, use for debugging */

void appMain() {
  DEBUG_PRINT("Waiting for activation ...\n");

  // Add the estimator task
  estimatorPendulumTaskInit(); 
  /*
  * estimatorKalmanTaskInit added kalmanTask already; estimatorKalmanTaskInit was called in
  * systemTask (system.c) since ESTIMATOR_KALMAN_ENABLE was selected in menuconfig.
  * ENABLE lets the calculations run in the background even if the estimator hasn't
  * been selected yet. We start to run both estimators in the init function (1) below
  */ 

  while(1) {
    vTaskDelay(M2T(2000));
    //DEBUG_PRINT("Hello World!\n");
  }
}

/* (1) Initialize estimator - stacks on existing kalman filter outputs */

void estimatorOutOfTreeInit() {
  estimatorKalmanInit(); // cascaded from

  uint32_t nowMs = T2M(xTaskGetTickCount());
  memset(&pendulumCoreData, 0, sizeof(pendulumCoreData_t));

  // Initialize states
  pendulumCoreData.S[THETA] = pendulumCoreParams.initialTheta;
  pendulumCoreData.S[THETA_DOT] = pendulumCoreParams.initialThetaDot;

  // Initialize initial state variances
  pendulumCoreData.P[THETA][THETA] = powf(pendulumCoreParams.stdDevInitialTheta, 2);
  pendulumCoreData.P[THETA_DOT][THETA_DOT] = powf(pendulumCoreParams.stdDevInitialThetaDot, 2);

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
  helperConstants.a1 = -(6*(2*mb + ms))/(L*(12*mb*mq + 4*mb*ms + 4*mq*ms + powf(ms,2))); 
  // B(1,0) = b1 - b2*sinf(theta) 
  // B(1,1) = -b1 - b2*sinf(theta)
  // In predict step: dt multiplied to B(1,0) and B(1,1)
  helperConstants.b1 = (L*ms*ms*r+12.0f*L*mb*mq*r+4.0f*L*mb*ms*r+4.0f*L*mq*ms*r)/(Ixx*L*(12.0f*mb*mq+4.0f*mb*ms+4.0f*mq*ms+ms*ms));
  helperConstants.b2 = (12.0f*Ixx*mb+6.0f*Ixx*ms)/(Ixx*L*(12.0f*mb*mq+4.0f*mb*ms+4.0f*mq*ms+ms*ms));
  // C(0,0) = c1 * [phi_d^2 + theta_d^2 + 2*phi_d*theta_d]*cosf(phi + theta) + c2 * [(Fl + Fr)*cosf(phi + 2.0*theta)]  
  // C(0,1) = c1 * 2*sinf(phi + theta)*(phi_d + theta_d) 
  // C(1,0) = c1 * [phi_d^2 + theta_d^2 + 2*phi_d*theta_d]*sinf(phi + theta) + c2 * [(Fl + Fr)*sinf(phi + 2.0*theta)]  
  // C(1,1) = -c1 * 2*cosf(phi + theta)*(phi_d + theta_d) 
  helperConstants.c1 = L*(2*mb + ms)/(2*(mb + mq + ms));
  helperConstants.c2 = 3*powf(2*mb + ms,2)/((mb + mq + ms)*(12*mb*mq + 4*mb*ms + 4*mq*ms + powf(ms,2)));

  // Prediction step helper constants
  helperConstants.pred1 = L*ms*ms*r - 12*L*mb*mq*r - 4*L*mb*ms*r - 4*L*mq*ms*r;
  helperConstants.pred2 = Ixx*L*(12*mb*mq + 4*mb*ms + 4*mq*ms + ms*ms);

  // Correction step helper constants for nonlinear function
  helperConstants.cphi2 = 12*mb*mb + 3*ms*ms + 12*mb*ms; // C_phi2
  helperConstants.cphi = 12*mb*mb + 5*ms*ms + 24*mb*mq + 20*mb*ms + 8*mq*ms; // C_phi
  helperConstants.vphi2 = ms*ms*ms
        + 24*mb*mb*mq
        + 6*mb*ms*ms
        + 8*mb*mb*ms
        + 4*mq*ms*ms
        + 20*mb*mq*ms; // V_phi2, velocity-related coefficient (multiply by L later)
  helperConstants.vphitheta = 2*ms*ms*ms
              + 48*mb*mb*mq
              + 12*mb*ms*ms
              + 16*mb*mb*ms
              + 8*mq*ms*ms
              + 40*mb*mq*ms; // V_phi_theta, velocity-related coefficient (multiply by L later)
  helperConstants.gblock = g*( 2*ms*ms*ms
            + 24*mb*mq*mq
            + 24*mb*mb*mq
            + 10*mb*ms*ms
            + 8*mb*mb*ms
            + 10*mq*ms*ms
            + 8*mq*mq*ms
            + 40*mb*mq*ms ); // G_block used only in yexp(2)
  helperConstants.expdenom = 2.0f*(mb + mq + ms)*(12.0f*mb*mq + 4.0f*mb*ms + 4.0f*mq*ms + ms*ms);

}

/* (2) Required test function for estimator */

void estimatorOutOfTreeTest() {
  return estimatorKalmanTest(); // cascaded from
}

/* (3) Estimator update function */

void estimatorOutOfTree(state_t *state, const stabilizerStep_t stabilizerStep) {
  estimatorKalman(state, stabilizerStep); // cascaded from

  // Signal pendulum task to update its state
  xSemaphoreGive(runTaskSemaphoreEP);
}

/* High level task functions */

// Called from app main - starts the pendulum task
void estimatorPendulumTaskInit() {
  runTaskSemaphoreEP = xSemaphoreCreateBinary();
  ASSERT(runTaskSemaphoreEP);

  dataMutexEP = xSemaphoreCreateMutexStatic(&dataMutexBufferEP);

  STATIC_MEM_TASK_CREATE(pendulumTask, pendulumTask, KALMAN_TASK_NAME, NULL, PENDULUM_TASK_PRI);

  isInit = true;
}

// For completion sake
bool estimatorPendulumTaskTest() {
  return isInit;
}

/* The Task */

// Excluded rate supervisor, estimator supervisor (check if in bounds), statistics loggers for simplicity
static void pendulumTask(void* parameters) {
  systemWaitStart();

  uint32_t nowMs = T2M(xTaskGetTickCount());
  uint32_t nextPredictionMs = nowMs;

  bool quadIsFlying = supervisorIsFlying();

  while (true) {
    xSemaphoreTake(runTaskSemaphoreEP, portMAX_DELAY);
    nowMs = T2M(xTaskGetTickCount()); 

    // 8bit PWM value (from 16bit)
    uint8_t m1 = (uint8_t) (motorsGetRatio(MOTOR_M1)>>8);
    uint8_t m2 = (uint8_t) (motorsGetRatio(MOTOR_M2)>>8);
    uint8_t m3 = (uint8_t) (motorsGetRatio(MOTOR_M3)>>8);
    uint8_t m4 = (uint8_t) (motorsGetRatio(MOTOR_M4)>>8);
    
    // Convert to Thrust (g, aka gram force) then N
    // https://www.bitcraze.io/documentation/repository/crazyflie-firmware/master/functional-areas/pwm-to-thrust/
    m1 = (0.000409*m1*m1 + 0.1405*m1 + -0.099) * 0.00980665;
    m2 = (0.000409*m2*m2 + 0.1405*m2 + -0.099) * 0.00980665;
    m3 = (0.000409*m3*m3 + 0.1405*m3 + -0.099) * 0.00980665;
    m4 = (0.000409*m4*m4 + 0.1405*m4 + -0.099) * 0.00980665;
    
    float Fl = m1 + m2; // N 
    float Fr = m3 + m4; // N

    // Predict state every OBSFREQ Hz
    if (nowMs >= nextPredictionMs) {
      // Filter forces if needed - moving average filter
      // Fl = Fl*(1-alpha) + Fl_prev*alpha

      /* PREDICT STEP */
      
      // Do not dequeue measurements- it will interfere with measurements that kalman filter uses
      /* 
      measurement_t m;
      estimatorDequeue(&m);
      axis3fSubSamplerAccumulate(&gyroSubSampler, &m.data.gyroscope.gyro);
      gyroLatest = m.data.gyroscope.gyro;
      */ 

      sensorsReadGyro(&gyroLatest); // take snapshot of measurement
      // queue is for collecting last bundle of measurements for integration purposes
      
      /* ADD PROCESS NOISE - in predict function */
      pendulumCorePredict(&pendulumCoreData, &pendulumCoreParams, &gyroLatest, Fl, Fr, nowMs, quadIsFlying);
      //nextPredictionMs = nowMs + PREDICTION_UPDATE_INTERVAL_MS; // moved
      
      // Process noise addition and correction step were moved into the OBSFREQ Hz loop so that 
      // entire observer runs at OBSFREQ Hz. Crazyflie splits by running prediction step at 100 Hz
      // and adding process noise and correcting faster (I think 1000 Hz is loop frequency)

      /* Get measurement */
      sensorsReadAcc(&accLatest);

      /* CORRECT */
      pendulumCoreCorrect(&pendulumCoreData, &pendulumCoreParams, &gyroLatest, &accLatest, Fl, Fr, nowMs, quadIsFlying);
      
      // Write to our storage variable so others (e.g. controller) can read
      xSemaphoreTake(dataMutexEP, portMAX_DELAY);
      pendulumEstimatorState.theta = pendulumCoreData.S[THETA];
      pendulumEstimatorState.theta_dot = pendulumCoreData.S[THETA_DOT];
      pendulumEstimatorState.timestamp = nowMs;
      xSemaphoreGive(dataMutexEP);

      //STATS_CNT_RATE_EVENT(&updateCounter);
      nextPredictionMs = nowMs + PREDICTION_UPDATE_INTERVAL_MS;
    }
  }
}

/* The magic */

// Main prediction step
void pendulumCorePredict(pendulumCoreData_t* this, const pendulumCoreParams_t *params, const Axis3f *gyro, float Fl, float Fr, const uint32_t nowMs, bool quadIsFlying) {
  float dt = (nowMs - this->lastPredictionMs) / 1000.0f; // seconds
  // The linearized update matrices
  NO_DMA_CCM_SAFE_ZERO_INIT static float A[STATE_DIM][STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Am = { STATE_DIM, STATE_DIM, (float *)A}; // linearized dynamics for covariance update;

  NO_DMA_CCM_SAFE_ZERO_INIT static float B[STATE_DIM][INPUT_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Bm = { STATE_DIM, INPUT_DIM, (float *)B};

  // Temp matrices for Q calculation
  NO_DMA_CCM_SAFE_ZERO_INIT static float BT[INPUT_DIM][STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 BTm = { INPUT_DIM, STATE_DIM, (float *)BT};

  NO_DMA_CCM_SAFE_ZERO_INIT static float Q[STATE_DIM][STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Qm = { STATE_DIM, STATE_DIM, (float *)Q};

  // Temporary matrices for the covariance updates
  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpNN1d[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpNN1m = { STATE_DIM, STATE_DIM, tmpNN1d};

  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpNN2d[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpNN2m = { STATE_DIM, STATE_DIM, tmpNN2d};

  //float dt2 = dt*dt;

  // ====== DYNAMICS LINEARIZATION ======
  // For helper constants: 
  // A(1,0) = a1*cosf(theta)*(Fl + Fr) 
  // B(1,0) = b1 + b2*sinf(theta) 
  // B(1,1) = b2*sinf(theta) - b1
  // C(0,0) = c1 * [phi_d^2 + theta_d^2 + 2*phi_d*theta_d]*cosf(phi + theta) + c2 * [(Fl + Fr)*cosf(phi + 2.0*theta)]  
  // C(0,1) = c1 * 2*sinf(phi + theta)*(phi_d + theta_d) 
  // C(1,0) = c1 * [phi_d^2 + theta_d^2 + 2*phi_d*theta_d]*sinf(phi + theta) + c2 * [(Fl + Fr)*sinf(phi + 2.0*theta)]  
  // C(1,1) = -c1 * 2*cosf(phi + theta)*(phi_d + theta_d) 

  // Populate
  float theta = this->S[THETA];
  float theta_d = this->S[THETA_DOT];
  state_t incomingState;
  // The second argument is the current tick, used for prediction/interpolation 
  // if necessary, but usually, passing the current tick is standard.
  estimatorGetState(&incomingState, xTaskGetTickCount());
  float roll_deg  = incomingState.attitude.roll;
  float phi = roll_deg * 3.1415/180;
  float phi_d = gyro->x;

  // Store for correction step
  this->phi_hold = phi;
  this->phi_d_hold = phi_d;

  A[0][0] = 1;
  A[0][1] = dt; // Observer at 100 Hz, not 500 Hz. Don't hardcode anyway
  A[1][0] = dt*helperConstants.a1*cosf(theta)*(Fl + Fr);
  A[1][1] = 1;

  B[0][0] = 0;
  B[0][1] = 0;
  B[1][0] = dt*(helperConstants.b1 + helperConstants.b2*sinf(theta));
  B[1][1] = dt*(helperConstants.b2*sinf(theta) - helperConstants.b1);

  C[0][0] = helperConstants.c1 * (powf(phi_d,2) + powf(theta_d,2) + 2*phi_d*theta_d)*cosf(phi + theta) + helperConstants.c2 * ((Fl + Fr)*cosf(phi + 2.0*theta));
  C[0][1] = helperConstants.c1 * 2*sinf(phi + theta)*(phi_d + theta_d);
  C[1][0] = helperConstants.c1 * (powf(phi_d,2) + powf(theta_d,2) + 2*phi_d*theta_d)*sinf(phi + theta) + helperConstants.c2 * ((Fl + Fr)*sinf(phi + 2.0*theta));
  C[1][1] = -helperConstants.c1 * 2*cosf(phi + theta)*(phi_d + theta_d);

  this->C[0][0] = C[0][0];
  this->C[0][1] = C[0][1];
  this->C[1][0] = C[1][0];
  this->C[1][1] = C[1][1];

  // ====== COVARIANCE UPDATE ======
  mat_mult(&Am, &this->Pm, &tmpNN1m); // A P
  mat_trans(&Am, &tmpNN2m); // A'
  mat_mult(&tmpNN1m, &tmpNN2m, &this->Pm); // A P A'
  
  // Calculate Q = B * B' * 0.0001
  mat_trans(&Bm, &BTm);                 // B'
  mat_mult(&Bm, &BTm, &Qm);             // B * B'
  mat_scale(&Qm, pendulumCoreParams.gamma, &Qm); // Scale by gamma = 0.0001

  // Add process noise
  // P = P + Q
  mat_add(&this->Pm, &Qm, &this->Pm);

  // ====== PREDICTION STEP ======
  float theta_prev = this->S[THETA];
  this->S[THETA] += dt * this->S[THETA_DOT];
  //this->S[THETA_DOT] += dt * -(Fr*L*ms^2*r - Fl*L*ms^2*r + 12*Fl*Ixx*mb*sinf(theta_prev) + 12*Fr*Ixx*mb*sinf(theta_prev) + 6*Fl*Ixx*ms*sinf(theta_prev) + 6*Fr*Ixx*ms*sinf(theta_prev) - 12*Fl*L*mb*mq*r - 4*Fl*L*mb*ms*r + 12*Fr*L*mb*mq*r + 4*Fr*L*mb*ms*r - 4*Fl*L*mq*ms*r + 4*Fr*L*mq*ms*r)/(Ixx*L*(12*mb*mq + 4*mb*ms + 4*mq*ms + ms^2));
  this->S[THETA_DOT] += dt *
  - (Fl - Fr) * (- 12*pendulumCoreParams.Ixx*pendulumCoreParams.mb*sinf(theta_prev) 
  - 6*pendulumCoreParams.Ixx*pendulumCoreParams.ms*sinf(theta_prev) + helperConstants.pred1) 
  / helperConstants.pred2;

  this->isUpdated = true;
  this->lastPredictionMs = nowMs;
}

void pendulumCoreCorrect(pendulumCoreData_t* this, const pendulumCoreParams_t *params, const Axis3f *gyro, const Axis3f *acc, float Fl, float Fr, const uint32_t nowMs, bool quadIsFlying)
{
  // The Kalman gain matrix
  NO_DMA_CCM_SAFE_ZERO_INIT static float L[STATE_DIM][MEAS_DIM];
  static arm_matrix_instance_f32 Lm = {STATE_DIM, MEAS_DIM, (float *)L};

  // The linearized update matrices
  NO_DMA_CCM_SAFE_ZERO_INIT static float C[MEAS_DIM][STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Cm = { MEAS_DIM, STATE_DIM, (float *)C};

  // Grab local C from data
  memcpy(C, this->C, sizeof(C));

  // Temporary matrices for math
  
  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpNN1d[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpNN1m = { STATE_DIM, STATE_DIM, tmpNN1d};

  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpNN2d[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpNN2m = { STATE_DIM, STATE_DIM, tmpNN2d};

  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpNN3d[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpNN3m = { STATE_DIM, STATE_DIM, tmpNN3d};

  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpNN4d[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpNN4m = { STATE_DIM, STATE_DIM, tmpNN4d};

  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpNN5d[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpNN5m = { STATE_DIM, STATE_DIM, tmpNN5d};

  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpNN6d[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpNN6m = { STATE_DIM, STATE_DIM, tmpNN6d};

  NO_DMA_CCM_SAFE_ZERO_INIT static float tmpNN7d[STATE_DIM * STATE_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 tmpNN7m = { STATE_DIM, STATE_DIM, tmpNN7d};

  // R covariance matrix
  NO_DMA_CCM_SAFE_ZERO_INIT static float R[MEAS_DIM][MEAS_DIM];
  static __attribute__((aligned(4))) arm_matrix_instance_f32 Rm = { MEAS_DIM, MEAS_DIM, (float *)R};

  R[0][0] = 1 * pendulumCoreParams.beta;
  R[0][1] = 0;
  R[1][0] = 0; 
  R[1][1] = 1 * pendulumCoreParams.beta;

  // ====== MEASUREMENT UPDATE ======

  /*
  % Compute Kalman gain (C-1)
    Lgain = (P*C')/(R + C*P*C');
    % Update (C-3)
    P = P - Lgain*C*P; 

    % Correction (C-2)
    %yexp = C*xbar + D*u + dmeas_const % replaced this- use nonlinear model
    yexp(1) = (12*Fl*mb^2*sinf(phi + 2*theta) - 12*Fr*mb^2*sinf(phi) - 5*Fl*ms^2*sinf(phi) - 5*Fr*ms^2*sinf(phi) - 12*Fl*mb^2*sinf(phi) + 12*Fr*mb^2*sinf(phi + 2*theta) + 3*Fl*ms^2*sinf(phi + 2*theta) + 3*Fr*ms^2*sinf(phi + 2*theta) + L*ms^3*theta_d^2*sinf(phi + theta) - 24*Fl*mb*mq*sinf(phi) - 20*Fl*mb*ms*sinf(phi) - 24*Fr*mb*mq*sinf(phi) - 20*Fr*mb*ms*sinf(phi) - 8*Fl*mq*ms*sinf(phi) - 8*Fr*mq*ms*sinf(phi) + 12*Fl*mb*ms*sinf(phi + 2*theta) + 12*Fr*mb*ms*sinf(phi + 2*theta) + L*ms^3*phi_d^2*sinf(phi + theta) + 2*L*ms^3*phi_d*theta_d*sinf(phi + theta) + 24*L*mb^2*mq*phi_d^2*sinf(phi + theta) + 6*L*mb*ms^2*phi_d^2*sinf(phi + theta) + 8*L*mb^2*ms*phi_d^2*sinf(phi + theta) + 4*L*mq*ms^2*phi_d^2*sinf(phi + theta) + 24*L*mb^2*mq*theta_d^2*sinf(phi + theta) + 6*L*mb*ms^2*theta_d^2*sinf(phi + theta) + 8*L*mb^2*ms*theta_d^2*sinf(phi + theta) + 4*L*mq*ms^2*theta_d^2*sinf(phi + theta) + 20*L*mb*mq*ms*phi_d^2*sinf(phi + theta) + 20*L*mb*mq*ms*theta_d^2*sinf(phi + theta) + 48*L*mb^2*mq*phi_d*theta_d*sinf(phi + theta) + 12*L*mb*ms^2*phi_d*theta_d*sinf(phi + theta) + 16*L*mb^2*ms*phi_d*theta_d*sinf(phi + theta) + 8*L*mq*ms^2*phi_d*theta_d*sinf(phi + theta) + 40*L*mb*mq*ms*phi_d*theta_d*sinf(phi + theta))/(2*(mb + mq + ms)*(12*mb*mq + 4*mb*ms + 4*mq*ms + ms^2));
    yexp(2) = -(2*g*ms^3 + 24*g*mb*mq^2 + 24*g*mb^2*mq + 10*g*mb*ms^2 + 8*g*mb^2*ms + 10*g*mq*ms^2 + 8*g*mq^2*ms - 12*Fl*mb^2*cosf(phi) - 12*Fr*mb^2*cosf(phi) - 5*Fl*ms^2*cosf(phi) - 5*Fr*ms^2*cosf(phi) + 12*Fl*mb^2*cosf(phi + 2*theta) + 12*Fr*mb^2*cosf(phi + 2*theta) + 3*Fl*ms^2*cosf(phi + 2*theta) + 3*Fr*ms^2*cosf(phi + 2*theta) + 40*g*mb*mq*ms - 24*Fl*mb*mq*cosf(phi) - 20*Fl*mb*ms*cosf(phi) - 24*Fr*mb*mq*cosf(phi) - 20*Fr*mb*ms*cosf(phi) - 8*Fl*mq*ms*cosf(phi) - 8*Fr*mq*ms*cosf(phi) + 12*Fl*mb*ms*cosf(phi + 2*theta) + 12*Fr*mb*ms*cosf(phi + 2*theta) + L*ms^3*phi_d^2*cosf(phi + theta) + L*ms^3*theta_d^2*cosf(phi + theta) + 2*L*ms^3*phi_d*theta_d*cosf(phi + theta) + 24*L*mb^2*mq*phi_d^2*cosf(phi + theta) + 6*L*mb*ms^2*phi_d^2*cosf(phi + theta) + 8*L*mb^2*ms*phi_d^2*cosf(phi + theta) + 4*L*mq*ms^2*phi_d^2*cosf(phi + theta) + 24*L*mb^2*mq*theta_d^2*cosf(phi + theta) + 6*L*mb*ms^2*theta_d^2*cosf(phi + theta) + 8*L*mb^2*ms*theta_d^2*cosf(phi + theta) + 4*L*mq*ms^2*theta_d^2*cosf(phi + theta) + 20*L*mb*mq*ms*phi_d^2*cosf(phi + theta) + 20*L*mb*mq*ms*theta_d^2*cosf(phi + theta) + 48*L*mb^2*mq*phi_d*theta_d*cosf(phi + theta) + 12*L*mb*ms^2*phi_d*theta_d*cosf(phi + theta) + 16*L*mb^2*ms*phi_d*theta_d*cosf(phi + theta) + 8*L*mq*ms^2*phi_d*theta_d*cosf(phi + theta) + 40*L*mb*mq*ms*phi_d*theta_d*cosf(phi + theta))/(2*(mb + mq + ms)*(12*mb*mq + 4*mb*ms + 4*mq*ms + ms^2));
    xhat = xbar + Lgain*(y_meas - yexp); 
  */

  // Calculate the Kalman gain and perform the state update
  mat_trans(&Cm, &tmpNN1m); // C'
  mat_mult(&this->Pm, &tmpNN1m, &tmpNN2m); // P C'
  mat_mult(&Cm, &tmpNN2m, &tmpNN3m); // C P C'
  mat_add(&Rm, &tmpNN3m, &tmpNN4m); // R + C P C'
  mat_inv(&tmpNN4m, &tmpNN5m); // (R + C P C')^-1
  mat_mult(&tmpNN2m, &tmpNN5m, &Lm); // 

  // ====== COVARIANCE UPDATE ======
  mat_mult(&Cm, &this->Pm, &tmpNN6m); // C P
  mat_mult(&Lm, &tmpNN6m, &tmpNN7m); // L C P
  mat_sub(&this->Pm, &tmpNN7m, &this->Pm); // P - L C P
  
  float theta = this->S[THETA];
  float theta_d = this->S[THETA_DOT];
  float phi = this->phi_hold;
  float phi_d = this->phi_d_hold;

  // Calculate expected measurement         
  float numerator_yexp1 =
    (Fl+Fr) * helperConstants.cphi2 * sinf(phi+2*theta)
  - (Fl+Fr) * helperConstants.cphi  * sinf(phi)
  + pendulumCoreParams.L * sinf(phi+theta) * ( helperConstants.vphi2*(powf(phi_d,2) + powf(theta_d,2)) + helperConstants.vphitheta*(phi_d*theta_d) );

  float numerator_yexp2 =
    -helperConstants.gblock
  + (Fl + Fr) * helperConstants.cphi * cosf(phi)
  - (Fl + Fr) * helperConstants.cphi2 * cosf(phi + 2.0f*theta)
  - pendulumCoreParams.L * cosf(phi + theta) * ( helperConstants.vphi2 * (powf(phi_d,2) + powf(theta_d,2)) + helperConstants.vphitheta * (phi_d*theta_d) );

  float yexp1 = numerator_yexp1 / helperConstants.expdenom;
  float yexp2 = numerator_yexp2 / helperConstants.expdenom;
  
  // Update state
  // xhat = xbar + Lgain*(y_meas - yexp); 
  float ymeas1 = acc->y;
  float ymeas2 = acc->z;

  this->S[THETA] += L[0][0]*(ymeas1 - yexp1) + L[0][1]*(ymeas2 - yexp2);
  this->S[THETA_DOT] += L[1][0]*(ymeas1 - yexp1) + L[1][1]*(ymeas2 - yexp2);

  this->isUpdated = true;
}

void pendulumEstimatorGetState(pendulum_state_t *out)
{
  // mutex makes sure others can't read while I write  
  xSemaphoreTake(dataMutexEP, portMAX_DELAY);
  memcpy(out, &pendulumEstimatorState, sizeof(pendulum_state_t));
  xSemaphoreGive(dataMutexEP);
}
/* 
Read with the getter above like this:

pendulum_state_t copy;
pendulumEstimatorGetState(&copy);
float theta = copy.theta;
*/

LOG_GROUP_START(pendEKF)

  /* ============
   *   STATES
   * ============ */

  /** @brief State: theta (angle) */
  LOG_ADD(LOG_FLOAT, theta, &pendulumCoreData.S[THETA])

  /** @brief State: theta_dot (angular velocity) */
  LOG_ADD(LOG_FLOAT, thetaDot, &pendulumCoreData.S[THETA_DOT])


  /* =====================
   *   STATE COVARIANCES
   * ===================== */

  /** @brief Var(theta) */
  LOG_ADD(LOG_FLOAT, varTheta, &pendulumCoreData.P[THETA][THETA])

  /** @brief Var(theta_dot) */
  LOG_ADD(LOG_FLOAT, varThetaDot, &pendulumCoreData.P[THETA_DOT][THETA_DOT])

LOG_GROUP_STOP(pendEKF)
