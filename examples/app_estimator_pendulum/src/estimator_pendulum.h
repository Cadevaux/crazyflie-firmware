/**
 * Estimator pendulum parameters and structs
 */

#pragma once

#include <stdint.h>
#include "stabilizer_types.h"
#include "cf_math.h" // for matrix math, includes arm_math.h

#define PENDULUM_TASK_NAME         "PENDULUM"
#define PENDULUM_TASK_STACKSIZE    (3 * configMINIMAL_STACK_SIZE)
#define PENDULUM_TASK_PRI          2

#define PENDULUM_DEFAULT_PARAMS_INIT \
  /* Initial params, states, and variances */ \
  .mq = 0.0314f,     /* kg */ \
  .mb = 0.00633f,    /* kg */ \
  .ms = 0.0001f,     /* kg */ \
  .Ixx = 6.410179e-06f,  /* kg m2*/ \
  .g = 9.81f,        /* m/s2 */ \
  .r = 0.028325f,   /* m */ \
  .L = 0.4f,         /* m */ \
  \
  .initialTheta = 3.14159265f,     /* radians */ \
  .initialThetaDot = 0.0f,  /* radians per second */ \
  \
  .stdDevInitialTheta = 0.01f,      /* radians (5 deg) */ \
  .stdDevInitialThetaDot = 0.01f,    /* radians per second */ \
  \
  .gamma = 0.00001f, /* Q tuning */ \
  .beta = 0.05f,  /* R tuning*/ \

// Struct to externalize state
typedef struct {
    float theta;
    float theta_dot;
    uint32_t timestamp;
} pendulum_state_t;

// Accessor
void pendulumEstimatorGetState(pendulum_state_t *out);

typedef struct {
  float mq;
  float mb;
  float ms;
  float Ixx;
  float g;
  float r;
  float L;

  float initialTheta;
  float initialThetaDot;

  float stdDevInitialTheta;
  float stdDevInitialThetaDot;

  float gamma;
  float beta;

} pendulumCoreParams_t;

typedef struct {
  float a1;
  float b1;
  float b2;
  float c1;
  float c2;
  float pred1;
  float pred2;
  float cphi2;
  float cphi;
  float vphi2;
  float vphitheta;
  float gblock;
  float expdenom;
} helperConstants_t;

typedef enum
{
  THETA, THETA_DOT, STATE_DIM
} pendulumCoreStateIdx_t;

#define INPUT_DIM 2
#define MEAS_DIM 2

typedef struct {
  // Pendulum state
  float S[STATE_DIM];

  // The covariance matrix
  __attribute__((aligned(4))) float P[STATE_DIM][STATE_DIM];
  arm_matrix_instance_f32 Pm;

  // Linearized C matrix used in Correct, but linearized with everything in Predict
  __attribute__((aligned(4))) float C[MEAS_DIM][STATE_DIM];
  arm_matrix_instance_f32 Cm;

  // Tracks whether an update to the state has been made, and the state therefore requires finalization
  bool isUpdated;

  // Variables to hold over from prediction to update
  float phi_hold;
  float phi_d_hold;

  uint32_t lastPredictionMs;
  uint32_t lastProcessNoiseUpdateMs;
} pendulumCoreData_t;

void estimatorPendulumTaskInit();
bool estimatorPendulumTaskTest();