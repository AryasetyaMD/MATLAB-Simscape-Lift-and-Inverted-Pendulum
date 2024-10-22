/* ==SOURCE file list of “PLL.c” with Structure B == */
#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME PLL

#include "simstruc.h"
#include <math.h>

// Constant
real_T pi = 3.141592653589793;
real_T factor1 = 0.8164965809;
real_T factor2 = 0.8660254038;

// IM Parameter
real_T Np = 2.0;
real_T Lm = 0.235,  Lr = 0.2473, Ls = 0.2442;
real_T Rs = 3.67, Rr = 2.32;

#define U(element) (*uPtrs[element]) /* Pointer to Input Port0 */

static void mdlInitializeSizes(SimStruct *S) {
    ssSetNumDiscStates(S, 8);
    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, 10);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortOverWritable(S, 0, 1);
    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, 4);
    ssSetNumSampleTimes(S, 1);
    ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE | SS_OPTION_DISCRETE_VALUED_OUTPUT));
}

static void mdlInitializeSampleTimes(SimStruct *S) {
    ssSetSampleTime(S, 0, 1e-4);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S) {
    real_T *X = ssGetRealDiscStates(S);
    int_T nXStates = ssGetNumDiscStates(S);
    int_T i;

    /* initialize the states to 0.0 */
    for (i = 0; i < nXStates; i++) {
        X[i] = 0.0;
    }
}

static void mdlOutputs(SimStruct *S, int_T tid) {
    real_T *Y = ssGetOutputPortRealSignal(S, 0);
    real_T *X = ssGetRealDiscStates(S);

    // Output
    real_T theta_e_est = X[2];
    real_T theta_r_est = X[3];
    real_T omega_e_est_filtered = X[6];
    real_T omega_r_est_filtered = X[7];

    Y[0] = omega_e_est_filtered;
    Y[1] = omega_r_est_filtered;
    Y[2] = theta_e_est;
    Y[3] = theta_r_est;
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid) {
    real_T *X = ssGetRealDiscStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);

    real_T dt = 1e-4;

    real_T Va = U(0);
    real_T Vb = U(1);
    real_T Vc = U(2);
    real_T theta_e = U(3);
    real_T omega_ref = U(4);
    // Type-3 PLL Gain
    real_T k1 = U(5);
    real_T k2 = U(6);
    real_T k3 = U(7);
    real_T fp = U(8);
    real_T alpha_omega = U(9);

    // State
    real_T integral_pll_old = X[0];
    real_T double_integral_pll_old = X[1];
    real_T theta_e_est_old = X[2];
    real_T theta_r_est_old = X[3];
    real_T Vsd_snsrl_filtered_old = X[4];
    real_T Vsq_snsrl_filtered_old = X[5];
    real_T omega_est_tuned_old = X[6];

    // Voltage ABC to Alpha-Beta
    real_T V_snsrl_alpha = factor1 * (Va - 0.5 * Vb - 0.5 * Vc);
    real_T V_snsrl_beta = factor1 * (factor2 * Vb - factor2 * Vc);

    // Sensorless Type-3 PLL
    real_T Vsd_snsrl = (V_snsrl_alpha * cos(theta_e)) + (V_snsrl_beta * sin(theta_e));
    real_T Vsq_snsrl = (-V_snsrl_alpha * sin(theta_e)) + (V_snsrl_beta * cos(theta_e));

    // Pre-filter (low-pass filter)
    real_T wp = 2 * pi * fp; // cutoff frequency (rad/s)
    real_T alpha = wp * dt / (1 + wp * dt);

    real_T Vsd_snsrl_filtered = alpha * Vsd_snsrl + (1 - alpha) * Vsd_snsrl_filtered_old;
    real_T Vsq_snsrl_filtered = alpha * Vsq_snsrl + (1 - alpha) * Vsq_snsrl_filtered_old;

    // Update filter state
    X[4] = Vsd_snsrl_filtered;
    X[5] = Vsq_snsrl_filtered;

    // Amplitude Normalization (AN)
    real_T amplitude = 1;
    real_T Vsd_snsrl_normalized = Vsd_snsrl_filtered / amplitude;
    real_T Vsq_snsrl_normalized = Vsq_snsrl_filtered / amplitude;

    real_T error = Np * omega_ref - Vsq_snsrl_normalized;
    real_T integral_pll = integral_pll_old + error * dt;
    real_T double_integral_pll = double_integral_pll_old + integral_pll_old * dt;
    real_T omega_est = k1 * error + k2 * integral_pll + k3 * double_integral_pll;

    // omega_est smoother (slower response)
    real_T omega_est_tuned = alpha_omega * omega_est + (1 - alpha_omega) * omega_est_tuned_old;

    real_T theta_e_est = theta_e_est_old + omega_est_tuned * dt;
    real_T theta_r_est = theta_r_est_old + (omega_est_tuned / Np) * dt;

    while (theta_e_est > 2 * pi) {
        theta_e_est -= 2 * pi;
    }

    while (theta_e_est < -2 * pi) {
        theta_e_est += 2 * pi;
    }

    while (theta_r_est > 2 * pi) {
        theta_r_est -= 2 * pi;
    }

    while (theta_r_est < -2 * pi) {
        theta_r_est += 2 * pi;
    }

    X[0] = integral_pll;
    X[1] = double_integral_pll;
    X[2] = theta_e_est;
    X[3] = theta_r_est;
    X[6] = omega_est_tuned; // omega_e_est
    X[7] = (omega_est_tuned / Np); // omega_r_est
}

static void mdlTerminate(SimStruct *S) {
    /* Keep this function empty since no memory is allocated */
}

#ifdef MATLAB_MEX_FILE
/* Is this file being compiled as a MEX-file? */
#include "simulink.c" /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /* Code generation registration function */
#endif
