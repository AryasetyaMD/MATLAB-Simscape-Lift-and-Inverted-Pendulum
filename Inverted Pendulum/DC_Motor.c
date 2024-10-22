#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME DC_Motor
#include "simstruc.h"
#include <math.h>

#define U(element) (*uPtrs[element])

static void mdlInitializeSizes(SimStruct *S) {
    ssSetNumContStates(S, 3);
    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, 2);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortOverWritable(S, 0, 1);
    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetNumSampleTimes(S, 1);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct *S) {
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS

static void mdlInitializeConditions(SimStruct *S) {
    real_T *X0 = ssGetContStates(S);
    int_T nStates = ssGetNumContStates(S);
    int_T i;
    for (i = 0; i < nStates; i++) {
        X0[i] = 0.0;
    }
}

static void mdlOutputs(SimStruct *S, int_T tid) {
    real_T *Y = ssGetOutputPortRealSignal(S, 0);
    real_T *X = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);
    Y[0] = X[2] * 0.027;
}

#define MDL_DERIVATIVES

static void mdlDerivatives(SimStruct *S) {
    real_T *dX = ssGetdX(S);
    real_T *X = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);

    real_T L = 0.032;
    real_T R = 1;
    real_T J = 0.0465;
    real_T b = 0.005;
    real_T K_T = 0.027;
    real_T K_e = 0.027;

    real_T V = U(0);
    real_T T_L = U(1);

    real_T theta = X[0];
    real_T omega = X[1];
    real_T I = X[2];

    real_T theta_dot = omega;
    real_T omega_dot = ((K_T / J) * I) - ((b / J) * omega) - ((1 / J) * T_L);
    real_T I_dot = -((R / L) * I) - ((K_e / L) * omega) + ((1 / L) * V);

    dX[0] = theta_dot;
    dX[1] = omega_dot;
    dX[2] = I_dot;
}

static void mdlTerminate(SimStruct *S) {}

#ifdef MATLAB_MEX_FILE
#include "simulink.c" /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /*Code generation registration function*/
#endif