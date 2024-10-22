/* == SOURCE file list of “RFOC.c” with Structure B == */ 

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME RFOC_Sensor

#include "simstruc.h"
#include <math.h>

// Constant
real_T pi = 22 / 7;
real_T factor1 = 0.8164965809;
real_T factor2 = 0.8660254038;

// IM Parameter
real_T Np = 2.0;
real_T Lm = 0.235,  Lr = 0.2473, Ls = 0.2442;
real_T Rs = 3.67, Rr = 2.32;

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/

static void mdlInitializeSizes(SimStruct *S) {
    ssSetNumDiscStates(S, 14);
    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, 7);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortOverWritable(S, 0, 1);
    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, 10);
    ssSetNumSampleTimes(S, 1);
    ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE | SS_OPTION_DISCRETE_VALUED_OUTPUT));
}

static void mdlInitializeSampleTimes(SimStruct *S) {
    ssSetSampleTime(S, 0, 1e-4);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S) {
    real_T *X0 = ssGetRealDiscStates(S);
    int_T nXStates = ssGetNumDiscStates(S);
    int_T i;

    /* initialize the states to 0.0 */
    for (i = 0; i < nXStates; i++) {
        X0[i] = 0.0;
    }
}

static void mdlOutputs(SimStruct *S, int_T tid) {
    real_T *Y = ssGetOutputPortRealSignal(S, 0);
    real_T *X = ssGetRealDiscStates(S);
    // Output
    real_T imr = X[2];
    real_T theta_e = X[3];
    real_T Te = X[6];
    real_T isd_ref = X[7];
    real_T isd_act = X[8];
    real_T isq_ref = X[9];
    real_T isq_act = X[10];
    real_T Va = X[11];
    real_T Vb = X[12];
    real_T Vc = X[13];
    
    Y[0] = Va;
    Y[1] = Vb;
    Y[2] = Vc;
    Y[3] = imr;
    Y[4] = Te;
    Y[5] = theta_e;
    Y[6] = isd_ref;
    Y[7] = isd_act;
    Y[8] = isq_ref;
    Y[9] = isq_act;
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid) {
    real_T *X = ssGetRealDiscStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);

    real_T dt = 1e-4;
    real_T Td = 0.001;
    // K Formula
    real_T sigma = 1 - ((Lm * Lm) / (Ls * Lr));
    real_T Kidp = sigma * Ls / Td ;
    real_T Kidi = Rs / Td;
    real_T Kiqp = sigma * Ls / Td;
    real_T Kiqi = Rs / Td;

    // Input
    real_T ia = U(0);
    real_T ib = U(1);
    real_T ic = U(2);
    real_T isd_ref = U(3);
    real_T isq_ref = U(4);
    real_T Vs_max = U(5);
    real_T omega_r = U(6);

    // State
    real_T isd_ref_old = X[0];
    real_T isq_ref_old = X[1];
    real_T imr_old = X[2];
    real_T theta_e_old = X[3];
    real_T integral_isd_old = X[4];
    real_T integral_isq_old = X[5];

    // Decoupling
    real_T isd_ref_new = isd_ref_old + ((isd_ref - isd_ref_old) * dt) / Td;
    real_T isq_ref_new = isq_ref_old + ((isq_ref - isq_ref_old) * dt) / Td;

    // Flux Model
    real_T omega_e = Np * omega_r + (Rr / Lr) * (isq_ref_new / isd_ref_new); // check isd_ref_new
    real_T theta_e = theta_e_old + omega_e * dt;

    // Three-Phase ABC to Alpha-Beta Current
    real_T i_alpha = factor1 * (ia - 0.5 * ib - 0.5 * ic);
    real_T i_beta = factor1 * (factor2 * ib - factor2 * ic);

    // Alpha-Beta Current to DQ Current
    real_T isd_act = (i_alpha * cos(theta_e)) + (i_beta * sin(theta_e));
    real_T isq_act = (-i_alpha * sin(theta_e)) + (i_beta * cos(theta_e));

    // Torque
    real_T Te = Np * (1 - sigma) * Ls * imr_old * isq_act;

    real_T imr_new = imr_old + ((Rr / Lr) * (isd_act - imr_old) * dt);

    real_T integral_isd_new = integral_isd_old + (isd_ref - isd_act) * dt;
    real_T integral_isq_new = integral_isq_old + (isq_ref - isq_act) * dt;

    real_T usd = (Kidp * (isd_ref - isd_act)) + (Kidi * integral_isd_new);
    real_T usq = (Kiqp * (isq_ref - isq_act)) + (Kiqi * integral_isq_new);

    real_T Vcd = -(omega_e * Ls * sigma * isq_act) + Ls * (1 - sigma) * imr_new;
    real_T Vcq = omega_e * Ls * sigma * isd_act + Ls * (1 - sigma) * omega_e * imr_old;

    real_T Vsd = usd + Vcd;
    real_T Vsq = usq + Vcq;

    // DQ Voltage to Alpha-Beta Voltage
    real_T V_alpha = (Vsd * cos(theta_e) - Vsq * sin(theta_e));
    real_T V_beta = (Vsd * sin(theta_e) + Vsq * cos(theta_e));

    // Alpha-beta Voltage to Three-Phase ABC Voltage
    real_T Va = factor1 * V_alpha;
    real_T Vb = factor1 * (-0.5 * V_alpha + factor2 * V_beta);
    real_T Vc = factor1 * (-0.5 * V_alpha - factor2 * V_beta);

    // Limit Voltage
    if(Va > Vs_max) {
        Va = Vs_max; 
    }
    else if(Va < (-Vs_max)) { 
        Va = -Vs_max; 
    }
    else if(Vb > (Vs_max)) { 
        Vb = Vs_max; 
    }
    else if(Vb < (-Vs_max)) { 
        Vb = -Vs_max; 
    }
    else if(Vc > Vs_max) { 
        Vc = Vs_max; 
    }
    else if(Vc < (-Vs_max)) { 
        Vc = -Vs_max; 
    }

    // Update State
    X[0] = isd_ref_new;
    X[1] = isq_ref_new;
    X[2] = imr_new;
    X[3] = theta_e;
    X[4] = integral_isd_new;
    X[5] = integral_isq_new;
    X[6] = Te;
    X[7] = isd_ref;
    X[8] = isd_act;
    X[9] = isq_ref;
    X[10] = isq_act;
    X[11] = Va;
    X[12] = Vb;
    X[13] = Vc;
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
