/*
 * MATLAB Compiler: 4.16 (R2011b)
 * Date: Sat Jul 07 17:13:54 2012
 * Arguments: "-B" "macro_default" "-W" "lib:libmatdyn_1_2" "-T" "link:lib"
 * "-d" "U:\work\MATLAB\c_dll\matdyn_1_2\libmatdyn_1_2\src" "-w"
 * "enable:specified_file_mismatch" "-w" "enable:repeated_file" "-w"
 * "enable:switch_ignored" "-w" "enable:missing_lib_sentinel" "-w"
 * "enable:demo_license" "-v" "U:\work\MATLAB\c_dll\matdyn_1_2\loadflow.m"
 * "U:\work\MATLAB\c_dll\matdyn_1_2\mdisplay.m"
 * "U:\work\MATLAB\c_dll\matdyn_1_2\Simulate.m" 
 */

#ifndef __libmatdyn_1_2_h
#define __libmatdyn_1_2_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libmatdyn_1_2
#define PUBLIC_libmatdyn_1_2_C_API __global
#else
#define PUBLIC_libmatdyn_1_2_C_API /* No import statement needed. */
#endif

#define LIB_libmatdyn_1_2_C_API PUBLIC_libmatdyn_1_2_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libmatdyn_1_2
#define PUBLIC_libmatdyn_1_2_C_API __declspec(dllexport)
#else
#define PUBLIC_libmatdyn_1_2_C_API __declspec(dllimport)
#endif

#define LIB_libmatdyn_1_2_C_API PUBLIC_libmatdyn_1_2_C_API


#else

#define LIB_libmatdyn_1_2_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libmatdyn_1_2_C_API 
#define LIB_libmatdyn_1_2_C_API /* No special import/export declaration */
#endif

extern LIB_libmatdyn_1_2_C_API 
bool MW_CALL_CONV libmatdyn_1_2InitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_libmatdyn_1_2_C_API 
bool MW_CALL_CONV libmatdyn_1_2Initialize(void);

extern LIB_libmatdyn_1_2_C_API 
void MW_CALL_CONV libmatdyn_1_2Terminate(void);



extern LIB_libmatdyn_1_2_C_API 
void MW_CALL_CONV libmatdyn_1_2PrintStackTrace(void);

extern LIB_libmatdyn_1_2_C_API 
bool MW_CALL_CONV mlxLoadflow(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

extern LIB_libmatdyn_1_2_C_API 
bool MW_CALL_CONV mlxMdisplay(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

extern LIB_libmatdyn_1_2_C_API 
bool MW_CALL_CONV mlxSimulate(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

extern LIB_libmatdyn_1_2_C_API 
long MW_CALL_CONV libmatdyn_1_2GetMcrID();



extern LIB_libmatdyn_1_2_C_API bool MW_CALL_CONV mlfLoadflow(int nargout, mxArray** baseMVA, mxArray** bus, mxArray** gen, mxArray** branch, mxArray** success, mxArray* baseMVAIn, mxArray* busIn, mxArray* genIn, mxArray* branchIn, mxArray* areasIn, mxArray* gencostIn, mxArray* qlim, mxArray* dc, mxArray* alg, mxArray* tol, mxArray* max_it);

extern LIB_libmatdyn_1_2_C_API bool MW_CALL_CONV mlfMdisplay(mxArray* in);

extern LIB_libmatdyn_1_2_C_API bool MW_CALL_CONV mlfSimulate(int nargout, mxArray** angles, mxArray** speeds, mxArray** eq_tr, mxArray** ed_tr, mxArray** efd, mxArray** PM, mxArray** voltages, mxArray** stepsize, mxArray** errest, mxArray** time, mxArray** success, mxArray* baseMVAIn, mxArray* busIn, mxArray* genIn, mxArray* branchIn, mxArray* areasIn, mxArray* genCostIn, mxArray* SSqlim, mxArray* SSdc, mxArray* SSalg, mxArray* SStol, mxArray* SSmax_it, mxArray* genDyn, mxArray* excDyn, mxArray* govDyn, mxArray* TDfreq, mxArray* TDstep, mxArray* TDstoptime, mxArray* TDmethod, mxArray* TDtol, mxArray* TDstepMin, mxArray* TDstepMax, mxArray* event, mxArray* buschange, mxArray* linechange);

#ifdef __cplusplus
}
#endif
#endif
