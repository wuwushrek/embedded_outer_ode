//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ode_dyn.h
//
// MATLAB Coder version            : 3.4
// C/C++ source code generated on  : 06-Aug-2018 12:52:40
//
#ifndef ODE_DYN_H
#define ODE_DYN_H

// Include Files
#include "aa_mod.h"

// Function Declarations
extern void ode_derivatives(const AF1 in1[2], AF1 ode_der[4]);
extern void ode_function(const AF1 in1[2], AF1 ode_fun[2]);
extern void ode_reminder(const AF1 in1[2], AF1 ode_rem[2]);

#endif

//
// File trailer for ode_dyn.h
//
// [EOF]
//
