
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

#ifndef __PIC_IMPL_H
#define __PIC_IMPL_H

#define Const_Boltzmann		1.38E-23		//	Boltzmann constant  [J/K]
#define Const_Permitt_0		8.854E-12		//	vacuum permittivity
#define Const_e				1.6E-19			//	elementary charge [C]
#define	Const_PI			3.14159

#define Const_AMU_Mass		1.66E-27		//	mass conversion of atomic mass unit in [kg]
#define Const_ERM_Mass		9.109E-31		//	electron rest mass in [kg]

extern double randf();
extern double Var_NormalDistri( double Mean, double Sigma );
extern double Velocity_MaxwellianDistri( double mass, double T, double mu );
extern void TDMA_solve( double* L, double* D, double* U, double* rhs, double* x );

#endif

