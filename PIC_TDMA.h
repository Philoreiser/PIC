#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "PIC_Impl.h"

using namespace std;

#ifndef __PIC_TDMA_H
#define __PIC_TDMA_H

class PIC_TDMA {

	public:

		int				rank;

		double*			L;
		double*			D;
		double*			U;
		double*			rhs;

		double*			X;

		PIC_TDMA();
		PIC_TDMA(int);
		~PIC_TDMA();

		void tdma_init(int);
		void tdma_mat_set_default();
		void tdma_rhs_set_default( double, double* , int );
		void tdma_solve( double* x );
		void tdma_solve_default();

		void tdma_verification();
};

#endif

