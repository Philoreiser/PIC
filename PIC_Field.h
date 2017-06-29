#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "PIC_Impl.h"
#include "PIC_TDMA.h"

using namespace std;

#if !defined(__PIC_FIELD_H)
#define __PIC_FIELD_H

class PIC_Field {

	protected:

	public:

		int			NN;		//	number of nodes
		double*		phi;	//	electrostatic potential
		double*		E;		//	electric field
		double*		rho;	//	volumetric charge density

		double		delta_x;

		PIC_TDMA*	PoissonTDMA;

		PIC_Field();
		~PIC_Field();

		void alloc( int );
		void solve_phi_by_TDMA_default();
		void resolve_field_default();

		void set_rho_zero();
};

#endif

