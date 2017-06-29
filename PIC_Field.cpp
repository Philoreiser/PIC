#include "PIC_Field.h"

PIC_Field::PIC_Field() { };
PIC_Field::~PIC_Field() { };

void PIC_Field::alloc( int num_of_nodes ) {
	
	NN = num_of_nodes;

	phi = new double [num_of_nodes];
	E = new double [num_of_nodes];
	rho = new double [num_of_nodes];

	for ( int i = 0; i < NN; i++ )
	{
		phi[i] = .0;
		E[i] = .0;
		rho[i] = .0;
	}

	PoissonTDMA = new PIC_TDMA( num_of_nodes - 2 ); // for two-side grounded case

};

void PIC_Field::solve_phi_by_TDMA_default() {

	int		i_sol = 0;

	PoissonTDMA->tdma_mat_set_default();
	PoissonTDMA->tdma_rhs_set_default( delta_x, rho, 1 );
	PoissonTDMA->tdma_solve_default();

	phi[0] = .0;
	phi[NN-1] = .0;

	for ( int i = 1; i < NN-1; i++ )
	{
		i_sol = i-1;
		phi[i] = PoissonTDMA->X[i_sol];
	}

};

void PIC_Field::resolve_field_default() {

	E[0] = -( phi[1] - phi[0] ) / delta_x;
	E[NN-1] = -( phi[NN-1] - phi[NN-2] ) / delta_x;
	
	for ( int i = 1; i < NN-1; i++ ) E[i] = -( phi[i+1] - phi[i-1] ) / delta_x / 2.0;

};

void PIC_Field::set_rho_zero() {

	for ( int i = 0; i < NN; i++ ) rho[i] = .0;
};

