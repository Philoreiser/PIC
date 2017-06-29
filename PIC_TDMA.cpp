#include "PIC_TDMA.h"

PIC_TDMA::PIC_TDMA() { };

PIC_TDMA::PIC_TDMA(int n) {

	tdma_init(n);
};

PIC_TDMA::~PIC_TDMA() { };

void PIC_TDMA::tdma_init(int n) {

	rank = n;
	L = new double [n];
	D = new double [n];
	U = new double [n];
	rhs = new double [n];

	X = new double [n];

	tdma_mat_set_default();

	for ( int i = 0; i < rank; i++ ) 
	{
		rhs[i] = .0;
		X[i] = .0;
	}

};

void PIC_TDMA::tdma_mat_set_default() {
	
	for ( int i = 0; i < rank; i++ )
	{
		L[i] = -1.0;
		D[i] = 2.0;
		U[i] = -1.0;
	}

};

void PIC_TDMA::tdma_rhs_set_default( double dx, double* rho, int i_start ) {

	double		eps = Const_Permitt_0;
	double		dim_factor = dx*dx/eps;

	for ( int i = 0; i < rank; i++ ) 
	{
		rhs[i] = dim_factor * rho[ i + i_start ]; // form right-hand-side
		X[i] = .0; // reset the solution vector as zero before solution algorithm
	}
};


void PIC_TDMA::tdma_solve( double* x ) {

	double		id = 1.0;

	// foward sweep
	U[0] = U[0] / D[0];
	rhs[0] = rhs[0] / D[0];


	for ( int i = 1; i < rank; i++ )
	{
		id = 1.0 / ( D[i] - U[i-1] * L[i] );
		U[i] = U[i]*id;
		rhs[i] = ( rhs[i] - L[i]*rhs[i-1] ) * id;

	}

	// back substitue
	x[rank - 1] = rhs[rank - 1];

	for ( int i = rank-2 ;  i >=0 ; i-- )
	{
		x[i] = rhs[i] - U[i] * x[i+1] ;

	}

};

void PIC_TDMA::tdma_solve_default() {

	tdma_solve(X);

};


void PIC_TDMA::tdma_verification() {

	int		n = rank;

	L[0] = 0.0; L[1] = -1.0; L[2] = -1.0; L[3] = -1.0;
	D[0] = 4.0; D[1] = 4.0; D[2] = 4.0; D[3] = 4.0;
	U[0] = -1.0; U[1] = -1.0; U[2] = -1.0; U[3] = 0.0;
	rhs[0] = 5.0; rhs[1] = 5.0; rhs[2] = 10.0; rhs[3] = 23.0;

	for ( int i = 0 ; i < rank ; i++ ) X[i] = .0;

	tdma_solve(X);

	for ( int i = 0; i < rank ; i++ ) cout << " rhs[" << i << "] = " << rhs[i] << endl;
	for ( int i = 0; i < rank ; i++ ) cout << " X[" << i << "] = " << X[i] << endl;

	// results: X = { 2, 3, 5, 7 }

}; // for rank = 4 !
