#include "PIC_Impl.h"

double randf() {

	double	var = .0;
	
	var = (double) rand()/RAND_MAX ;

	return var;
};

double Var_NormalDistri( double Mean, double Sigma ){

	double	RandomVar = .0;

	double	x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if ( use_last )
	{
		y1 = y2;
		use_last = 0;

	} else
	{
		do {
			x1 = 2.0 * randf() - 1.0;
			x2 = 2.0 * randf() - 1.0;
			w = x1*x1 + x2*x2;

		} while ( w >= 1.0 );

		w = sqrt( ( -2.0 * log(w) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;

		use_last = 1;
	}

	//cout << "y1 = " << y1 << endl;
	RandomVar = Mean + y1*Sigma ;

	return RandomVar;
};

double Velocity_MaxwellianDistri( double mass, double T, double mu ) {
	
	double	MB_distri_v = .0;
	
	//double	Sigma = sqrt( Const_Boltzmann * T / mass );
	double	Sigma = sqrt( Const_e * T / mass );
	double	standard_variance = Var_NormalDistri( .0, 1.0 );

//	cout << "Sigma = " << Sigma << endl;
//	cout << "random_MB = " << random_MB << endl;

	MB_distri_v = Sigma * standard_variance + mu;

	return MB_distri_v;
};


void TDMA_solve( double* L, double* D, double* U, double* rhs, double* x ) {

};
