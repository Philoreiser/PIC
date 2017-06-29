#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "PIC_Impl.h"

using namespace std;

#if !defined(__PIC_PARTICLE_H)
#define __PIC_PARTICLE_H

class PIC_Particle {

	public:

		//bool	is_first_timestep;

		int		cell_index;	//	to index the control cell which particle belongs to
		bool	is_act;		//	to specify whether particle still acts, not reaching boundary walls
		double	WN;			//	weighting number

		int		sgn_q;		//	sign of charge
		double	mass;		//	mass of each particle
		double	x;			//	position
		double	v;			//	velocity
		double	F;			//	external force

		// for foward leap-frog time marching scheme
		double	v_minus;
		double	v_plus;

		PIC_Particle();
		PIC_Particle( double m, int sgn, double weight_N );
		~PIC_Particle();

		void specify_property( double , int , double );
		void specify_motion_init( double , double );
		void specify_motion_init( double , double , double );

		void motion_electric_force( double E );

		void motion_force( double dt, double force );
		void motion_move( double dt );
		void motion_update_v( double dt, bool is_first_timestep );

};

#endif

