#include "PIC_Particle.h"

PIC_Particle::PIC_Particle() {

	is_act = true;
	WN = 1.0;
	sgn_q = 0;

	mass = .0;
	x = .0;
	v = .0;
	F = .0;

	v_minus = .0;
	v_plus = .0;

};


PIC_Particle::PIC_Particle( double m, int sgn, double weight_N ) {
	
	is_act = true;
	sgn_q = sgn;
	mass = m;
	WN = weight_N;
	
	x = .0;
	v = .0;
	F = .0;

	v_minus = .0;
	v_plus = .0;
};

PIC_Particle::~PIC_Particle() { };

void PIC_Particle::specify_property( double m, int sgn, double weight_N )
{
	mass = m;
	sgn_q = sgn;
	WN = weight_N;
};

void PIC_Particle::specify_motion_init( double position, double velocity )
{
	x = position;
	v = velocity;
}

void PIC_Particle::specify_motion_init( double position, double velocity, double force)
{
	x = position;
	v = velocity;
	F = force;
};

void PIC_Particle::motion_electric_force( double E ) {

	F = sgn_q * Const_e * E ;
};

void PIC_Particle::motion_force( double dt, double force ) {

	F = force;

	// update temporal velocity according to leap-frog scheme
	v_minus = v_plus ;
	v_plus = v_plus + force/mass * dt;

};

void PIC_Particle::motion_move( double dt )
{
	x = x + v_plus * dt;
};

void PIC_Particle::motion_update_v( double dt, bool is_first_timestep )
{

	if ( !is_first_timestep )
	{
		v_minus = v_plus ;
		v_plus = v_plus + F / mass * dt;

		v = (v_minus + v_plus) * 0.5;

		//cout << "j = " << j << endl;
	} else 
	{
		//v_minus = v;
		v_plus = v + F / mass * dt * 0.5;
		is_first_timestep = false;
		//cout << "first timestep j = " << j << endl;
	}

};

