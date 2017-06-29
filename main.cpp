#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "PIC.h"

using namespace std;

int main( int argc, char** args )
{

	PIC*		pic_simulation;

	int			num_of_cells = 50;
	int			num_of_nodes;
	double		length = 0.05;

	double		num_density = 1E15;
	double		num_WN = 1E8;
	double		mass_e = Const_ERM_Mass;
	int			sign_electron = -1;
	double		mass_ion = 39.948 * Const_AMU_Mass; // argon
	int			sign_ion = 1;

	double		dt = 1e-12;
	double		phy_time = 1e-5;
	int			output_mod = 1000;
	int			output_bulk_mod = 1000;


	pic_simulation = new PIC();

	num_of_nodes = num_of_cells + 1;

	srand( (unsigned int) time(0) );


	// main process
	// (1) initialization: 
	//		1. the simulation domain
	pic_simulation->build_domain_uniform( length, num_of_cells );
	pic_simulation->TS_init( dt, phy_time );
	//		2. particles
	pic_simulation->particle_alloc( num_density, num_WN, mass_e, sign_electron, mass_ion, sign_ion );
	pic_simulation->show_init();
	//		3. electrostatics property
	pic_simulation->field_alloc( num_of_nodes, pic_simulation->Domain->delta_x );
	//		4. memory storage for counting number density of particles
	pic_simulation->N_alloc( num_of_nodes );
	//
	//
	// (2) initial sampling of particle motions: 
	//		1. position 
	//		2. velocity with maxwellian distribution
	pic_simulation->particle_sampling();
	//
	//
	// (3) loop of PIC simulation processes:

	do {

		pic_simulation->particle_tracking();
	//
	//		1. weight the charge density on each node
		pic_simulation->compute_rho_on_nodes();
	//		2. solve Poisson equation for electric potential, and resolve electric field on nodes
		pic_simulation->resolve_electro_stat_prop();
	//		3. estimate the electric field affected on each particle
		pic_simulation->compute_electric_force_on_particles();
	//		4. do the changes of each particle motion, by algorithm of leap-frog scheme 
		pic_simulation->update_particle_motion();
		pic_simulation->compute_particle_number_density_on_nodes();
	//		5. return to process 1. 
	//
		pic_simulation->TS_elapse_dt();

		pic_simulation->output_data_profile(output_mod);
		pic_simulation->output_data_bulk_values(output_bulk_mod);

	} while ( !pic_simulation->TS_reach_terminal() );


	return 0;
}

