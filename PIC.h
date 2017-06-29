#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "PIC_Impl.h"
#include "PIC_TDMA.h"
#include "PIC_Domain.h"
#include "PIC_Particle.h"
#include "PIC_Field.h"

using namespace std;

#if !defined(__PIC_H)
#define __PIC_H


class PIC {


	public:

		double			dt;
		double			phy_time_elapse;
		double			phy_time_terminal;

		double			init_ND;			//	initial number density of particles (electron-ion pairs )
		double			WN;					//	weighting number for number of particles
		int				Num_of_Parti;		//	number of simulated weighting particles
		int				Num_of_Parti_act;	//	number of active weighting particles

		double			IonMass;			//	mass of ion particle
		int				Sgn_ion;			//	charge sign of ion
		double			Ti;					//	initial ion temperature in unit [eV]

		double			ElectronMass;		//	mass of electron
		int				Sgn_e;				//	charge sign of electron (should be given as -1 )
		double			Te;					//	initial electron termperature in unit [eV]

		PIC_Particle*	Electron;			//	electrons
		PIC_Particle*	Ion;				//	ions

		PIC_Field*		ElectroStat;		//	electrostatics properties
		PIC_Domain*		Domain;				//	domain of simulation
		
		double*			N_e;				//	number density of electrons
		double*			N_i;				//	number density of ions


		PIC();
		~PIC();

		void show_init();

		void build_domain_uniform( double , int );
		void field_alloc( int );
		void field_alloc( int, double );
		void N_alloc(int);
		void particle_alloc( double , double , double , int , double , int );
		void particle_sampling();


		void TS_init( double , double );
		void TS_elapse_dt();
		bool TS_reach_terminal();

		void init_assembly();

		void particle_tracking();

		void weight_charge_to_nodes( PIC_Particle* , int NP );
		void compute_rho_on_nodes();

		void resolve_electro_stat_prop();

		void weight_electric_force_to_particles( PIC_Particle* , int NP );
		void compute_electric_force_on_particles();

		void particle_motion_processing( PIC_Particle* , int NP );
		void update_particle_motion();

		void N_set_zero();
		void weight_particle_to_nodes( PIC_Particle* Particle, int NP, double* N, double weight_N );
		void compute_particle_number_density_on_nodes();

		void output_electrostatics();

		void output_data_profile(int);
		void output_data_bulk_values(int);

};

#endif
