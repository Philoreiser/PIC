#include "PIC.h"

PIC::PIC() 
{ 
	cout << "**************************************" << endl;
	cout << "**                                  **" << endl;
	cout << "** Created a PIC Simulation Context **" << endl;
	cout << "**                                  **" << endl;
	cout << "**************************************" << endl;

	dt = 1e-12;
	phy_time_elapse = .0;
	phy_time_terminal = 1e-6;

	Ti = 0.0259;
	//Ti = 2.0;
	Te = 2.0;

};

PIC::~PIC() { };

void PIC::show_init() {
	/*
		double			init_ND;		//	initial number density of particles (electron-ion pairs )
		int				Num_of_Parti;	//	number of simulated weighting particles

		double			IonMass;		//	mass of ion particle
		int				Sgn_ion;		//	charge sign of ion
		double			ElectronMass;	//	mass of electron
		int				Sgn_e;			//	charge sign of electron (should be given as -1 )

	*/

	cout << "Initial number density of electrons/ions = " << init_ND << endl;
	cout << "Wieghting number for particles = " << WN << endl;
	cout << "Number of simulated weighting particles = " << Num_of_Parti << endl;
	cout << "Mass of one ion particle = " << IonMass << " kg" << " (" << IonMass/Const_AMU_Mass << " AMU)" << endl;
	cout << "Ion temperature = " << Ti << " eV" << " (" << Ti*Const_e/Const_Boltzmann << " K)" <<  endl;
	cout << "Electron rest mass = " << ElectronMass << " kg" << endl; 
	cout << "Electron temperature = " << Te << " eV" << endl;

};

void PIC::build_domain_uniform( double L, int num_of_cells) {

	Domain = new PIC_Domain( L, num_of_cells );
	Domain->mesh_uniform_partition();
};

void PIC::field_alloc( int num_of_nodes ) {

	ElectroStat = new PIC_Field();
	ElectroStat->alloc( num_of_nodes );
};

void PIC::field_alloc( int num_of_nodes, double dx ) {

	field_alloc( num_of_nodes );
	ElectroStat->delta_x = dx;
};

void PIC::N_alloc(int num_of_nodes) {

	N_e = new double[num_of_nodes];
	N_i = new double[num_of_nodes];

	for ( int i = 0; i < num_of_nodes; i++ )
	{
		N_e[i] = .0;
		N_i[i] = .0;
	}

};

void PIC::N_set_zero() {

	int		num_of_nodes = Domain->NodeNum ;

	for ( int i = 0; i < num_of_nodes; i++ )
	{
		N_e[i] = .0;
		N_i[i] = .0;
	}

};

void PIC::particle_alloc( double number_density, double weighting_number, double mass_electron, int sign_electron, double mass_ion, int sign_ion ) {

	init_ND = number_density;
	WN = weighting_number;
	ElectronMass = mass_electron;
	Sgn_e = sign_electron;
	IonMass = mass_ion;
	Sgn_ion = sign_ion;

//	double		ND = number_density;
//	double		WN = weighting_number;
//	int			nswp;	// number of simulated weighting particles

//	nswp = (int) ( ( ND * Domain->L ) / (WN * Domain->CellNum) ) ;
	Num_of_Parti = (int) ( ( number_density * Domain->L ) / ( weighting_number * Domain->CellNum) );

//	cout << "ND = " << ND << endl;
//	cout << "Total number of particles = " << ND * Domain->L << endl;
//	cout << "Total weighting number = " << WN * Domain->CellNum << endl;
//	cout << "nswp = " << nswp << endl;

	// allocation of particles
	Electron = new PIC_Particle [Num_of_Parti] ;
	Ion = new PIC_Particle [Num_of_Parti] ;

	for ( int i = 0; i < Num_of_Parti; i++ )
	{
		Electron[i].specify_property( ElectronMass, Sgn_e, weighting_number );
		Electron[i].is_act = true;

		Ion[i].specify_property( IonMass, Sgn_ion, weighting_number );
		Ion[i].is_act = true;
	}

	
};

void PIC::particle_sampling() {

	double		rf = .0;	// uniform random from 0 to 1
	double		x = .0;
	double		v = .0;

	for ( int i = 0; i < Num_of_Parti; i++ )
	{
		// random position distribution for particles
		rf = randf();
		x = rf * Domain->L ;

		// maxwellian velocity distribution for particles
		v = Velocity_MaxwellianDistri( IonMass, Ti, .0 );

		Ion[i].specify_motion_init( x, v);
		//cout << "x(" << i << ") = " << x << ", v(" << i << ") = " << v << endl;
	}


	for ( int i = 0; i < Num_of_Parti; i++ )
	{
		// random position distribution for particles
		rf = randf();
		x = rf * Domain->L;

		// maxwellian velocity distribution for particles
		v = Velocity_MaxwellianDistri( ElectronMass, Te, .0 );

		Electron[i].specify_motion_init( x, v);
	}

};

void PIC::TS_init( double timestep_size , double time_terminal ) {

	dt = timestep_size;
	phy_time_terminal = time_terminal;
	phy_time_elapse = .0;
};

void PIC::TS_elapse_dt() {

	phy_time_elapse += dt;

	//printf("phy_time_elapse = %e\n", phy_time_elapse );
};

bool PIC::TS_reach_terminal() {

	if ( phy_time_elapse >= phy_time_terminal ) return true;
	else return false;
};

void PIC::init_assembly() {
/*
	
	// prepare the algorithms
	// (1) initialization: 
	//		1. the simulation domain
	pic_simulation->build_domain_uniform( length, num_of_cells );
	pic_simulation->TS_init( dt, phy_time );
	//		2. particles
	pic_simulation->particle_alloc( num_density, num_WN, mass_e, sign_electron, mass_ion, sign_ion );
	pic_simulation->show_init();
	//		3. electrostatics property
	pic_simulation->field_alloc( num_of_nodes );
	//
	//
	// (2) initial sampling of particle motions: 
	//		1. position 
	//		2. velocity with maxwellian distribution
	pic_simulation->particle_sampling();
*/
};


void PIC::particle_tracking() {

	//cout << "time = " << phy_time_elapse << endl; 

	for ( int i = 0; i < Num_of_Parti; i++ )
	{
		if ( Electron[i].is_act ) 
		{
			// to track whether the particle reaches the boundary wall
			Electron[i].is_act = Domain->is_inside( Electron[i].x ); 

			// to track the cell of each particle moves to
			if ( Electron[i].is_act )
			{
				Electron[i].cell_index = Domain->specify_cell_index( Electron[i].cell_index, Electron[i].x );
				//cout << "i = " << i << " cell_index = " << Electron[i].cell_index << ", x = " << Electron[i].x << ", v = " << Electron[i].v << ", F = " << Electron[i].F  << endl;
			}

		}

		if ( Ion[i].is_act ) 
		{
			// to track whether the particle reaches the boundary wall
			Ion[i].is_act = Domain->is_inside( Ion[i].x );

			// to track the cell of each particle moves to
			if ( Ion[i].is_act )
			{
				Ion[i].cell_index = Domain->specify_cell_index( Ion[i].cell_index, Ion[i].x );
				//cout << "i = " << i << " cell_index = " << Ion[i].cell_index << ", x = " << Ion[i].x << ", v = " << Ion[i].v << ", F = " << Ion[i].F  << endl;
			}
			
		}
	}

};


void PIC::weight_charge_to_nodes( PIC_Particle* Particle , int NP ) {

	int		i_cell = 0;

	double	delta_x = Domain->delta_x;
	double	dx_plus = .0;
	double	dx_minus = .0;

	for ( int i = 0; i < NP; i++ )
	{
		if ( Particle[i].is_act )
		{
			i_cell = Particle[i].cell_index;

			//delta_x = Domain->Cell[i_cell].dx;
			dx_plus = Domain->Cell[i_cell].node_R - Particle[i].x ;
			dx_minus = Particle[i].x - Domain->Cell[i_cell].node_L ;

			ElectroStat->rho[ i_cell ] += dx_plus / delta_x * Particle[i].sgn_q ;
			ElectroStat->rho[ i_cell + 1 ] += dx_minus / delta_x * Particle[i].sgn_q ;
		}

	}

};

void PIC::compute_rho_on_nodes() {

	double	dx = Domain->delta_x;

	ElectroStat->set_rho_zero();

	// to weight number of charge particles
	weight_charge_to_nodes( Electron, Num_of_Parti );
	weight_charge_to_nodes( Ion, Num_of_Parti );

	// to multiply the unit of charge density
	for ( int i = 0; i < ElectroStat->NN; i++ ) ElectroStat->rho[i] *= WN * Const_e / dx;

};

void PIC::resolve_electro_stat_prop() {
	
		ElectroStat->solve_phi_by_TDMA_default();
		ElectroStat->resolve_field_default();

};

void PIC::weight_electric_force_to_particles( PIC_Particle* Particle, int NP ) {

	int			i_cell = 0;

	double		dx = Domain->delta_x;
	double		dx_plus = .0;
	double		dx_minus = .0;

	double		E = .0;

	for ( int i = 0; i < NP; i++ )
	{
		if ( Particle[i].is_act )
		{
			i_cell = Particle[i].cell_index;

			dx_plus = Domain->Cell[i_cell].node_R - Particle[i].x ;
			dx_minus = Particle[i].x - Domain->Cell[i_cell].node_L ;

			E =  ( ElectroStat->E[i_cell] * dx_plus + ElectroStat->E[i_cell+1] * dx_minus ) / dx ;

			Particle[i].motion_electric_force( E );

		}

	}

};

void PIC::compute_electric_force_on_particles() {

	weight_electric_force_to_particles( Electron, Num_of_Parti );
	weight_electric_force_to_particles( Ion, Num_of_Parti );
};

void PIC::particle_motion_processing( PIC_Particle* Particle, int NP ) {

	static bool is_first_timestep = true;

	for ( int i = 0; i < NP; i++ )
	{
		if ( Particle[i].is_act )
		{

			Particle[i].motion_update_v( dt, is_first_timestep );
			Particle[i].motion_move( dt );
		}

	}

	if ( is_first_timestep ) is_first_timestep = false;
};

void PIC::update_particle_motion() {

	particle_motion_processing( Electron, Num_of_Parti );
	particle_motion_processing( Ion, Num_of_Parti );

};

void PIC::weight_particle_to_nodes( PIC_Particle* Particle, int NP, double* N, double weight_N ) {

	int		i_cell = 0;

	double	delta_x = Domain->delta_x;
	double	dx_plus = .0;
	double	dx_minus = .0;

	for ( int i = 0; i < NP; i++ )
	{
		if ( Particle[i].is_act )
		{
			i_cell = Particle[i].cell_index;

			//delta_x = Domain->Cell[i_cell].dx;
			dx_plus = Domain->Cell[i_cell].node_R - Particle[i].x ;
			dx_minus = Particle[i].x - Domain->Cell[i_cell].node_L ;

			N[ i_cell ] += dx_plus / delta_x * weight_N / delta_x ;
			N[ i_cell + 1 ] += dx_minus / delta_x * weight_N / delta_x ;
		}

	}

};

void PIC::compute_particle_number_density_on_nodes() {

	N_set_zero();

	weight_particle_to_nodes( Electron, Num_of_Parti, N_e, WN );
	weight_particle_to_nodes( Ion, Num_of_Parti, N_i, WN );
};

void PIC::output_electrostatics() {

	ofstream		file;

	//FileOutput.open( ( ss_buffer.str() ).c_str() , ios::out | ios::app );
	file.open( "output_electrostatics.dat" , ios::out | ios::app );

	for ( int i = 0; i < Domain->NodeNum; i++ )
	{
		file << Domain->Node[i] << "\t" << ElectroStat->phi[i] << "\t" << ElectroStat->E[i] << "\t" << ElectroStat->rho[i]  << endl;
	}


	file.close();
};

void PIC::output_data_profile( int i_mod ) {

	static int		i_profile = 0;
	static int		i_call = 0;

	ofstream		file;
	stringstream	ss_buffer;

	ss_buffer.str( string() );
	ss_buffer.clear();
	
	i_call += 1;

	if ( i_call % i_mod == 0 )
	{
		i_profile += 1;

		ss_buffer << "pic_profile_" << i_profile << "_" << phy_time_elapse << ".dat" ;

		file.open( ( ss_buffer.str() ).c_str(), ios::out | ios::app );

		for ( int i = 0; i < Domain->NodeNum; i++ )
		{
			file << Domain->Node[i] << "\t" 
									<< ElectroStat->phi[i] << "\t" 
									<< ElectroStat->E[i] << "\t" 
									<< ElectroStat->rho[i] << "\t" 
									<< N_e[i] << "\t" 
									<< N_i[i]  << endl;
		}

		file.close();
	}

	//cout << "i_call = " << i_call << ", i_profile = " << i_profile << " : " << (ss_buffer.str() ).c_str() << endl;
};

void PIC::output_data_bulk_values(int i_mod) {
	

	static int		i_profile = 0;
	static int		i_call = 0;

	int				i_bulk = Domain->NodeNum / 2;

	ofstream		file;
	stringstream	ss_buffer;

	i_call += 1;

	ss_buffer.str( string() );
	ss_buffer.clear();

	if ( i_call % i_mod == 0 )
	{
		i_profile += 1;

		ss_buffer << "pic_bulk_values.dat" ;

		file.open( ( ss_buffer.str() ).c_str(), ios::out | ios::app );

		file << phy_time_elapse << "\t" << ElectroStat->phi[i_bulk] << "\t" << N_e[i_bulk] << "\t" << N_i[i_bulk] << endl;

		file.close();
	}
};

