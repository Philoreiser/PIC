#include "PIC_Domain.h"


PIC_Cell::PIC_Cell(){ 

	node_L = .0;
	node_R = .0;
	dx = 1.0;
};

PIC_Cell::~PIC_Cell(){ 

};

void PIC_Cell::set( double pt_L, double pt_R, double delta_L ) {

	node_L = pt_L;
	node_R = pt_R;
	dx = delta_L;
};

int PIC_Cell::is_outside( double coor_x ) {

	int		feedback = 0;
/*
	if ( coor_x >= node_L && coor_x < node_R ) feedback = 0;
	else if ( coor_x < node_L ) feedback = 1;
	else if ( coor_x >= node_R ) feedback = 2;
*/

	if ( coor_x < node_L ) feedback = 1;
	else if ( coor_x > node_R ) feedback = 2;
	else feedback = 0;

	return feedback;
};


PIC_Domain::PIC_Domain() { 
	L = 0.05;
};


PIC_Domain::PIC_Domain( double length, int num_of_cells) {

	L = length;
	CellNum = num_of_cells;
	NodeNum = num_of_cells + 1;

	Cell = new PIC_Cell [num_of_cells];
	Node = new double [ num_of_cells + 1 ];

	delta_x = L / num_of_cells;
};

PIC_Domain::~PIC_Domain() { };


void PIC_Domain::mesh_uniform_partition() {

	int		i = 0;
	double	dx = L / CellNum ;

	//for ( i = 0; i < CellNum - 1; i++ )
	for ( i = 0; i < CellNum; i++ )
	{
		Cell[i].set( (double) i*dx, (double) (i+1)*dx, dx );

	}

	for ( i = 0; i < NodeNum ; i++ ) Node[i] = i*dx;
	

	delta_x = dx;
};

bool PIC_Domain::is_inside( double coor_x ) {

	if ( coor_x > .0 && coor_x < L ) return true;
	else return false;

};

int PIC_Domain::specify_cell_index( int seed, double coor_x ) {

//	int		i = seed;
	int		feedback = 0;
	int		index = seed;

	feedback = Cell[seed].is_outside(coor_x);

	if ( feedback == 0 )
	{
		// just do nothing

	} else if ( feedback == 1 ) // coor_x < node_L
	{
		for ( int i = seed-1; i >=0 ; i-- )
		{
			feedback = Cell[i].is_outside(coor_x);

			if ( !feedback )
			{
				index = i;
				break;
			}

		}

	} else if ( feedback == 2 ) // coor_x > node_R
	{
		for ( int i = seed+1; i < CellNum; i++ )
		{
			feedback = Cell[i].is_outside(coor_x);

			if ( !feedback )
			{
				index = i;
				break;
			}

		}

	}

	return index;

};

