#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#if !defined(__PIC_DOMAIN_H)
#define __PIC_DOMAIN_H

class PIC_Cell {
	public:
		double	node_L;
		double	node_R;
		double	dx;

		PIC_Cell();
		~PIC_Cell();

		void set( double , double , double );
		int is_outside( double coor_x );		// to classify whether a particle is outside cell or not
};

class PIC_Domain {
	
	public:

		double			L;			//	scale length in unit [m]

		int				CellNum;	//	number of cells
		int				NodeNum;	//	number of nodes

		PIC_Cell*		Cell;		//	cell
		double*			Node;

		double			delta_x;

		PIC_Domain();
		PIC_Domain( double , int );
		~PIC_Domain();

		void mesh_uniform_partition();
		bool is_inside( double );

		int specify_cell_index( int , double );
};



#endif
