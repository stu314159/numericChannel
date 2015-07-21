/*
 * WallMountedBrick.h
 *
 *  Created on: Jun 11, 2015
 *      Author: sblair
 */

#ifndef WALLMOUNTEDBRICK_H_
#define WALLMOUNTEDBRICK_H_

#include "OpenChannel3D.h"

using namespace std;

class WallMountedBrick : public OpenChannel3D
{
public:
	WallMountedBrick(const int rank, const int size, const string input_file,
			const string obstacle_file);
	~WallMountedBrick();

	string obstFileName;
	int numObstNd; // how many new solid nodes to add
	int * obst_list; // array for holding new solid nodes corresponding to the obstacle
	void read_obst_file();
	void update_snl_list();

private:
	int get_partition_node_num(const int gNd); // get the local (partition) node number from a global node number
	// this local number includes the halo nodes (as does all of the node lists...)


};



#endif /* WALLMOUNTEDBRICK_H_ */
