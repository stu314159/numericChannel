#include "WallMountedBrick.h"
#include "OpenChannel3D.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

WallMountedBrick::WallMountedBrick(const int rk, const int sz, const string input_file, const string obst_file) :
					OpenChannel3D(rk,sz,input_file)
{
	obstFileName = obst_file;
	read_obst_file();
	update_snl_list();

}

WallMountedBrick::~WallMountedBrick(){

	delete [] obst_list;
}

void WallMountedBrick::read_obst_file(){
	// open the file for reading
	fstream obstFile(obstFileName.c_str(),ios::in);
	obstFile >> numObstNd;
	obst_list = new int[numObstNd];
	for(int i=0;i<numObstNd;i++){
		obstFile >> obst_list[i];
	}

	obstFile.close();
}

void WallMountedBrick::update_snl_list(){

	// for each node listed in the obst_list, the object must:

	// a) determine if that node is on its partition; and
	// b) its corresponding entry in the snl_list

	/*
	 * The above tasks are carried out by get_partition_node_num.
	 * if get_partition_node_num returns a positive number, then that is
	 * the appropriate node number for the solid node on that partition.
	 *
	 * if it returns -1, then the solid node is not on the partion and
	 * the snl need not be updated to reflect it.
	 */

	int l_nd;

	int count_l_nd = 0;
	for (int nd = 0; nd<numObstNd;nd++){
		l_nd = get_partition_node_num(obst_list[nd]);

		if(l_nd >= 0){
			count_l_nd++;
			snl[l_nd]=1;
		}


	}



}

int WallMountedBrick::get_partition_node_num(const int gNd){

	// if the gNd does not fall on this rank's partition, return -1.

	int physicalNdSlice = gNd/(Nx*Ny);
	int slicePos = gNd%(Nx*Ny);
	int l_nd = -1;



	/*
	 * in most cases, physicalNdSlice will be between firstSlice and lastSlice.
	 * In that case: this is simple.  It's the logically periodic streaming domain
	 * and the first/last rank that are of issue.
	 */

	if ((physicalNdSlice >= (firstSlice-HALO)) & (physicalNdSlice <= (lastSlice+HALO))){
		l_nd = (physicalNdSlice-firstSlice+HALO)*Nx*Ny + slicePos;

	}

	/*
	 * For WMB, this is enough since the obstruction will never be placed on the inlet or outlet
	 * slice.  I will leave this as a to-do, until I confirm the proper operation of all components
	 * of this class.  When done, I will make the class more general by properly dealing with the
	 * wrap-around halo slice issue.
	 */

	return l_nd;


}
