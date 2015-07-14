#include <mpi.h>
#include <fstream>
#include "dw_dense_matrix.hpp"
#include "mpi_parameter.hpp"

using namespace std; 

void master_deploying(const string &write_directory, int nNode, int nGroup, bool if_initial_parameter, const TDenseVector &initial_p)
{
	double *sPackage= new double [N_MESSAGE];
        double *rPackage= new double [N_MESSAGE];
        MPI_Status status;

	// start_point; 
	vector<TDenseVector> start_point(nGroup); 
	if (if_initial_parameter)
	{
		for (int i=0; i<nGroup; i++)
			start_point[i] = initial_p; 
	}
	else 
	{
		for (int i=0; i<nGroup; i++) 
			start_point[i].RandomNormal(initial_p.Dimension()); 
	}

	string filename = write_directory + string("/") + START_POINT_FILE; 
	ofstream output_file(filename.c_str(), ios::out); 
	if (!output_file)
	{
		cerr << "Error in opening " << output_file << endl; 
		exit(-1); 
	}
	for (int i=0; i<(int)start_point.size(); i++)
	{
		for (int j=0; j<start_point[i].Dimension(); j++)
			output_file << start_point[i][j] << "\t"; 
		output_file << endl; 
	}
	output_file.close(); 

	// Telling slave node the file is ready to read. 
	sPackage[NGROUP_INDEX] = nGroup;
	sPackage[NPARAMETER_INDEX] = initial_p.Dimension(); 
	for (int i=1; i<nNode; i++)
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, START_POINT_FILE_READY_TAG, MPI_COMM_WORLD);	


	// deploying task
	vector<int> availableNode(nNode-1); 
	for (int i=0; i<(int)availableNode.size(); i++)
		availableNode[i] = i+1;

	for (int i=0; i<nGroup; i++)
	{
		sPackage[IGROUP_INDEX] = i; 
		if (availableNode.empty())
		{
			MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, OPTIMIZATION_TAG, MPI_COMM_WORLD, &status);
			availableNode.push_back(status.MPI_SOURCE); 
		}
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, availableNode.back(), OPTIMIZATION_TAG, MPI_COMM_WORLD);
		availableNode.pop_back();
	} 

	for (int i=(int)availableNode.size(); i<nNode-1; i++)
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, OPTIMIZATION_TAG, MPI_COMM_WORLD, &status);
	
	// destroy
	for (int i=1; i<nNode; i++)
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);	

	delete [] sPackage; 
	delete [] rPackage; 
}
