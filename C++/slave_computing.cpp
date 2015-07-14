#include <mpi.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "mpi_parameter.hpp"
#include "dw_dense_matrix.hpp"
#include "prcsn.h"
#include "CAR.hpp"

using namespace std; 

void slave_computing(CAR &model, const string &write_directory, const TDenseMatrix &bound, int perturbation_iteration, double perturbation_scale, int maximization_iteration)
{
	int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Status status;

	double *rPackage = new double [N_MESSAGE];
	double *sPackage = new double [N_MESSAGE];

	TDenseMatrix start_point; 	

	while (1)
        {
                MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (status.MPI_TAG == END_TAG)
                {
                  	delete []rPackage;
                  	delete []sPackage;
			return; 
		}
		else if (status.MPI_TAG == START_POINT_FILE_READY_TAG)
		{
			int nGroup = (int)rPackage[NGROUP_INDEX];
			int nParameter = (int)rPackage[NPARAMETER_INDEX]; 
			start_point.Resize(nGroup, nParameter); 
			string filename = write_directory + string("/") + START_POINT_FILE; 
			ifstream input_file(filename.c_str()); 
			if (!input_file)
			{
				cerr << "Error in opening " << filename << endl; 
				exit(-1); 
			}
			input_file >> start_point; 
			input_file.close(); 
		}
		else if (status.MPI_TAG == OPTIMIZATION_TAG)
		{
			int iGroup = (int)rPackage[IGROUP_INDEX]; 
			TDenseVector parameter = start_point.RowVector(iGroup); 
		
			// initial 
			TDenseVector logLSeries;
       			TDenseMatrix state;
        		double log_likelihood = model.LogLikelihood(logLSeries, state, parameter);
			TDenseVector initial_parameter; 
			initial_parameter.CopyContent(parameter); 	

			// optimal 
        		double optimal_log_likelihood = model.MaximizeLogLikelihood(parameter, bound.ColumnVector(0), bound.ColumnVector(1), std::vector<std::string>(0), perturbation_iteration, perturbation_scale, maximization_iteration);

			// write out
			if (log_likelihood > MINUS_INFINITY || optimal_log_likelihood > MINUS_INFINITY)
			{
				stringstream convert; 
				convert << iGroup; 
				string filename = write_directory + string("/solution.") + convert.str(); 
				ofstream output_file(filename.c_str(), ios::out); 
				if (!output_file)
				{
					cerr << "Error in opening " << filename << endl; 
					exit(-1); 
				}
				if (log_likelihood > MINUS_INFINITY)
				{ 
        				output_file << setprecision(20) << log_likelihood;
        				for (int i=0; i<initial_parameter.Dimension(); i++)
                				output_file<< "\t" << initial_parameter[i];
        				output_file << endl;
				}
				if (optimal_log_likelihood > MINUS_INFINITY)
				{
        				output_file << setprecision(20) << optimal_log_likelihood;
        				for (int i=0; i<parameter.Dimension(); i++)
                				output_file << "\t" << parameter[i];
        				output_file << endl;
				}
				output_file.close(); 
			}
		}
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD);
	}
}
