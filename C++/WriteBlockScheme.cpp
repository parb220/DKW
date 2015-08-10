#include <fstream>
#include <vector>
#include "dw_dense_matrix.hpp"
using namespace std; 

bool WriteBlockScheme(const string &file_name, const vector<TIndex> &blocks)
{
	ofstream output; 
	output.open(file_name.c_str(), ios::out|ios::binary); 
	if (!output.is_open())
		return false; 
	int n_blocks = (int)(blocks.size()), n_indices, index; 
	output.write((char*)&(n_blocks),sizeof(int)); 
	for (int i_block=0; i_block<n_blocks; i_block++)
	{
		n_indices = blocks[i_block].size; 
		output.write((char*)&(n_indices), sizeof(int)); 
		for (int i_index=0; i_index<n_indices; i_index++)
		{
			index = blocks[i_block][i_index]; 
			output.write((char*)&(index), sizeof(int)); 
		}
	}
	output.close(); 
	return true; 
}
