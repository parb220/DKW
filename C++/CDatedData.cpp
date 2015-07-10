#include <algorithm>
#include <fstream>
#include <sstream>
#include "CDatedData.hpp"
#include "dw_dense_matrix.hpp"

CDatedData::CDatedData(int _N) : 
date(),
data(TDenseVector(_N,0.0))
{}

CDatedData::CDatedData(const CDate &_date, const TDenseVector &_data, const TIndex &selectedDim) : 
date(_date)
{
	if (selectedDim.Size() <= 0)
		data.CopyContent(_data); 
	else 
	{
		TIndex validDim; 
		for (int i=0; i<selectedDim.Size(); i++)
		{
			if (selectedDim[i] >=0 && selectedDim[i] < _data.Dimension())
				validDim += selectedDim[i]; 
		}
		data = _data.SubVector(validDim); 
	}
}

CDatedData::CDatedData(const CDatedData &right, const TIndex &selectedDim) : 
date(right.date)
{
	if ( selectedDim.Size() <= 0 )
		data.CopyContent(right.data); 
	else 
	{
		TIndex validDim; 
                for (int i=0; i<selectedDim.Size(); i++)
                {
                        if (selectedDim[i] >=0 && selectedDim[i] < right.DataDimension())
                                validDim += selectedDim[i];
                }
                data = right.data.SubVector(validDim);
	}
}

CDatedData::CDatedData(const TDenseVector &date_data, const TIndex &selectedDim) :
date(date_data)
{
	if (date_data.Dimension() > 3) // there are data besides date	
	{
		if (selectedDim.Size() <= 0)
			data = date_data.SubVector(3,date_data.Dimension()-1); 
		else 
		{
			TIndex validDim; 
			for (int i=0; i<selectedDim.Size(); i++)
			{
				if (selectedDim[i] >= 3 && selectedDim[i] < date_data.Dimension())
					validDim += selectedDim[i]; 
			}
			data = date_data.SubVector(validDim); 
		}
	}
}

int CDatedData::DataDimension() const
{
	return data.Dimension(); 
}
 
const CDatedData& CDatedData::operator = (const CDatedData &right)
{
	this->date = right.date; 
	this->data.CopyContent(right.data); 
	return *this; 
}

const CDatedData &CDatedData::operator = (const TDenseVector &date_data)
{
	CDatedData temp(date_data); 
	this->date = temp.date;
        this->data.CopyContent(temp.data);
        return *this;
}

CDatedData CDatedData::GetDatedData(const TIndex &selectedDim) const
{
	return CDatedData(*this, selectedDim); 
}

TDenseVector CDatedData::GetData(const TIndex &selectedDim) const
{
	CDatedData temp = GetDatedData(selectedDim); 
	return temp.data; 
}

CDatedDataSorter::CDatedDataSorter(int _index) : 
index(_index)
{}

bool CDatedDataSorter::operator()(const CDatedData &left, const CDatedData &right)
{
	if (index < 0 || index >= (left.DataDimension() <= right.DataDimension() ? left.DataDimension() : right.DataDimension() ) )
		return left.date < right.date;
        else
             	return left.data[index] < right.data[index];
}

CDatedDataSeries::CDatedDataSeries(int _T, int _N) :
std::vector<CDatedData>(_T, CDatedData(_N))
{}

CDatedDataSeries::CDatedDataSeries(const CDatedDataSeries &right, const TIndex &selectedTime, const TIndex &selectedDim)
{
	if (selectedTime.Size() <= 0)
	{
		for (int i=0; i<right.ObservationLength(); i++)
			this->push_back(CDatedData(right[i],selectedDim));
	}
	else 
	{
		for (int i=0; i<selectedTime.Size(); i++)
		{
			if (selectedTime[i] >=0 && selectedTime[i] < (int)(right.size()) )
				this->push_back(CDatedData(right[selectedTime[i]],selectedDim)); 
		}
	}
}

CDatedDataSeries::CDatedDataSeries(const std::vector<CDate> &date_series, const std::vector<TDenseVector> &data_series, const TIndex &selectedTime, const TIndex &selectedDim)
{
	int N = (int)(date_series.size() <= data_series.size() ? date_series.size() : data_series.size()); 
	if (selectedTime.Size() <= 0)
	{
		for (int i=0; i<N; i++)
			this->push_back(CDatedData(date_series[i],data_series[i],selectedDim)); 
	}
	else 
	{
		for (int i=0; i<selectedTime.Size(); i++)
		{
			if (selectedTime[i] >= 0 && selectedTime[i] < N)
				this->push_back(CDatedData(date_series[selectedTime[i]], data_series[selectedTime[i]], selectedDim)); 
		}
	}
}

CDatedDataSeries::CDatedDataSeries(const std::vector<CDate> &date_series, const TDenseMatrix &data_matrix, const TIndex &selectedTime, const TIndex &selectedDim)
{
	int N = (int)date_series.size() <= data_matrix.rows ? (int)date_series.size() : data_matrix.rows; 
	if (selectedTime.Size() <= 0)
	{
		for (int i=0; i<N; i++)
			this->push_back(CDatedData(date_series[i],data_matrix.RowVector(i,selectedDim))); 
	}
	else 
	{
		for (int i=0; i<selectedTime.Size(); i++)
		{
			if (selectedTime[i] >= 0 && selectedTime[i] < N)
				this->push_back(CDatedData(date_series[selectedTime[i]],data_matrix.RowVector(selectedTime[i],selectedDim) )); 	
		}
	}
}

CDatedDataSeries::CDatedDataSeries(const TDenseMatrix &date_data_matrix, const TIndex &selectedTime, const TIndex &selectedDim)
{
	if (selectedTime.Size() <= 0)
	{
		for (int i=0; i<date_data_matrix.rows; i++)
			this->push_back(CDatedData(date_data_matrix.RowVector(i),selectedDim)); 
	}
	else 
	{
		for (int i=0; i<selectedTime.Size(); i++)
		{
			if (selectedTime[i] >= 0 && selectedTime[i] < date_data_matrix.rows)
				this->push_back(CDatedData(date_data_matrix.RowVector(selectedTime[i]), selectedDim)); 
		}
	}	
}

int CDatedDataSeries::DataDimension() const
{
	if (ObservationLength())
		return this->operator[](0).DataDimension(); 
	else 
		return 0; 
}

void CDatedDataSeries::Sort(int _index)
{
	CDatedDataSorter comparator(_index); 
	sort(this->begin(), this->end(), comparator); 
}


std::vector<TDenseVector> CDatedDataSeries::GetDataSeries_Vector(const TIndex &selectedTime, const TIndex &selectedDim) const
{
	std::vector<TDenseVector> data_series;
	if (selectedTime.Size() <= 0)
	{
		for (int i=0; i<ObservationLength(); i++)
			data_series.push_back(this->operator[](i).GetData(selectedDim)); 
	}
	else 
	{
		for (int i=0; i<selectedTime.Size(); i++)
		{
			if (selectedTime[i] >= 0 &&  selectedTime[i] < ObservationLength())
				data_series.push_back(this->operator[](selectedTime[i]).GetData(selectedDim)); 
		}
	}
	return data_series; 
}


TDenseMatrix CDatedDataSeries::GetDataSeries_Matrix(const TIndex &selectedTime, const TIndex &selectedDim) const
{
	std::vector<TDenseVector> data_series = GetDataSeries_Vector(selectedTime, selectedDim); 
	if (data_series.empty())
		return TDenseMatrix(0,0); 
	int M = (int)data_series.size(), N=data_series[0].Dimension(); 
	if (N <= 0)
		return TDenseMatrix(M,0); 
	TDenseMatrix data_matrix(M,N); 
	for (int i=0; i<M; i++)
		data_matrix.InsertRowMatrix(i,0,data_series[i]); 
	return data_matrix; 
}

std::vector<CDate> CDatedDataSeries::GetDateSeries(const TIndex &selectedTime) const
{
	std::vector<CDate> date_series; 
	if (selectedTime.Size() <= 0)
	{
		for (int i=0; i<ObservationLength(); i++)
			date_series.push_back(this->operator[](i).GetDate()); 
	}
	else 
	{
		for (int i=0; i<selectedTime.Size(); i++)
		{
			if (selectedTime[i] >= 0 || selectedTime[i] < ObservationLength())
				date_series.push_back(this->operator[](selectedTime[i]).GetDate()); 
		}
	}
	return date_series; 
}

const CDatedDataSeries& CDatedDataSeries::operator=(const CDatedDataSeries &right)
{
	this->clear(); 
	for (int i=0; i<right.ObservationLength(); i++)
		this->push_back(right[i]); 
	return *this; 	
}

const CDatedDataSeries& CDatedDataSeries::operator=(const TDenseMatrix &date_data_matrix)
{
	CDatedDataSeries temp(date_data_matrix); 
	this->clear();
        for (int i=0; i<temp.ObservationLength(); i++)
                this->push_back(temp[i]);
        return *this;
}

bool CDatedDataSeries:: LoadFromFile(const std::string &filename)
{
	std::ifstream input_file(filename.c_str()); 
	if (!input_file)
		return false; 
	
	std::string line; 

	int nRows = 0, nCols = 0; 
	std::istringstream line_stream; 
	while (!input_file.eof())
	{
		line.clear(); 
		getline(input_file,line); 
		if (!line.empty())
		{
			if (nRows == 0)
			{
				line_stream.clear(); 
				line_stream.str(line); 
				double tempV; 
				while (line_stream >> tempV)
					nCols ++; 
			}
			nRows ++; 
		}
	}
	input_file.close(); 
	
	if (nRows == 0 || nCols == 0)
		return false; 

	input_file.open(filename.c_str()); 
	if (!input_file.is_open())
		return false; 
	TDenseMatrix date_data_matrix(nRows, nCols); 
	input_file>> date_data_matrix; 
	input_file.close(); 
	this->operator=(date_data_matrix); 
	return true; 
}
