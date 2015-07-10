#ifndef _CDATED_DATA_HEADER_
#define _CDATED_DATA_HEADER_

#include <vector>
#include <string>
#include "dw_dense_matrix.hpp"
#include "CDate.hpp"

class CDate; 
class TDenseVector; 
class TDenseMatrix;
class TIndex; 

class CDatedData; 
class CDatedDataSorter; 
class CDatedDataSeries; 

class CDatedData
{
private:
	CDate date; 
	TDenseVector data;
public: 
	CDatedData(int _N=0); 
	CDatedData(const CDate &_date, const TDenseVector &_data, const TIndex &selectedDim=TIndex()); 
	CDatedData(const CDatedData &right, const TIndex &selectedDim=TIndex()); 
	CDatedData(const TDenseVector &date_data, const TIndex &selectedDim=TIndex()); 
	~CDatedData() {}
	int DataDimension() const; 
	CDatedData GetDatedData(const TIndex &selectedDim=TIndex()) const; 
	TDenseVector GetData(const TIndex &selectedDim=TIndex()) const; 
	CDate GetDate() const { return date; } 
	
	TDenseVector &Data() { return data; }
	CDate &Date() { return date; }

	const CDatedData &operator = (const CDatedData &right);
	const CDatedData &operator = (const TDenseVector &date_data);
	
	friend class CDatedDataSorter; 
};

class CDatedDataSorter
{
private:
	int index; 
public:
	CDatedDataSorter(int _index=-1); 
	bool operator() (const CDatedData &left, const CDatedData &right); 
}; 

class CDatedDataSeries : public std::vector<CDatedData>
{
public:
	CDatedDataSeries(int _T=0, int _N=0); 
	CDatedDataSeries(const CDatedDataSeries &right, const TIndex &selectedTime=TIndex(), const TIndex &selectedDim=TIndex()); 
	CDatedDataSeries(const std::vector<CDate> &date_series, const std::vector<TDenseVector> &data_series, const TIndex &selectedTime=TIndex(), const TIndex &selectedDim=TIndex()); 
	CDatedDataSeries(const std::vector<CDate> &date_series, const TDenseMatrix &data_matrix, const TIndex &selectedTime=TIndex(), const TIndex &selectedDim=TIndex());
	CDatedDataSeries(const TDenseMatrix &date_data_matrix, const TIndex &selectedTime=TIndex(), const TIndex &selectedDim=TIndex()); 
	~CDatedDataSeries() {}

	int ObservationLength() const { return (int)this->size(); }
	int DataDimension() const; 
	void Sort(int _index=-1); 
	TDenseMatrix GetDataSeries_Matrix(const TIndex &selectedTime=TIndex(), const TIndex &selectedDim=TIndex()) const; 
	std::vector<TDenseVector> GetDataSeries_Vector(const TIndex &selectedTime=TIndex(), const TIndex &selectedDim=TIndex()) const; 
	std::vector<CDate> GetDateSeries(const TIndex &selectedTime=TIndex()) const;  

	const CDatedDataSeries& operator=(const CDatedDataSeries &right); 
	const CDatedDataSeries& operator=(const TDenseMatrix &date_data_matrix);

	bool LoadFromFile(const std::string &filename); 
}; 


#endif
