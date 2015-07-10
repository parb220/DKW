#ifndef _CDATA_FRBA_HEADER
#define _CDATA_FRBA_HEADER

#include <string>
#include <vector>
#include "dw_dense_matrix.hpp"
#include "CDate.hpp"

class TDenseVector; 
class TDenseMatrix; 
class CDate; 

class CData_FRBA 
{
public: // Data
	double Dt; 
	std::vector<CDate> date; 
        TDenseMatrix y_vec; // yield 
        TDenseMatrix ytips_vec; //tips 
        TDenseMatrix p_vec; // cpi 
        TDenseMatrix caps_vec;      // caps 
        TDenseMatrix floors_vec;    // floors
        TDenseMatrix IE_options;
        TDenseMatrix oil_vec;
        TDenseMatrix horLT;
        TDenseMatrix bcf_vec;
        TDenseMatrix bcfLT_vec;
        TDenseMatrix ieS_vec;
        TDenseMatrix ieL_vec;
        TDenseMatrix yields;
        TDenseMatrix alltips;
        TDenseMatrix swap_vec;
	TDenseMatrix yS_vec; 
        TDenseMatrix oil_vec2;

	TDenseVector MATgrid; 
	TDenseVector TIPSgrid;
	TDenseVector MATgrid_caps;
	TDenseVector STRIKEgrid_caps;
	TDenseVector MATgrid_floors; 
	TDenseVector STRIKEgrid_floors;
	TDenseVector MATgrid_options; 
	TDenseVector STRIKEgrid_options; 
	TDenseMatrix MATgrid_oil; 
	TDenseVector MATgrid_allyield; 
	TDenseVector MATgrid_alltips; 
	TDenseVector MATgrid_swap; 
	TDenseMatrix MATgrid_oil2; 
public:
        CData_FRBA(double _Dt=0.0);
        CData_FRBA(const CData_FRBA &right);
        ~CData_FRBA(){}

 	int NumberObservations() const { return (int)date.size(); }
        double TimeIncrement() const { return Dt; }

	bool LoadDate_FRBA(std::string &file, const CDate &start_date, const CDate &end_date); 
        bool LoadDataSeries_FRBA(std::vector<std::string> &file, const CDate &start_date, const CDate &end_date);
};

#endif
