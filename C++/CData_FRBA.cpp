#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include "prcsn.h"
#include "dw_dense_matrix.hpp"
#include "CDatedData.hpp"
#include "CData_FRBA.hpp"

CData_FRBA::CData_FRBA(double _Dt) :
Dt(_Dt),
date(),
y_vec(),
ytips_vec(),
p_vec(),
caps_vec(),
floors_vec(),
IE_options(),
oil_vec(),
horLT(),
bcf_vec(),
bcfLT_vec(),
ieS_vec(),
ieL_vec(),
yields(),
alltips(),
swap_vec(),
yS_vec(),
oil_vec2(),
MATgrid(),
TIPSgrid(),
MATgrid_caps(),
STRIKEgrid_caps(), 
MATgrid_floors(),
STRIKEgrid_floors(),
MATgrid_options(),
STRIKEgrid_options(),
MATgrid_oil(),
MATgrid_allyield(),
MATgrid_alltips(),
MATgrid_swap(),
MATgrid_oil2()
{}

CData_FRBA::CData_FRBA(const CData_FRBA &right) :
Dt(right.Dt),
date(right.date),
y_vec(right.y_vec),
ytips_vec(right.ytips_vec),
p_vec(right.p_vec),
caps_vec(right.caps_vec),
floors_vec(right.floors_vec),
IE_options(right.IE_options),
oil_vec(right.oil_vec),
horLT(right.horLT),
bcf_vec(right.bcf_vec),
bcfLT_vec(right.bcfLT_vec),
ieS_vec(right.ieS_vec),
ieL_vec(right.ieL_vec),
yields(right.yields),
alltips(right.alltips),
swap_vec(right.swap_vec),
yS_vec(right.yS_vec),
oil_vec2(right.oil_vec2),
MATgrid(right.MATgrid),
TIPSgrid(right.TIPSgrid),
MATgrid_caps(right.MATgrid_caps),
STRIKEgrid_caps(right.STRIKEgrid_caps),
MATgrid_floors(right.MATgrid_floors),
STRIKEgrid_floors(right.STRIKEgrid_floors),
MATgrid_options(right.MATgrid_options),
STRIKEgrid_options(right.STRIKEgrid_options),
MATgrid_oil(right.MATgrid_oil), 
MATgrid_allyield(right.MATgrid_allyield),
MATgrid_alltips(right.MATgrid_alltips),
MATgrid_swap(right.MATgrid_swap),
MATgrid_oil2(right.MATgrid_oil2)
{}

TDenseMatrix LoadDataMatrix(const std::string &file_name, const CDate &start_date, const CDate &end_date, const TIndex &selected_dim, double normalization=1.0, int missing_value_index = -1, double missing_value=0.0, bool take_log=false)
{
	CDatedDataSeries data_series; 
	if (!data_series.LoadFromFile(file_name))
		return TDenseMatrix(0,0); 
	
	CDatedDataSeries return_data_series; 	
	for (int t=0; t<(int)data_series.size(); t++)
	{
		if (data_series[t].GetDate() >= start_date && data_series[t].GetDate() <= end_date)
		{
			TDenseVector dataVec = data_series[t].GetData(selected_dim) * (1.0/normalization); 
			if (missing_value_index >= 0 && missing_value_index < dataVec.Dimension())
			{
				if (dataVec[missing_value_index] == missing_value)
					dataVec = Ones(dataVec.Dimension()) * MINUS_INFINITY;  
				else if (take_log)
				{
					for (int j=0; j<dataVec.Dimension(); j++)
						dataVec[j] = log(dataVec[j]); 
				}
			}
			else if (take_log)
			{
				for (int j=0; j<dataVec.Dimension(); j++)
					dataVec[j] = log(dataVec[j]); 
			}
			return_data_series.push_back(CDatedData(data_series[t].GetDate(), dataVec)); 
		}
	}
	return_data_series.Sort(); 	// sort according to date
	return return_data_series.GetDataSeries_Matrix(); 
}

bool CData_FRBA::LoadDate_FRBA(std::string &file_name, const CDate &start_date, const CDate &end_date)
{
	CDatedDataSeries data_series; 
	if (!data_series.LoadFromFile(file_name)) 
		return false; 
	date.clear(); 
	for (int t=0; t<(int)data_series.size(); t++)
	{
		if (data_series[t].GetDate() >= start_date && data_series[t].GetDate() <= end_date)
			date.push_back(data_series[t].GetDate()); 
	}
	return true; 
}

bool CData_FRBA::LoadDataSeries_FRBA(std::vector<std::string> &file, const CDate &start_date, const CDate &end_date)
{
	// y_vec 
	// FMAyield.txt, file[0]
	MATgrid.Resize(7);
	MATgrid[0] = 0.25; MATgrid[1] = 0.5; MATgrid[2] = 1.0; MATgrid[3] = 2.0; 
	MATgrid[4] = 4.0;  MATgrid[5] = 7.0; MATgrid[6] = 10.0;	
	// TIndex(0,6) correspond to 4:10 in data_FRBA.m because (1) 0 indexed in C++ and (2) first 3 columns correspond to date
	y_vec = LoadDataMatrix(file[0], start_date, end_date, TIndex(0,6), 100.0);
	if (y_vec.rows == 0 || y_vec.cols == 0)
		return false; 

	// ytips_vec
	// FMAtips.txt, file[1]
	TIPSgrid.Resize(3); 
	TIPSgrid[0] = 5; TIPSgrid[1] = 7; TIPSgrid[2] = 10; 
	// TIndex(0,2) corresponds to 4:6 in data_FRBA.m
	ytips_vec = LoadDataMatrix(file[1], start_date, end_date, TIndex(0,2), 100.0, 0, 0.0, false); // missing value at [0] == 0.0, no log taking 
	if (ytips_vec.rows == 0 || ytips_vec.cols == 0)
		return false; 
	// p_vec
	// FMAcpi.txt, file[2]
	// TIndex(0) corresponds to 4 in data_FRBA.m
	p_vec = LoadDataMatrix(file[2], start_date, end_date, TIndex(0), 1.0, 0, 0.0, true); // missing value at [0] == 0.0, log taking
	if (p_vec.rows == 0 || p_vec.cols == 0)
		return false; 
	
	//caps_vec
	// caps.txt, file[3]
	MATgrid_caps.Resize(2); 
	MATgrid_caps[0] = 1.0; MATgrid_caps[1] = 3.0; 
	STRIKEgrid_caps.Resize(3); 
	STRIKEgrid_caps[0] = 0.01; STRIKEgrid_caps[1] = 0.02; STRIKEgrid_caps[2] = 0.03; 
	// TIndex(0,5) corresponds to 4:9 in data_FRBA.m
	caps_vec = LoadDataMatrix(file[3], start_date, end_date, TIndex(0,5), 10000.0, 0, 0.0, true); // missing value at [0] = 0.0, log takeing 
	if (caps_vec.rows == 0 || caps_vec.cols == 0)
		return false; 
	
	//floors_vec
	//floors.txt, file[4]
	MATgrid_floors.Resize(2); 
	MATgrid_floors[0] = 1.0; MATgrid_floors[1] = 3.0; 
	STRIKEgrid_floors.Resize(3);
	STRIKEgrid_floors[0] = 0.01; STRIKEgrid_floors[1] = 0.02; STRIKEgrid_floors[2] = 0.03; 	
	// TIndex(0,5) corresponds to 4:9 in data_FRBA.m
	floors_vec = LoadDataMatrix(file[4], start_date, end_date, TIndex(0,5), 10000.0, 0, 0.0, true); // missing value at [0] = 0.0, log taking 
	if (floors_vec.rows == 0 || floors_vec.cols == 0)
		return false; 

	// IE_options
	// allyield.txt, file[5]
	MATgrid_options.Resize(2); 
	MATgrid_options[0] = 1.0; MATgrid_options[1] = 3.0; 
	STRIKEgrid_options.Resize(3); 
	STRIKEgrid_options[0] = 0.01; STRIKEgrid_options[1] = 0.02; STRIKEgrid_options[2] = 0.03; 
	// TIndex(2)(4) corresponds to [6,8] in data_FRBA.m
	TDenseMatrix yvec_options = LoadDataMatrix(file[5], start_date, end_date, TIndex(2)(4), 100.0); 
	if (yvec_options.rows == 0 || yvec_options.cols == 0)
		return false; 
	IE_options.Resize(yvec_options.rows, MATgrid_options.Dimension()*STRIKEgrid_options.Dimension()); 
	for (int t=0; t<IE_options.rows; t++)
	{
		if (caps_vec(t,0) == MINUS_INFINITY || floors_vec(t,0) == MINUS_INFINITY)
			IE_options.InsertRowMatrix(t,0,Ones(IE_options.cols)*MINUS_INFINITY); 
		else 
		{
			int tempID; 
			double MAT, STRIKE, tmp_IE_options; 
			for (int i=0; i<MATgrid_options.Dimension(); i++)
			{
				for (int j=0; j<STRIKEgrid_options.Dimension(); j++)
				{
					tempID= i*STRIKEgrid_options.Dimension()+j; 
					MAT = MATgrid_options(i); 
					STRIKE=STRIKEgrid_options(j); 
					tmp_IE_options = (exp(caps_vec(t,tempID)) - exp(floors_vec(t,tempID))) / exp(-MAT*yvec_options(t,i) ) + pow(1.0+STRIKE, MAT);
					IE_options(t,tempID) = 1.0/MAT * log(tmp_IE_options);   
				}
			}
		}
	}

	// oil_vec
	// oil.txt, file[6]
	TDenseMatrix oil_mat_grid = LoadDataMatrix(file[6], start_date, end_date, TIndex(5,7), 360.0); // TIndex(5,7) corresponds to [9:11] in data_FRBA.m
	if (oil_mat_grid.rows == 0 || oil_mat_grid.cols == 0)
		return false; 
	MATgrid_oil = Transpose(oil_mat_grid);
	// TIndex(1,3) corresponds to [5:7] in data_FRBA.m, [F2, F4, F13] contracts
	oil_vec = LoadDataMatrix(file[6], start_date, end_date, TIndex(1,3), 1.0, 0, 0.0, true); // miss value [0] = 0.0, log taking 
	if (oil_vec.rows == 0 || oil_vec.cols == 0)
		return false; 

	// horLT, bcf_vec and bcfLT_vec
	// FMAbcf.txt, file[7]
	// FMAbcfLT.txt, file[8]
	TDenseMatrix bcf05 = LoadDataMatrix(file[7], start_date, end_date, TIndex(0)); // TIndex(0) corresponds to 4 in data_FRBA.m, forecast of 5-year yield; 
	if (bcf05.rows == 0 || bcf05.cols == 0)
		return false; 
	TDenseMatrix bcf10 = LoadDataMatrix(file[7], start_date, end_date, TIndex(1)); // TIndex(1) corresponds to 5 in data_FRBA.m, forecast of 10-year yield;
	if (bcf10.rows == 0 || bcf10.cols == 0)
		return false; 
	TDenseMatrix bcf6_11 = LoadDataMatrix(file[8], start_date, end_date, TIndex(0)); // TIndex(0) corresponds to 4 in data_FRBA.m
	if (bcf6_11.rows == 0 || bcf6_11.cols == 0)
		return false; 
	TDenseMatrix bcf6_11_year = LoadDataMatrix(file[8], start_date, end_date, TIndex(1)); // TIndex(1) corresponds to 5 in data_FRBA.m
	if (bcf6_11_year.rows == 0 || bcf6_11_year.cols == 0)
		return false; 

	// horLT; 
	horLT.Resize(bcf6_11.rows, bcf6_11.cols);
	for (int t=0; t<horLT.rows; t++)
	{
		if (bcf6_11(t,0) == 0.0)
			horLT.InsertRowMatrix(t,0,Ones(horLT.cols)*MINUS_INFINITY); 
		else 
		{
			CDate bcf_11_date((int)bcf6_11_year(t,0), 12, 31); 
			double hortmp = (bcf_11_date-date[t])/365.25-5.0; 
			horLT(t,0) = hortmp; 
		}	
	}

	// bcf_vec; 
	bcf_vec.Resize(bcf05.rows, bcf05.cols+bcf10.cols); 
	bcf_vec.Insert(0, 0, bcf05*(1.0/100)); 
	bcf_vec.Insert(0, bcf05.cols, bcf10*(1.0/100)); 
	// checking missing value
	for (int t=0; t<bcf_vec.rows; t++)
	{
		if (bcf_vec(t,0) == 0.0)
			bcf_vec.InsertRowMatrix(t,0,Ones(bcf_vec.cols)*MINUS_INFINITY); 
	}

	// bcfLT_vec; 
	bcfLT_vec.Resize(bcf6_11.rows, bcf6_11.cols); 
	bcfLT_vec.Insert(0, 0, bcf6_11*(1.0/100)); 
	// checking missing value
	for (int t=0; t<bcfLT_vec.rows; t++)
	{
		if (bcfLT_vec(t,0) == 0.0)
			bcfLT_vec.InsertRowMatrix(t,0,Ones(bcfLT_vec.cols)*MINUS_INFINITY); 
	}

	// ieS_vec
	// FMAieS.txt, file[9]
	TDenseMatrix ieS_tmp = LoadDataMatrix(file[9], start_date, end_date, TIndex(), 100.0); // no missing value check, no log taking
	if (ieS_tmp.rows == 0 || ieS_tmp.cols == 0)
		return false; 
	ieS_vec = ieS_tmp.SubMatrix(0, ieS_tmp.rows-1, ieS_tmp.cols-1, ieS_tmp.cols-1); 
	for (int t=0; t<ieS_vec.rows; t++)
	{
		if (ieS_vec(t,0) == 0.0)
			ieS_vec.InsertRowMatrix(t,0,Ones(ieS_vec.cols)*MINUS_INFINITY); 
	}

	// ieL_vec
	// FMAieL.txt, file[10]
	TDenseMatrix ieL_tmp = LoadDataMatrix(file[10], start_date, end_date, TIndex(), 100.0); // no missing value check, no log taking
	if (ieL_tmp.rows == 0 || ieL_tmp.cols == 0)
		return false; 
	ieL_vec = ieL_tmp.SubMatrix(0, ieL_tmp.rows-1, ieL_tmp.cols-1, ieL_tmp.cols-1); 
	for (int t=0; t<ieL_vec.rows; t++)
	{
		if (ieL_vec(t,0) == 0.0)
			ieL_vec.InsertRowMatrix(t,0,Ones(ieL_vec.cols)*MINUS_INFINITY); 
	}

	// yields
	// allyield.txt, file[5]
	yields = LoadDataMatrix(file[5], start_date, end_date, TIndex(), 100.0); // no missing value check, no log taking
	MATgrid_allyield.Resize(32); 
	MATgrid_allyield[0] = 0.25; MATgrid_allyield[1] = 0.5; 
	for (int i=2; i<MATgrid_allyield.Dimension(); i++)
		MATgrid_allyield[i] = i-1.0; 

	// alltips
	// alltips.txt, file[11]
	alltips = LoadDataMatrix(file[11], start_date, end_date, TIndex(), 100.0); // no missing value check, no log taking
	MATgrid_alltips.Resize(29); 
	for (int i=0; i<MATgrid_alltips.Dimension(); i++)
		MATgrid_alltips[i] = i+2.0; 

	// swap_vec
	// swap.txt, file[12]
	MATgrid_swap.Resize(10); 
	for (int i=0; i<MATgrid_swap.Dimension(); i++)
		MATgrid_swap[i] = i+1.0;  
	swap_vec = LoadDataMatrix(file[12], start_date, end_date, TIndex(4)(6)(9), 1.0); // no missing value check, no log taking yet;
	for (int j=0; j<swap_vec.cols; j++)
	{
		for (int t=0; t<swap_vec.rows; t++)
		{
			if (swap_vec(t,j) == 999)
				swap_vec(t,j) = MINUS_INFINITY; 
		}
	}
	
	
	yS_vec.Resize(swap_vec.rows, swap_vec.cols); 
	for (int i=0; i<yS_vec.cols; i++)
	{
		double MAT = TIPSgrid(i); 
		TIndex id_allyield; 
		for (int j=0; j<MATgrid_allyield.Dimension(); j++)
		{
			if (MATgrid_allyield(j) == MAT)
				id_allyield += j; 
		}
		TIndex tmp_id_swap; 
		for (int t=0; t<swap_vec.rows; t++)
		{
			if (swap_vec(t,i) != MINUS_INFINITY)
				tmp_id_swap += t; 
		}
		swap_vec.Insert(tmp_id_swap, i, swap_vec.SubMatrix(tmp_id_swap, TIndex(i))*(1.0/100.0)); 
		yS_vec.Insert(tmp_id_swap, i, yields.SubMatrix(tmp_id_swap, id_allyield)-swap_vec.SubMatrix(tmp_id_swap,TIndex(i))); 	
	}

	// oil_vec2
	// oil.txt file[6]
	oil_vec2 = LoadDataMatrix(file[6], start_date, end_date, TIndex(0,3), 1.0, -1, 0.0, true); // no missing value check, taking log
	if (oil_vec2.rows == 0 || oil_vec2.cols == 0)
		return false; 
	oil_mat_grid = LoadDataMatrix(file[6], start_date, end_date, TIndex(4,7), 360.0); 
        if (oil_mat_grid.rows == 0 || oil_mat_grid.cols == 0)
                return false;
        MATgrid_oil2 = Transpose(oil_mat_grid);
	return true; 
}

