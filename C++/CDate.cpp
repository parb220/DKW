#include <cmath>
#include "dw_dense_matrix.hpp"
#include "CDate.hpp"

CDate::CDate(int _year, int _month, int _day) :
year(_year),
month(_month),
day(_day)
{}

CDate::CDate(const CDate &right) :
year(right.year),
month(right.month),
day(right.day)
{}

CDate::CDate(const TDenseVector &ymd) : 
year(0),
month(0),
day(0)
{
	if (ymd.Dimension() >= 1)
		year = (int)ymd[0];
	if (ymd.Dimension() >= 2)
		month = (int)ymd[1];
	if (ymd.Dimension() >= 3)
		day = (int)ymd[2]; 
}

bool CDate::operator == (const CDate &right) const
{
	if (this->year == right.year && this->month == right.month && this->day == right.day)
		return true; 
	else 
		return false; 
}

bool CDate::operator < (const CDate &right) const
{
	if (this->year < right.year || (this->year == right.year && this->month < right.month) || (this->year == right.year && this->month == right.month && this->day < right.day))
		return true; 
	else 
		return false; 
}

bool CDate::operator > (const CDate &right) const
{
	if (this->year > right.year || (this->year == right.year && this->month > right.month) || (this->year == right.year && this->month == right.month && this->day > right.day))
		return true; 
	else 
		return false; 
}

bool CDate::operator <= (const CDate &right) const 
{
	return !(operator>(right));
}

bool CDate::operator >= (const CDate &right) const
{
	return !(operator<(right)); 
}

const CDate& CDate::operator = (const CDate &right)
{
	this->year = right.year; 
	this->month = right.month; 
	this->day = right.day; 
	return *this; 	
}

std::istream& operator >> (std::istream &input_str, CDate &date)
{
	input_str >> date.year >> date.month >> date.day; 
	return input_str; 	
}

std::ostream& operator << (std::ostream &output_str, const CDate &date)
{
	output_str << date.year << "\t" << date.month << "\t" << date.day << "\t"; 
	return output_str; 
}

const int monthDays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

int CDate::CountLeapYears() const
{
	int years = this->year; 
	if (this->month <= 2)
		years --; 

	return years/4 - years/100 + years/400; 
}

int CDate::operator- (const CDate &right) const
{
	long int n_this = this->year*365 + this->day; 
	for (int i=0; i<this->month-1; i++)
		n_this += monthDays[i]; 
	n_this += this->CountLeapYears(); 

	long int n_right = right.year*365 + right.day; 
	for (int i=0; i<right.month-1; i++)
		n_right += monthDays[i]; 
	n_right += right.CountLeapYears(); 

	return (n_this - n_right); 
}
