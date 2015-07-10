#ifndef _CDATE_HEADER
#define _CDATE_HEADER

#include <iostream>
#include <istream>
#include <ostream>
class TDenseVector; 
class CDate
{
private:
	int year; 
	int month; 
	int day; 

public:
	bool operator < (const CDate &right) const; 
	bool operator <= (const CDate &right) const; 
	bool operator > (const CDate &right) const; 
	bool operator >= (const CDate &right) const; 
	bool operator == (const CDate &right) const;
	const CDate &operator = (const CDate &right); 
	
	friend std::istream& operator >> (std::istream &, CDate &); 
	friend std::ostream& operator << (std::ostream &, const CDate &); 

	CDate(int _year=0, int _month=0, int _day=0); 
	CDate(const TDenseVector &);
	CDate(const CDate &right); 
	~CDate(){} 

	int CountLeapYears() const; 
	int operator-(const CDate &right) const; 
};

#endif
