// Calculate.h: interface for the CCalculate class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CALCULATE_H__B04C6BB5_0274_4B70_8764_2FD5B08D8446__INCLUDED_)
#define AFX_CALCULATE_H__B04C6BB5_0274_4B70_8764_2FD5B08D8446__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Information.h"
#include "Alphabet.h"
#include "stdlib.h"

//#include "windows.h"

class CCalculate  
{
public:
	double doCalculation(CInformation* info);
	CCalculate();
	virtual ~CCalculate();
private:
	void twoHitsCase1(int regionindex, int hits1, int hits2);
	void twoHitsCase0(int regionindex, int hits1, int hits2);
	double calSingleCase0(int regionindex, int hits1, int hits2, int prefixindex, int setindex, int q);
	double calSingleCase1(int regionindex, int hits1, int hits2, int prefixindex, int setindex);
	void calMultiProb();
	void calcase0(int i, int h, int j, int q);
	void calcase1(int i, int h, int j);
	void calProb();
	CInformation* m_info;
	double* prob_c;
	double* prob_trans;
	//test variable
};

#endif // !defined(AFX_CALCULATE_H__B04C6BB5_0274_4B70_8764_2FD5B08D8446__INCLUDED_)
