 // ResultInfo.h: interface for the CResultInfo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RESULTINFO_H__5ED8F716_427F_4023_9928_38207E08A2C9__INCLUDED_)
#define AFX_RESULTINFO_H__5ED8F716_427F_4023_9928_38207E08A2C9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Alphabet.h"
#include "MotifInfo.h"
#include "Prefix.h"
#include "stdlib.h"

class CResultInfo  
{
public:
	void printself(CMotifInfo* motif);
	void inittotal(int motiflen, int prefixtotalsize, int hittimes);
	void inittotal(int motiflen, int prefixtotalsize, int* hittimes);
	CResultInfo();
	virtual ~CResultInfo();	
	//member variable
	double**** prob;
	double***** multiprob;
	int length[3];//prob is 3 dimension: length[0]*length[1]*others
	int multilength[4];
private:
	void releaseMulti();
	void release();
	void release(int index);
};

#endif // !defined(AFX_RESULTINFO_H__5ED8F716_427F_4023_9928_38207E08A2C9__INCLUDED_)
