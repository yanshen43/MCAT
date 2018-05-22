// Rank.h: interface for the CRank class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RANK_H__39B90BA6_A899_4210_8F52_769D5DF6F85D__INCLUDED_)
#define AFX_RANK_H__39B90BA6_A899_4210_8F52_769D5DF6F85D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Alphabet.h"
#include "Information.h"
#include "Calculate.h"
#define MAX_MOTIFLEN 20

class CRank  
{
public:
	void rankresult(double* prob, int* index,CInformation* info);
	void doRanking(CInformation* info);
	void duplicate(CAlphabet::ALPHABET* destination, CAlphabet::ALPHABET* source, int size);
	CRank();
	virtual ~CRank();
	CAlphabet* alp;
private:
	void printmp(FILE* outputfile, CMotifInfo* m, double p);
};

#endif // !defined(AFX_RANK_H__39B90BA6_A899_4210_8F52_769D5DF6F85D__INCLUDED_)
