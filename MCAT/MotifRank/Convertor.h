// Convertor.h: interface for the CConvertor class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CONVERTOR_H__26634343_7267_4449_89BF_DCA5326088E7__INCLUDED_)
#define AFX_CONVERTOR_H__26634343_7267_4449_89BF_DCA5326088E7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Alphabet.h"
#include "stdlib.h"
#include "stdio.h"
#include "MotifSet.h"
#include "RegionInfo.h"
#include "math.h"
#include "List.h"
#include "Node.h"

#define derive 0.01
#define MAXSETSIZE 1024

class CConvertor  
{
public:
	void transformation(CMotifSet* motifinfo, CRegionInfo* regioninfo, int index,bool brc);
	void duplicate(int index,int* results,int alphbetindex);
	void initAll(int index, int length);
	CConvertor(int motifnum);
	virtual ~CConvertor();
	
	int** PWM[ALPHA_NUM]; //PWM is the position weigth matrix of the motif
	                     //PWM[j][k] is the weight of alphabet j on postion k for motif
	double** PSSM[ALPHA_NUM];
	int** orderPSSM[ALPHA_NUM];//orderPSSM[j][k] is the index of j-th char at postion k of motif 
	int* motiflen;
	double* maxscore;
	int* hittime;
	char* motifset[MAXSETSIZE];
	char* rcmotifset[MAXSETSIZE];
	int rcmotifsetsize;
	int numofmotifs;
	int motifsetsize;
	int maxlen;
	int motifsetvalue[MAXSETSIZE];
	bool rc;
private:
	void showResult(int type);
	void initMotifInfo(CMotifSet* motifinfo, int index,int setindex);
	bool isZero(int index, int *currvalue);
	bool IsAdd(int value);
	void Motif2Int(int index);
	void addCompl(int index);
	void increase(int index, int *currvalue);
	void skipincrease(int index, int *currvalue);
	void addMotifWord(int index, int *currvalue);
	void BranchAndBound(int index);
	bool isSame(double d1, double d2);
	char convert(int one);
	int convert(char one);
	void calMaxScoreOne(char *promoter, int promoterlen, double *maxscore,int* hittime, CList** hitposition);
	void calMaxScore(CRegionInfo* regioninfo, CMotifSet* motifinfo);
	void initMotifSet();
	void orderedPSSM();
	void PWM2PSSM(CRegionInfo* region);
};

#endif // !defined(AFX_CONVERTOR_H__26634343_7267_4449_89BF_DCA5326088E7__INCLUDED_)
