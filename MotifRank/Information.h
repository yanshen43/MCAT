 // Information.h: interface for the CInformation class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_INFORMATION_H__FDF4C55F_CE16_41B1_B590_D3F32DA1ADDA__INCLUDED_)
#define AFX_INFORMATION_H__FDF4C55F_CE16_41B1_B590_D3F32DA1ADDA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "MotifInfo.h"
#include "MotifSet.h"
#include "RegionInfo.h"
#include "ResultInfo.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "Convertor.h"

#ifdef WIN32
#include "ERRNO.h"
#else
#include "error.h"
#endif

#define FILENAME 64
#define LINEMAX 512

class CInformation  
{
public:
	int initMatrix(FILE* hfile);
	int initMotif(FILE* hfile);
	int initSequences(FILE* hfile);

	void preProcess();
	CInformation();
	virtual ~CInformation();
	CMotifInfo* motif;
	CMotifSet* motifset;
	CRegionInfo* region;
	CResultInfo* result;
	char m_filename[FILENAME];
	CAlphabet* alp;
	bool rc;//rc=true if considering reverse string
private:
	int readMotifCase4(FILE* hfile);
	int parserLine(int* results, int maxresult, char* buffer, int maxlen, int* alphbetindex);
	int readSize(FILE* hfile, int index, int setindex);
	int readName(FILE* hfile, int index);
	int readLine(FILE* hfile, char line[], int maxlinelen, char endtoken, int readtype);
	int readMotifCase3(FILE* hfile);
	int readMotifCase2(FILE* hfile);
	int readMotifCase1(FILE* hfile);
	int readMotifCase0(FILE* hfile);
	char readChar(FILE* hfile, bool* readflag,bool* endflag);
	int readInt(FILE* hfile, bool* readflag, bool* endflag);
	double readDouble(FILE* hfile, bool* readflag, bool* endflag);
	void readLine(FILE* hfile, bool*readflag, bool* endflag, int* len, char* oneline);
	int calHit(CAlphabet::ALPHABET* m, int mlen);
	int calHitonSeq(char* seq, int seqlen, CAlphabet::ALPHABET* m, int mlen,CList* hitpostion);
	int calHit(CAlphabet::ALPHABET** m, int mnum, int mlen);
	int calHitonSeq(char* seq, int seqlen, CAlphabet::ALPHABET** m, int mnum, int mlen,CList* hitpostion);
};

#endif // !defined(AFX_INFORMATION_H__FDF4C55F_CE16_41B1_B590_D3F32DA1ADDA__INCLUDED_)
