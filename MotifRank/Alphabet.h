// Alphabet.h: interface for the CAlphabet class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ALPHABET_H__594C20E7_B5BD_470D_BC71_705BADEE7BAE__INCLUDED_)
#define AFX_ALPHABET_H__594C20E7_B5BD_470D_BC71_705BADEE7BAE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "string.h"
#include "stdio.h"

#define RELATION_NUM 4
#define ALPHA_NUM 4
#define EXALPHA_NUM 9
#define ALPHA_NUM_SQUARE ALPHA_NUM*ALPHA_NUM
class CAlphabet  
{
public:
	int getIndex(char one);
	enum ALPHABET{
		A=0,//Adenine
		C,	//Cytosine
		G,	//Guanine
		T,	//Thymine
		R,	//A or G
		W,  //A or T
		S,  //C or G
		Y,	//C or T
		N,	//spacer
		ERROR 
	};
	CAlphabet();
	virtual ~CAlphabet();
	void printalp(CAlphabet::ALPHABET a);
	bool checkValid(int type,char one);
	ALPHABET relation[EXALPHA_NUM][RELATION_NUM];
	int extendnum[EXALPHA_NUM];
	bool compatible[EXALPHA_NUM][ALPHA_NUM];
	int compatibleindex[EXALPHA_NUM][ALPHA_NUM];
};


#endif // !defined(AFX_ALPHABET_H__594C20E7_B5BD_470D_BC71_705BADEE7BAE__INCLUDED_)
