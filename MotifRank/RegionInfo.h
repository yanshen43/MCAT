// RegionInfo.h: interface for the CRegionInfo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_REGIONINFO_H__DEB17E5A_E40A_423D_BE9E_F7B7DB24A6F7__INCLUDED_)
#define AFX_REGIONINFO_H__DEB17E5A_E40A_423D_BE9E_F7B7DB24A6F7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Alphabet.h"
#define MAXREGION 10000
#define MAXSEQUENCE 100

class CRegionInfo  
{
public:
	void initRCpromoter();
	CRegionInfo();
	virtual ~CRegionInfo();
	int m_regiontype;	//0 means single region,1 means concatnated
	int m_alphabet;		//size of alphabet
	/*
	 *	basic_initprobability中一次存储的是前2位的联合概率值，也就是说
	 *	Pr(R[0]=0,R[1]=0),Pr(R[0]=0,R[1]=1),...,Pr(R[0]=0,R[1]=ALPHA_NUM-1)
	 *  Pr(R[0]=1,R[1]=0),Pr(R[0]=1,R[1]=1),...,Pr(R[0]=1,R[1]=ALPHA_NUM-1)
	 *  .....
	 *  Pr(R[0]=ALPHA_NUM-1,R[1]=0),Pr(R[0]=ALPHA_NUM-1,R[1]=1),...,Pr(R[0]=ALPHA_NUM-1,R[1]=ALPHA_NUM-1)
	 *  上面所有ALPHA_NUM^2项的总和是1
	 *  Pr(R[0]=x,R[1]=y)对应着basic_initprobability[x][y]
	 */
	double initprobability[ALPHA_NUM];				//initial probability
	/*
	 *	transitionmatrix中存储的是由前两位决定下一位的转移概率，也就是说
	 *  Pr(R[i+1]=0|R[i]=0,R[i-1]=0),Pr(R[i+1]=1|R[i]=0,R[i-1]=0),...,Pr(R[i+1]=ALPHA_NUM-1|R[i]=0,R[i-1]=0),
	 *	Pr(R[i+1]=0|R[i]=1,R[i-1]=0),Pr(R[i+1]=1|R[i]=1,R[i-1]=0),...,Pr(R[i+1]=ALPHA_NUM-1|R[i]=1,R[i-1]=0),
	 *  ...
	 *	Pr(R[i+1]=0|R[i]=ALPHA_NUM-1,R[i-1]=0),Pr(R[i+1]=1|R[i]=ALPHA_NUM-1,R[i-1]=0),...,Pr(R[i+1]=ALPHA_NUM-1|R[i]=ALPHA_NUM-1,R[i-1]=0),
	 *  
	 *	Pr(R[i+1]=0|R[i]=0,R[i-1]=1),Pr(R[i+1]=1|R[i]=0,R[i-1]=1),...,Pr(R[i+1]=ALPHA_NUM-1|R[i]=0,R[i-1]=1),
	 *	Pr(R[i+1]=0|R[i]=1,R[i-1]=1),Pr(R[i+1]=1|R[i]=1,R[i-1]=1),...,Pr(R[i+1]=ALPHA_NUM-1|R[i]=1,R[i-1]=1),
	 *  ...
	 *	Pr(R[i+1]=0|R[i]=ALPHA_NUM-1,R[i-1]=1),Pr(R[i+1]=1|R[i]=ALPHA_NUM-1,R[i-1]=1),...,Pr(R[i+1]=ALPHA_NUM-1|R[i]=ALPHA_NUM-1,R[i-1]=1),
	 *  
	 *  .......
	 *
	 *	Pr(R[i+1]=0|R[i]=0,R[i-1]=ALPHA_NUM-1),Pr(R[i+1]=1|R[i]=0,R[i-1]=ALPHA_NUM-1),...,Pr(R[i+1]=ALPHA_NUM-1|R[i]=0,R[i-1]=ALPHA_NUM-1),
	 *	Pr(R[i+1]=0|R[i]=1,R[i-1]=ALPHA_NUM-1),Pr(R[i+1]=1|R[i]=1,R[i-1]=ALPHA_NUM-1),...,Pr(R[i+1]=ALPHA_NUM-1|R[i]=1,R[i-1]=ALPHA_NUM-1),
	 *  ...
	 *	Pr(R[i+1]=0|R[i]=ALPHA_NUM-1,R[i-1]=ALPHA_NUM-1),Pr(R[i+1]=1|R[i]=ALPHA_NUM-1,R[i-1]=ALPHA_NUM-1),...,Pr(R[i+1]=ALPHA_NUM-1|R[i]=ALPHA_NUM-1,R[i-1]=ALPHA_NUM-1),
	 *  
	 *  Pr(R[i+1]=z|R[i]=y,R[i-1]=x)对应着basic_transitionmatrix[x][y][z]
	 *  这里的0~ALPHA_NUM-1和enum ALPHABET中0,1,2,3一一对应
	 */
	double transitionmatrix[ALPHA_NUM][ALPHA_NUM];	//transition probability matrix
	char* sequences[MAXSEQUENCE];//sequences[m_sequences][m_regionlen]
	char* rcsequences[MAXSEQUENCE];
	int m_sequences;	//number of sequences
	int m_regionlen;	//length of region
};

#endif // !defined(AFX_REGIONINFO_H__DEB17E5A_E40A_423D_BE9E_F7B7DB24A6F7__INCLUDED_)
