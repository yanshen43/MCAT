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
	 *	basic_initprobability��һ�δ洢����ǰ2λ�����ϸ���ֵ��Ҳ����˵
	 *	Pr(R[0]=0,R[1]=0),Pr(R[0]=0,R[1]=1),...,Pr(R[0]=0,R[1]=ALPHA_NUM-1)
	 *  Pr(R[0]=1,R[1]=0),Pr(R[0]=1,R[1]=1),...,Pr(R[0]=1,R[1]=ALPHA_NUM-1)
	 *  .....
	 *  Pr(R[0]=ALPHA_NUM-1,R[1]=0),Pr(R[0]=ALPHA_NUM-1,R[1]=1),...,Pr(R[0]=ALPHA_NUM-1,R[1]=ALPHA_NUM-1)
	 *  ��������ALPHA_NUM^2����ܺ���1
	 *  Pr(R[0]=x,R[1]=y)��Ӧ��basic_initprobability[x][y]
	 */
	double initprobability[ALPHA_NUM];				//initial probability
	/*
	 *	transitionmatrix�д洢������ǰ��λ������һλ��ת�Ƹ��ʣ�Ҳ����˵
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
	 *  Pr(R[i+1]=z|R[i]=y,R[i-1]=x)��Ӧ��basic_transitionmatrix[x][y][z]
	 *  �����0~ALPHA_NUM-1��enum ALPHABET��0,1,2,3һһ��Ӧ
	 */
	double transitionmatrix[ALPHA_NUM][ALPHA_NUM];	//transition probability matrix
	char* sequences[MAXSEQUENCE];//sequences[m_sequences][m_regionlen]
	char* rcsequences[MAXSEQUENCE];
	int m_sequences;	//number of sequences
	int m_regionlen;	//length of region
};

#endif // !defined(AFX_REGIONINFO_H__DEB17E5A_E40A_423D_BE9E_F7B7DB24A6F7__INCLUDED_)
