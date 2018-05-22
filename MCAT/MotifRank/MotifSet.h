// MotifSet.h: interface for the CMotifSet class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MOTIFSET_H__6AC8D33A_35DB_4419_BCCA_DA8D31BF616C__INCLUDED_)
#define AFX_MOTIFSET_H__6AC8D33A_35DB_4419_BCCA_DA8D31BF616C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "MotifInfo.h"

class CMotifSet  
{
public:
	CMotifSet(int n);
	virtual ~CMotifSet();
	int numberofmotif;
	int type;
	CMotifInfo** motifs;
};

#endif // !defined(AFX_MOTIFSET_H__6AC8D33A_35DB_4419_BCCA_DA8D31BF616C__INCLUDED_)
