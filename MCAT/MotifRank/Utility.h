// Utility.h: interface for the CUtility class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_UTILITY_H__9C8CC0D2_D367_40E4_AD75_85A910EEB87E__INCLUDED_)
#define AFX_UTILITY_H__9C8CC0D2_D367_40E4_AD75_85A910EEB87E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Information.h"
#include "stdlib.h"
#include "string.h"
#define MAX_FILENAME 64

class CUtility  
{
public:
	CUtility();
	virtual ~CUtility();
	int readFile(CInformation* info,int argc, char* argv[]);
};

#endif // !defined(AFX_UTILITY_H__9C8CC0D2_D367_40E4_AD75_85A910EEB87E__INCLUDED_)
