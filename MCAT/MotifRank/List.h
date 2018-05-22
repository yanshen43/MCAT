// List.h: interface for the CList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LIST_H__14B75DD6_F238_4F9A_B5D5_140AAE775D76__INCLUDED_)
#define AFX_LIST_H__14B75DD6_F238_4F9A_B5D5_140AAE775D76__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "Node.h"
#include "stdlib.h"
#include "stdio.h"

class CList  
{
public:
	void insert(CNode* nodeptr);
	CList();
	virtual ~CList();
	CNode* head;
	CNode* end;
};

#endif // !defined(AFX_LIST_H__14B75DD6_F238_4F9A_B5D5_140AAE775D76__INCLUDED_)
