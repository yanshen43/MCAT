// Node.h: interface for the CNode class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_NODE_H__4940709E_EBD7_4006_9C1F_02E60E9E0DCE__INCLUDED_)
#define AFX_NODE_H__4940709E_EBD7_4006_9C1F_02E60E9E0DCE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CNode  
{
public:
	CNode(int p);
	virtual ~CNode();
	int position;
	CNode* next;
};

#endif // !defined(AFX_NODE_H__4940709E_EBD7_4006_9C1F_02E60E9E0DCE__INCLUDED_)
