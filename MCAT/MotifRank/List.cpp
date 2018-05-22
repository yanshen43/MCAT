// List.cpp: implementation of the CList class.
//
//////////////////////////////////////////////////////////////////////

#include "List.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CList::CList()
{
	head = NULL;
	end = NULL;
}

CList::~CList()
{
	CNode* delptr = NULL;
	while (head!=NULL) {
		delptr = head;
		head = head->next;
		delete delptr;
		delptr = NULL;
	}
}


void CList::insert(CNode *nodeptr)
{
	if (nodeptr==NULL) {
		return;
	}
	if (head==NULL) {
		head = end = nodeptr;
		nodeptr->next = NULL;
	} else {
		end->next = nodeptr;
		end = nodeptr;
		nodeptr->next = NULL;
	}
}
