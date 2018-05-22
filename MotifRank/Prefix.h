// Prefix.h: interface for the CPrefix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PREFIX_H__F2449C5E_DBE7_4D40_A7BA_488AC0BBF7D2__INCLUDED_)
#define AFX_PREFIX_H__F2449C5E_DBE7_4D40_A7BA_488AC0BBF7D2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Alphabet.h"
#include "stdlib.h"
#include "stdio.h"

class CPrefix  
{
public:
	void removeFirst(CPrefix* old);
	int getLastTwoSum();
	int getLastTwo();
	int getLast();
	void printprefix();
	void append(CAlphabet::ALPHABET* prefix, CAlphabet::ALPHABET end);//append end to prefix to init m_prefix
	void duplicate(CAlphabet::ALPHABET* motif, int len);
	bool isSame(CAlphabet::ALPHABET* motif, int s, int t);
	bool isSame(CAlphabet::ALPHABET* motif, int selfs, int selft, int s, int t);
	CPrefix(int len);
	virtual ~CPrefix();
	//member variable
	CAlphabet::ALPHABET* m_prefix;	//prefix of motif
	int length;						//length of m_prefix
};

#endif // !defined(AFX_PREFIX_H__F2449C5E_DBE7_4D40_A7BA_488AC0BBF7D2__INCLUDED_)
