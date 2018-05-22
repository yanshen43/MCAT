  // Prefix.cpp: implementation of the CPrefix class.
//
//////////////////////////////////////////////////////////////////////

#include "Prefix.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPrefix::CPrefix(int len)
{
	if (len>0) {
		length = len;
		m_prefix = new CAlphabet::ALPHABET[len];
		if (m_prefix==NULL) {
			printf("fail to allocate memory.\n");
			exit(-1);
		}
	} else {
		length = 0;
		m_prefix = NULL;
	}
}

CPrefix::~CPrefix()
{
	if (m_prefix!=NULL) {
		delete m_prefix;
	}
}

void CPrefix::append(CAlphabet::ALPHABET* prefix, CAlphabet::ALPHABET end)
{
	if (m_prefix==NULL) {//null point
		return;
	}
	if (prefix==NULL) {
		m_prefix[0] = end;
	} else {
		memcpy(m_prefix, prefix, (length-1)*sizeof(CAlphabet::ALPHABET));
		m_prefix[length-1] = end;
	}
}

void CPrefix::duplicate(CAlphabet::ALPHABET* motif, int len)
{
	if (m_prefix==NULL || motif==NULL || len>this->length) {
		return;
	}
	memcpy(m_prefix, motif, len*sizeof(CAlphabet::ALPHABET));

}

bool CPrefix::isSame(CAlphabet::ALPHABET* motif, int s, int t)
{
	if (m_prefix==NULL && motif==NULL) {
		return true;
	} else if (m_prefix==NULL && motif!=NULL) {
		return false;
	} else if (m_prefix!=NULL && motif==NULL) {
		return false;
	}
	if (length!=t-s) {
		return false;
	}
	int index = s;
	while (index<t && motif[index]==m_prefix[index-s]) {
		index++;
	}
	if (index==t) {
		return true;
	} else 
		return false;
}

bool CPrefix::isSame(CAlphabet::ALPHABET* motif, int selfs, int selft, int s, int t)
{
	if (m_prefix==NULL && motif==NULL) {
		return true;
	} else if (m_prefix==NULL && motif!=NULL) {
		return false;
	} else if (m_prefix!=NULL && motif==NULL) {
		return false;
	}
	if ( selfs<0 || selft>length || (selft-selfs) != (t-s)) {
		return false;
	}
	int index = s;
	while (index<t && motif[index]==m_prefix[index-s+selfs]) {
		index++;
	}
	if (index==t) {
		return true;
	} else 
		return false;
}

void CPrefix::printprefix()
{
	for(int i=0; i<length; i++) {
		switch(m_prefix[i]) {
		case CAlphabet::A :
			printf("A");
			break;
		case CAlphabet::C :
			printf("C");
			break;
		case CAlphabet::G :
			printf("G");
			break;
		case CAlphabet::T :
			printf("T");
			break;
		default:
			printf("?");
		}
	}
}

int CPrefix::getLast()
{
	if (m_prefix!=NULL && length>0) {
		return m_prefix[length-1];
	} else 
		return 0;
}

int CPrefix::getLastTwoSum()
{
	if (m_prefix!=NULL && length>1) {
		return m_prefix[length-2]*ALPHA_NUM+m_prefix[length-1];
	} else 
		return 0;
}

int CPrefix::getLastTwo()
{
	if (m_prefix!=NULL && length>1) {
		return m_prefix[length-2];
	} else
		return 0;
}

void CPrefix::removeFirst(CPrefix *old)
{
	if (m_prefix==NULL || old==NULL) {
		return;
	}
	if (old->m_prefix==NULL) {
		return;
	}
	for(int i=0; i<length; i++) {
		m_prefix[i] = old->m_prefix[i+1];
	}
}
