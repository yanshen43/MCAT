// RegionInfo.cpp: implementation of the CRegionInfo class.
//
//////////////////////////////////////////////////////////////////////

#include "RegionInfo.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CRegionInfo::CRegionInfo()
{
	for(int i=0; i<MAXSEQUENCE; i++) {
		sequences[i] = NULL;
		rcsequences[i] = NULL;
	}
	m_sequences = 0;
	m_regionlen = 0;
}

CRegionInfo::~CRegionInfo()
{
	for(int i=0; i<m_sequences; i++) {
		if (sequences[i]!=NULL) {
			delete  sequences[i];
			sequences[i] = NULL;
		}
		if (rcsequences[i] != NULL) {
			delete rcsequences[i];
			rcsequences[i] = NULL;
		}
	}
}

void CRegionInfo::initRCpromoter()
{
	for(int j=0; j<m_sequences; j++) {
		rcsequences[j] = new char[m_regionlen];
		for(int i=0; i<m_regionlen; i++) {
			switch(sequences[j][i]) {
			case 'A':
			case 'a':
				rcsequences[j][m_regionlen-1-i] = 'T';
				break;
			case 'C':
			case 'c':
				rcsequences[j][m_regionlen-1-i] = 'G';
				break;
			case 'G':
			case 'g':
				rcsequences[j][m_regionlen-1-i] = 'C';
				break;
			case 'T':
			case 't':
				rcsequences[j][m_regionlen-1-i] = 'A';
				break;
			default:
				return;
			}
		}
	}
}
