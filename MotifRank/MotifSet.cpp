// MotifSet.cpp: implementation of the CMotifSet class.
//
//////////////////////////////////////////////////////////////////////

#include "MotifSet.h"
#include "MotifInfo.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMotifSet::CMotifSet(int n)
{
	if (n>0) {
		motifs = (CMotifInfo**)malloc(sizeof(CMotifInfo*)*n);
		numberofmotif = n;
		for(int i=0; i<numberofmotif; i++) {
			motifs[i] = new CMotifInfo();
		}
	} else 
		exit(-1);
}

CMotifSet::~CMotifSet()
{
	if (motifs!=NULL) {
		for(int i=0; i<numberofmotif; i++){
			delete motifs[i];
		}
	}
	free(motifs);
}
