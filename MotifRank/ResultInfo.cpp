  // ResultInfo.cpp: implementation of the CResultInfo class.
//
//////////////////////////////////////////////////////////////////////

#include "ResultInfo.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CResultInfo::CResultInfo()
{
	prob = NULL;
	multiprob = NULL;
}

CResultInfo::~CResultInfo()
{
	if (prob!=NULL) {
		release();
	}
	if (multiprob!=NULL) {
		releaseMulti();
	}
}

void CResultInfo::release()
{
	for(int i=0; i<length[0]; i++) {
		for(int j=0; j<length[1];j++) {
			for(int k=0; k<length[2]; k++)
				delete prob[i][j][k];
			free(prob[i][j]);
		}
		free(prob[i]);
	}
	free(prob);
	prob = NULL;
}

void CResultInfo::releaseMulti()
{
	for(int i=0; i<multilength[0]; i++) {
		for(int j=0; j<multilength[1];j++) {
			for(int k=0; k<multilength[2]; k++){
				for(int s=0; s<multilength[3];s++)
					delete multiprob[i][j][k][s];
				free(multiprob[i][j][k]);
			}
			free(multiprob[i][j]);
		}
		free(multiprob[i]);
	}
	free(multiprob);
	multiprob = NULL;
}

void CResultInfo::inittotal(int motiflen, int prefixtotalsize, int hittimes)
{
	if (motiflen<1 || prefixtotalsize<1 || hittimes<1) {
		return;
	}
	length[0] = motiflen+1;
	length[1] = hittimes;
	length[2] = prefixtotalsize;
	prob = (double****)malloc(sizeof(double***)*(length[0]));//for each postion in the region
	for(int i=0; i<length[0]; i++) {
		prob[i] = (double***)malloc(sizeof(double**)*length[1]);//for each prefix of motif
		for(int j=0; j<length[1]; j++) {
			prob[i][j] = (double**)malloc(sizeof(double*)*length[2]);
			for(int k=0; k<length[2]; k++)
				prob[i][j][k] = new double[ALPHA_NUM];
		}
	}
}

void CResultInfo::inittotal(int motiflen, int prefixtotalsize, int* hittimes) 
{
	if (motiflen<1 || prefixtotalsize<1 || hittimes==NULL) {
		return;
	}
	multilength[0] = motiflen+1;
	multilength[1] = hittimes[0]+1;
	multilength[2] = hittimes[1]+1;
	multilength[3] = prefixtotalsize;
	multiprob = (double*****)malloc(sizeof(double****)*multilength[0]);
	for(int i=0; i<multilength[0]; i++) {
		multiprob[i] = (double****)malloc(sizeof(double***)*multilength[1]);
		for(int j=0; j<multilength[1]; j++) {
			multiprob[i][j] = (double***)malloc(sizeof(double**)*multilength[2]);
			for(int k=0; k<multilength[2]; k++) {
				multiprob[i][j][k] = (double**)malloc(sizeof(double*)*multilength[3]);
				for(int s=0; s<multilength[3]; s++) {
					multiprob[i][j][k][s] = new double[ALPHA_NUM];
				}
			}
		}
	}
}

void CResultInfo::printself(CMotifInfo* motif)
{
	int k;
	for(int i=0; i<length[0]; i++) {
		for(int p=0; p<motif->hittimes; p++) {
			for(int j=motif->m_prefixtotalsize-1; j>-1;j--) {
				switch(motif->prefixset[j]->length) {
					case 0:
						for(k=0;k<ALPHA_NUM;k++) {
							motif->prefixset[j]->printprefix();
							printf("  length %d, hits %d, %.15f\n",i,p+1,prob[i][p][j][k]);	
							printf("\n");
						}
						break;
					default:
						motif->prefixset[j]->printprefix();
						printf("  length %d, hits %d, %.15f\n",i,p+1,prob[i][p][j][0]);
						printf("\n");
				}
			}
		}
		
	}
}


