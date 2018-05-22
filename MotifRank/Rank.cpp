  // Rank.cpp: implementation of the CRank class.
//
//////////////////////////////////////////////////////////////////////

#include "Rank.h"
#include "Information.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CRank::CRank()
{
	alp = new CAlphabet();
}

CRank::~CRank()
{
	delete alp;
}

void CRank::doRanking(CInformation *info)
{
	//DWORD start, end;
	double* resultprob = new double[info->motifset->numberofmotif];
	int* resultindex = new int[info->motifset->numberofmotif];//record the result of ranking
	int i=0;
	memset(resultindex,0,sizeof(int)*info->motifset->numberofmotif);
	memset(resultprob,0,sizeof(double)*info->motifset->numberofmotif);
	//start = GetTickCount();
	info->preProcess();//pre precess the motifset
	//end  = GetTickCount();
	//printf("preProcess Time: %lums\n", end-start);
	CCalculate* cal = new CCalculate();
	for(i=0; i<info->motifset->numberofmotif; i++) {
		//info->result->init(info->region->m_regionlen,info->motifset->motifs[i]);
		if (info->motifset->motifs[i]->fdo == true) {
			info->motif = info->motifset->motifs[i];
			resultprob[i] =  cal->doCalculation(info);
		} else {
			resultprob[i] = 1;
		}
		
		//printmp(info->motifset->motifs[i],resultprob[i]);
	}	
	rankresult(resultprob,resultindex,info);
	//printf("result:\n");
	int index = 0;
	char outputfilename[FILENAME];
	memset(outputfilename,0,sizeof(char)*FILENAME);
	while (index<FILENAME && info->m_filename[index]!='.') {
		outputfilename[index] = info->m_filename[index];
		index++;
	}
	outputfilename[index++] = '.';
	outputfilename[index++] = 'p';
	outputfilename[index++] = 'v';
	outputfilename[index++] = 'a';
	outputfilename[index++] = 'l';
	//FILE* outputfile  = fopen(outputfilename,"w");
	for(i=0; i<info->motifset->numberofmotif; i++) {
		printmp(NULL, info->motifset->motifs[resultindex[i]],resultprob[resultindex[i]]);
	}
	//fclose(outputfile);
	delete cal;
	delete resultprob;
	delete resultindex;
}

void CRank::duplicate(CAlphabet::ALPHABET* destination, CAlphabet::ALPHABET* source, int size)
{
	if (destination==NULL || source==NULL) {
		return;
	}
	for(int i=0; i<size; i++)
		destination[i] = source[i];
}

void CRank::printmp(FILE* outputfile, CMotifInfo *m, double p)
{

	/*for(int i=0; i<m->m_motiflen; i++) {
		alp->printalp(m->m_motif[i]);
	}*/
	/*fprintf(outputfile,m->m_motifname);
	fprintf(outputfile, "\t%.40f\n",p);*/
	printf(m->m_motifname);
	printf("\t%.40f\n",p);
}

void CRank::rankresult(double *prob, int *index, CInformation *info)
{
	index[0] = 0;
	int pivot = 0;
	for(int i=1; i<info->motifset->numberofmotif; i++) {
		//insert the i-th motif, the probability is prob[i]
		pivot = 0;
		while (pivot<i && prob[index[pivot]]<prob[i]) {
			pivot++;
		}
		for(int j=i; j>pivot; j--)
			index[j] = index[j-1];
		index[pivot] = i;
	}
}
