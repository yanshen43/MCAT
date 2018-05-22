                               // Convertor.cpp: implementation of the CConvertor class.
//
//////////////////////////////////////////////////////////////////////

#include "Convertor.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CConvertor::CConvertor(int motifnum)
{
	numofmotifs = motifnum;
	int i;
	for(i=0; i<ALPHA_NUM; i++) {
		PWM[i] = (int**)malloc(sizeof(int*)*motifnum);
		PSSM[i] = (double**)malloc(sizeof(double*)*motifnum);
		orderPSSM[i] = (int**)malloc(sizeof(int*)*motifnum);
		for(int j=0; j<motifnum; j++) {
			PWM[i][j] = NULL;
			PSSM[i][j] = NULL;
			orderPSSM[i][j] = NULL;
		}
	}
	motiflen = new int[motifnum];
	maxscore = new double[motifnum];
	hittime = new int[motifnum];
	for(i=0; i<MAXSETSIZE; i++) {
		motifset[i] = NULL;
		rcmotifset[i] = NULL;
	}
}

CConvertor::~CConvertor()
{
	int i,j;
	for(i=0; i<ALPHA_NUM; i++) {
		if (PWM[i]!=NULL) {
			for(j=0; j<numofmotifs;j++) {
				if (PWM[i][j]!=NULL) {
					delete PWM[i][j];
					PWM[i][j] = NULL;
				}
			}
			free(PWM[i]);
			PWM[i] = NULL;
		}
		if (PSSM[i]!=NULL) {
			for(j=0; j<numofmotifs;j++) {
				if (PSSM[i][j]!=NULL) {
					delete PSSM[i][j];
					PSSM[i][j] = NULL;
				}
			}
			free(PSSM[i]);
			PSSM[i] = NULL;
		}
		if (orderPSSM[i]!=NULL) {
			for(j=0; j<numofmotifs;j++) {
				if (orderPSSM[i][j]!=NULL) {
					delete orderPSSM[i][j];
					orderPSSM[i][j] = NULL;
				}
			}
			free(orderPSSM[i]);
			orderPSSM[i] = NULL;
		}
	}
	if (motiflen!=NULL) {
		delete motiflen;
		motiflen = NULL;
	}
	if (maxscore!=NULL) {
		delete maxscore;
		maxscore = NULL;
	}
	if (hittime!=NULL) {
		delete hittime;
		hittime = NULL;
	}
	for(i=0; i<MAXSETSIZE; i++) {
		if (motifset[i]!=NULL) {
			delete motifset[i];
			motifset[i] = NULL;
		}
		if (rc) {
			if (rcmotifset[i]!=NULL) {
				delete rcmotifset[i];
				rcmotifset[i] = NULL;
			}
		}		
	}
}

void CConvertor::initAll(int index, int length)
{
	motiflen[index] = length;
	for(int i=0; i<ALPHA_NUM; i++) {
		PWM[i][index] = new int[length];
		memset(PWM[i][index],0,sizeof(int)*length);
		PSSM[i][index] = new double[length];
		memset(PSSM[i][index],0,sizeof(double)*length);
		orderPSSM[i][index] = new int[length];
		memset(orderPSSM[i][index],0,sizeof(int)*length);
	}
}

void CConvertor::duplicate(int index, int *results, int alphbetindex)
{
	if (alphbetindex<0 || alphbetindex>3) {
		return;
	}
	memcpy(PWM[alphbetindex][index],results,sizeof(int)*motiflen[index]);	
}

void CConvertor::transformation(CMotifSet* motifinfo, CRegionInfo* regioninfo, int index, bool brc)
{
	if (motifinfo==NULL || regioninfo==NULL) {
		return;
	}
	rc = brc;
	PWM2PSSM(regioninfo);
	//showResult(0);
	//showResult(1);
	orderedPSSM();
	initMotifSet();
	calMaxScore(regioninfo,motifinfo);
	for(int i=0; i<numofmotifs; i++) {
		BranchAndBound(i);
		initMotifInfo(motifinfo,i,index);
	}
}

void CConvertor::PWM2PSSM(CRegionInfo *region)
{
	int wsum,ch;
	double ratio;
	for(int mindex=0; mindex<numofmotifs; mindex++) {
		wsum = 0;
		for(ch=0; ch<ALPHA_NUM; ch++) {
			wsum += PWM[ch][mindex][0];
		}
		for(int pos=0; pos<motiflen[mindex]; pos++) {
			for(ch=0; ch<ALPHA_NUM; ch++) {
				//calculate PSSM[mindex][ch][pos]
				ratio = ((double)PWM[ch][mindex][pos]/(double)wsum+derive)/(1+4*derive);
				PSSM[ch][mindex][pos] = log(ratio/region->initprobability[ch]);
			}
		}
	}
}

void CConvertor::orderedPSSM()
{
	int mindex,ch,pos,comindex;
	double ordered[ALPHA_NUM];
	int orderedindex[ALPHA_NUM];
	for(mindex=0; mindex<numofmotifs; mindex++) {
		for(pos=0; pos<motiflen[mindex]; pos++) {
			comindex = 0;
			for(int i=0; i<ALPHA_NUM; i++) {
				ordered[i] = 100;
			}
			while (comindex<ALPHA_NUM) {
				//insert PSSM[comindex][mindex][pos]
				int pivot = 0;
				while (pivot<comindex && PSSM[comindex][mindex][pos]<ordered[pivot]) {
					pivot++;
				}
				//insert at pivot
				for(int i=comindex; i>pivot; i--) {
					ordered[i] = ordered[i-1];
					orderedindex[i] = orderedindex[i-1];
				}
				ordered[pivot] = PSSM[comindex][mindex][pos];
				orderedindex[pivot] = comindex;
				comindex++;
			}
			for(ch=0; ch<ALPHA_NUM; ch++) {
				orderPSSM[ch][mindex][pos] = orderedindex[ch];
			}
		}
	}	
}

void CConvertor::initMotifSet()
{
	//char* motifset[MAXSETSIZE];
	int i;
	maxlen = 0;
	for(i=0; i<numofmotifs; i++) {
		if ( maxlen<motiflen[i] ) {
			maxlen = motiflen[i];
		}
	}
	for(i=0; i<MAXSETSIZE; i++) {
		motifset[i] = new char[maxlen+1];
	}
	if (rc) {
		for(i=0; i<MAXSETSIZE; i++) {
			rcmotifset[i] = new char[maxlen+1];
		}
	}
	
}

void CConvertor::calMaxScore(CRegionInfo *regioninfo, CMotifSet *motifinfo)
{
	if (rc==false) {
		calMaxScoreOne(regioninfo->sequences[0], regioninfo->m_regionlen, maxscore,hittime,NULL);
	} else {
		double* maxscores[2];
		int* hittimes[2];
		CList** hitpositions[2];
		for(int k1=0; k1<2; k1++) {
			maxscores[k1] = new double[numofmotifs];
			hittimes[k1] = new int[numofmotifs];
			hitpositions[k1] = (CList**)malloc(sizeof(CList*)*numofmotifs);//for each motif there is a postion list
			for(int k2=0; k2<numofmotifs; k2++)
				hitpositions[k1][k2] = new CList();
			
		}
		calMaxScoreOne(regioninfo->sequences[0], regioninfo->m_regionlen,maxscores[0],hittimes[0],hitpositions[0]);
		calMaxScoreOne(regioninfo->rcsequences[0], regioninfo->m_regionlen, maxscores[1],hittimes[1],hitpositions[1]);
		int larger;
		for(int i=0; i<numofmotifs; i++) {
			larger = maxscores[0][i]>maxscores[1][i]?0:1;
			if (maxscores[0][i]==maxscores[1][i]) {
				this->maxscore[i] = maxscores[larger][i];
				this->hittime[i] = hittimes[0][i];
				CNode* pivot1 = hitpositions[1][i]->head;
				for(int j=0; j<hittimes[1][i]; j++) {//check whether add hittims[1]
					int index = 0;
					CNode* pivot0 = hitpositions[0][i]->head;
					while (index<hittimes[0][i] && (pivot0->position+pivot1->position+motiflen[i] != regioninfo->m_regionlen)) {
						index++;
						pivot0 = pivot0->next;
					}
					if (index==hittimes[0][i]) {//add it
						this->hittime[i]++;
					}
					pivot1 = pivot1->next;
				}
			} else {
				this->maxscore[i] = maxscores[larger][i];
				this->hittime[i] = hittimes[larger][i];
			}
			//printf("%s\t%f\t%d\n",motifs->motifname[i],motifs->maxscore[i],motifs->hittime[i]);
		}
		for(int d1=0; d1<2;d1++) {
			delete maxscores[d1];
			delete hittimes[d1];
			for(int d2=0; d2<numofmotifs; d2++) {
				delete hitpositions[d1][d2];
			}
			free(hitpositions[d1]);
		}
	}
}

void CConvertor::calMaxScoreOne(char *promoter, int promoterlen, double *maxscore, int *hittime, CList** hitposition)
{
	double currscore = 0;
	int ch = 0;
	for(int mindex=0; mindex<numofmotifs; mindex++) {//calculate maxcore for motif[mindex] on current promoter
		maxscore[mindex] = -100;
		hittime[mindex] = 0;
		for(int start=0; start<promoterlen-motiflen[mindex]+1; start++) {
			currscore = 0;
			for(int pos=0; pos<motiflen[mindex]; pos++) {
				ch = convert(promoter[start+pos]);
				currscore+= PSSM[ch][mindex][pos];
			}
			if (currscore>maxscore[mindex]) {
				maxscore[mindex] = currscore;
				if (hitposition!=NULL) {
					hitposition[mindex]->insert(new CNode(start));
				}
				hittime[mindex] = 1;
			} else if (isSame(currscore,maxscore[mindex])) {
				if (hitposition!=NULL) {
					hitposition[mindex]->insert(new CNode(start));
				}
				hittime[mindex]++;
			}
		}

	}
}

int CConvertor::convert(char one)
{
	switch(one) {
	case 'A':
	case 'a':
		return 0;
	case 'C':
	case 'c':
		return 1;
	case 'G':
	case 'g':
		return 2;
	case 'T':
	case 't':
		return 3;
	default:
		return -1;
	}
}

char CConvertor::convert(int one)
{
	switch(one) {
	case 0:
		return 'A';
	case 1:
		return 'C';
	case 2:
		return 'G';
	case 3:
		return 'T';
	default:
		return 'N';
	}
}

bool CConvertor::isSame(double d1, double d2)
{
	if (d1>d2 && (d1-d2)<0.000000001) {
		return true;
	} else if (d1<d2 && (d2-d1)<0.000000001) {
		return true;
	} else {
		return false;
	}
}

void CConvertor::BranchAndBound(int index)
{
	double currscore = 0;
	int pivot = 0;
	int realindex = 0;
	motifsetsize = 0;
	int* currvalue = new int[motiflen[index]];
	memset(currvalue, 0, sizeof(int)*motiflen[index]);//init to all one	
	do {
		//calculate the currscore
		pivot = 0;
		currscore = 0;
		while (pivot<motiflen[index]) {
			realindex = orderPSSM[currvalue[pivot]][index][pivot];
			currscore += PSSM[realindex][index][pivot];
			pivot++;
		}
		if (isSame(currscore,maxscore[index])) {
			addMotifWord(index,currvalue);
			skipincrease(index,currvalue);
		} else if (currscore<maxscore[index]) {//not add, skip increase
			skipincrease(index,currvalue);
			//increase(index,currvalue);
		} else {//add and currvalue increase one
			addMotifWord(index,currvalue);
			increase(index, currvalue);
		}
		if (motifsetsize==MAXSETSIZE) {//the upperbound
			motifsetsize = 0;
			delete currvalue;
			return;
		}
	} while(isZero(index,currvalue));
	delete currvalue;
	if (rc) {//should add compulsory
		addCompl(index);
	}
}

void CConvertor::addMotifWord(int index, int *currvalue)
{
	memset(motifset[motifsetsize],0,sizeof(char)*(maxlen+1));//all zero
	for(int i=0; i<motiflen[index]; i++) {
		motifset[motifsetsize][i] = convert(orderPSSM[currvalue[i]][index][i]);
	}
	motifsetsize++;
}

void CConvertor::skipincrease(int index, int *currvalue)
{
	//find the first pivot that currvalue[pivot]>0
	int pivot = 0;
	while (pivot<motiflen[index] && currvalue[pivot]==0) {
		pivot++;
	}
	if (pivot==motiflen[index]) {
		return;
	}
	int carry = 1 ;
	currvalue[pivot++] = 0;
	while (pivot<motiflen[index]) {
		currvalue[pivot] = (currvalue[pivot]+carry)%ALPHA_NUM;
		if (currvalue[pivot]==0) {
			carry = 1;
			pivot++;
		} else 
			return;
	}
}

void CConvertor::increase(int index, int *currvalue)
{
	int carry = 1;
	int pivot = 0;
	while (pivot<motiflen[index]) {
		currvalue[pivot] = (currvalue[pivot]+carry)%ALPHA_NUM;
		if (currvalue[pivot]==0) {
			carry = 1;
			pivot++;
		} else 
			return;
	}
}

void CConvertor::addCompl(int index)
{
	if (motifsetsize==0) {//zero case
		rcmotifsetsize = 0;
		return;
	}
	Motif2Int(index);
	rcmotifsetsize = 0;
	char* tmp = new char[maxlen+1];
	memset(tmp,0,sizeof(char)*(maxlen+1));
	int pivot, value;
	for(int i=0; i<motifsetsize; i++) {
		//reverse the i-th motif,char* rcmotifset[MAXSETSIZE]
		memset(tmp,0,sizeof(char)*maxlen);
		for(int j=0; j<motiflen[index]; j++) {
			if (motifset[i][j]=='A' || motifset[i][j]=='a') {
				tmp[motiflen[index]-j-1] = 'T';
			} else if (motifset[i][j]=='C' || motifset[i][j]=='c') {
				tmp[motiflen[index]-j-1] = 'G';
			} else if (motifset[i][j]=='G' || motifset[i][j]=='g') {
				tmp[motiflen[index]-j-1] = 'C';
			} else {
				tmp[motiflen[index]-j-1] = 'A';
			}
		}
		//calculate the value of this motifs
		pivot = 0;
		value = 0;
		while (pivot<motiflen[index]) {
			value = value*4+convert(tmp[pivot]);
			pivot++;
		}
		//check whether value is in motifsetvalue
		if (IsAdd(value)) {//should add
			memcpy(rcmotifset[rcmotifsetsize],tmp, sizeof(char)*(maxlen+1));
			rcmotifsetsize++;
		}
	}
	delete tmp;
}

void CConvertor::Motif2Int(int index)
{
	memset(motifsetvalue,0,sizeof(int)*MAXSETSIZE);
	int pivot,swap;
	for(int i=0; i<motifsetsize; i++) {//calculate 	motifs->motifsetvalue[i]
		pivot = 0;
		while (pivot<motiflen[index]) {
			motifsetvalue[i] = motifsetvalue[i]*4+convert(motifset[i][pivot]);
			pivot++;
		}		
	}
	for(int s = 1; s<motifsetsize; s++) {
		pivot = s;
		while (pivot>0 && motifsetvalue[pivot]<motifsetvalue[pivot-1]) {
			swap = motifsetvalue[pivot];
			motifsetvalue[pivot] = motifsetvalue[pivot-1];
			motifsetvalue[pivot-1] = swap;
			pivot--;
		}
	}
}

bool CConvertor::IsAdd(int value)
{
	//binary search in motifs->motifsetvalue
	int s = 0;
	int e = motifsetsize;
	int pivot;
	while (!(s>e)) {
		pivot = motifsetvalue[(s+e)/2];
		if (value==pivot) {
			return false;
		} else if (value>pivot) {
			s = (s+e)/2+1;
		} else if (value<pivot) {
			e = (s+e)/2-1;
		}
	}
	return true;
}	

bool CConvertor::isZero(int index, int *currvalue)
{
	int i = 0;
	while (i<motiflen[index] && currvalue[i]==0) {
		i++;
	}
	if (i==motiflen[index]) {//all zero
		return false;
	} else //not all zero, we need it
		return true;
}

void CConvertor::initMotifInfo(CMotifSet* motifinfo, int index, int setindex)
{
	//copy motifset and rcmotifset,hittime to motifinfo
	//init the size of set
	if (setindex==-1) {//one set
		if (rc) {
			motifinfo->motifs[index]->m_multinum = motifsetsize+rcmotifsetsize;
		} else {
			motifinfo->motifs[index]->m_multinum = motifsetsize;
		}
		motifinfo->motifs[index]->m_motiflen = motiflen[index];
		motifinfo->motifs[index]->type = 1;
	} else {
		if (rc) {
			motifinfo->motifs[index]->m_motifsetsnum[setindex] = motifsetsize+rcmotifsetsize;
		} else {
			motifinfo->motifs[index]->m_motifsetsnum[setindex] = motifsetsize;
		}
		motifinfo->motifs[index]->m_motifsetslen[setindex] = motiflen[index];
		motifinfo->motifs[index]->type = 2;
	}
	
	if (motifinfo->motifs[index]->m_multinum==0) {
		motifinfo->motifs[index]->fdo = false;
	} else {
		motifinfo->motifs[index]->initMultiMotifs(setindex);
		int i,k;
		for(i=0; i<motifsetsize; i++) {
			//insert motifset[i]
			motifinfo->motifs[index]->m_index = 0;
			for(k=0; k<motiflen[index];k++) {//copy from buffer to motif
				motifinfo->motifs[index]->insert(motifset[i][k],i,setindex);	
			}
		}
		if (rc) {
			for(i=0; i<rcmotifsetsize; i++) {
				motifinfo->motifs[index]->m_index = 0;	
				for(k=0; k<motiflen[index];k++) {//copy from buffer to motif
					motifinfo->motifs[index]->insert(rcmotifset[i][k],i+motifsetsize,setindex);	
				}
			}
		}
		motifinfo->motifs[index]->hittimes = hittime[index];
	}
}

void CConvertor::showResult(int type)
{
	for(int i=0; i<numofmotifs; i++) {
		for(int j=0; j<ALPHA_NUM; j++) {
			for(int k=0; k<motiflen[i]; k++) {
				if (type==0) {
					printf("%d\t",PWM[j][i][k]);
				} else if (type==1) {
					printf("%f\t",PSSM[j][i][k]);
				}
			}
			printf("\n");
		}
	}
}
