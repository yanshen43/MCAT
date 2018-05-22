      // Calculate.cpp: implementation of the CCalculate class.
//
//////////////////////////////////////////////////////////////////////

#include "Calculate.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCalculate::CCalculate()
{
	prob_c = new double[ALPHA_NUM];
	prob_trans = new double[ALPHA_NUM];
	if (prob_c==NULL || prob_trans==NULL) {
		printf("fail to allocate memory.\n");
		exit(-1);
	}
}

CCalculate::~CCalculate()
{
	if (prob_c!=NULL) {
		delete prob_c;
	}
	if (prob_trans!=NULL) {
		delete prob_trans;
	}
}
/************************************************************************/
/* calculate the probability                                            */
/* according to information												*/
/************************************************************************/
double CCalculate::doCalculation(CInformation *info)
{
	if (info == NULL) {
		return 0;
	}
	m_info = info;
	if (m_info->motifset->numberofmotif==0) {
		return 1;
	}
	if (m_info->motifset->type<2) {
		calProb();
		if (m_info->region->m_regiontype==0) {
			return m_info->result->prob[m_info->region->m_regionlen%(m_info->motif->m_motiflen+1)][m_info->motif->hittimes-1][0][0];
		} else 
			return m_info->result->prob[(m_info->region->m_regionlen*m_info->region->m_sequences)%(m_info->motif->m_motiflen+1)][m_info->motif->hittimes-1][0][0];	
	} else {
		calMultiProb();
		if (m_info->region->m_regiontype==0) {
			return m_info->result->multiprob[m_info->region->m_regionlen%(m_info->motif->m_motiflen+1)][m_info->motif->motifsetshittimes[0]][m_info->motif->motifsetshittimes[1]][0][0];
		} else 
			return m_info->result->multiprob[(m_info->region->m_regionlen*m_info->region->m_sequences)%(m_info->motif->m_motiflen+1)][m_info->motif->motifsetshittimes[0]][m_info->motif->motifsetshittimes[1]][0][0];
	}
}

void CCalculate::calProb()
{
	int q;
	int regionlen = 0;
	if (m_info->region->m_regiontype==0) {
		regionlen = m_info->region->m_regionlen+1;
	} else 
		regionlen = m_info->region->m_regionlen*m_info->region->m_sequences+1;
	for(int i=0; i<regionlen; i++) {//for each position,0 is the last position
		for(int h=0; h<m_info->motif->hittimes; h++) {//for each hits
			for(int j = m_info->motif->m_prefixtotalsize-1; j>-1; j--) {//for each prefix
				switch(m_info->motif->prefixset[j]->length) {
				case 0://the length of prefix == 0
					if (i<m_info->motif->m_motiflen) {
						for(q=0; q<ALPHA_NUM; q++)
							m_info->result->prob[i%(m_info->motif->m_motiflen+1)][h][j][q] = 0;
					} else {
							for(q=0; q<ALPHA_NUM; q++)
								calcase0(i,h,j,q);
					}
					break;
				default://the length of prefix >= 1
					if (i<m_info->motif->m_motiflen) {// the length of the prefix is greater than i
						m_info->result->prob[i%(m_info->motif->m_motiflen+1)][h][j][0] = 0;
					} else {
							calcase1(i,h,j);
					}
				}

			}

		}		
	}
}


// it is the case that prefix length >=1
// prefix|?
void CCalculate::calcase1(int i, int h, int j)
{
	int j_prime, i_prime,h_prime;
	if (m_info->motif->prefixset[j]->length==m_info->motif->m_motiflen) {// hit a motif
		j_prime = m_info->motif->skiplist[0][j];
		i_prime = i - m_info->motif->prefixset[j]->length+m_info->motif->prefixset[j_prime]->length;
		if (i>m_info->region->m_regionlen && i%(m_info->region->m_regionlen)>0 && i%(m_info->region->m_regionlen)<m_info->motif->m_motiflen) {
			//the motif over the joint, so it should be omit
			h_prime = h;
		} else {
			if(h==0) {
				m_info->result->prob[i%(m_info->motif->m_motiflen+1)][h][j][0] = 1;
				return;
			}
			h_prime = h-1;
		}
		double next;
		j_prime = m_info->motif->skiplist[0][j];
		i_prime = i - m_info->motif->prefixset[j]->length+m_info->motif->prefixset[j_prime]->length;
		switch(m_info->motif->prefixset[j_prime]->length) {
		case 0:
			next = m_info->result->prob[i_prime%(m_info->motif->m_motiflen+1)][h_prime][j_prime][m_info->motif->prefixset[j]->getLast()];
			break;
		default:
			next = m_info->result->prob[i_prime%(m_info->motif->m_motiflen+1)][h_prime][j_prime][0]; 
		}
		m_info->result->prob[i%(m_info->motif->m_motiflen+1)][h][j][0] = next;
		return;
	}
	int k;
	for(k=0; k<ALPHA_NUM; k++) {
		//after appending k to the end of prefix j, the prefix skips to skiplist[k][j]
		int j_prime = m_info->motif->skiplist[k][j];
		int i_prime = i - (m_info->motif->prefixset[j]->length+1)+m_info->motif->prefixset[j_prime]->length;
		switch(m_info->motif->prefixset[j_prime]->length) {
		case 0:
			prob_c[k] = m_info->result->prob[i_prime%(m_info->motif->m_motiflen+1)][h][j_prime][k];
			break;
		default://the left size >= 1
			prob_c[k] = m_info->result->prob[i_prime%(m_info->motif->m_motiflen+1)][h][j_prime][0];
		}
		if(i>m_info->region->m_regionlen && (i-m_info->motif->prefixset[j]->length)%(m_info->region->m_regionlen)==0) {
			prob_trans[k] = m_info->region->initprobability[k];
		} else {
			prob_trans[k] = m_info->region->transitionmatrix[m_info->motif->prefixset[j]->getLast()][k];
		}
	}
	m_info->result->prob[i%(m_info->motif->m_motiflen+1)][h][j][0] = 0;
	for(k=0; k<ALPHA_NUM; k++) {
		m_info->result->prob[i%(m_info->motif->m_motiflen+1)][h][j][0] +=prob_c[k]*prob_trans[k];
	}
}

/************************************************************************/
/* q is the char before the prefix                                      */
/*
 * q|prefix|?
 */
/************************************************************************/
void CCalculate::calcase0(int i, int h, int j, int q)
{
	int k;
	for(k=0; k<ALPHA_NUM; k++) {
		//after appending k to the end of prefix j, the prefix skips to skiplist[k][j]
		int j_prime = m_info->motif->skiplist[k][j];
		int i_prime = i - (m_info->motif->prefixset[j]->length+1)+m_info->motif->prefixset[j_prime]->length;
		switch(m_info->motif->prefixset[j_prime]->length) {
		case 0:
			prob_c[k] = m_info->result->prob[i_prime%(m_info->motif->m_motiflen+1)][h][j_prime][k];
			break;
		default://the left size >= 1
			prob_c[k] = m_info->result->prob[i_prime%(m_info->motif->m_motiflen+1)][h][j_prime][0];
		}
		if (i>0 && i%m_info->region->m_regionlen==0) {
			prob_trans[k] = m_info->region->initprobability[k];
		} else
			prob_trans[k] = m_info->region->transitionmatrix[q][k];
	}
	m_info->result->prob[i%(m_info->motif->m_motiflen+1)][h][j][q] = 0;
	for(k=0; k<ALPHA_NUM; k++) {
		m_info->result->prob[i%(m_info->motif->m_motiflen+1)][h][j][q] +=prob_c[k]*prob_trans[k];
	}
	
}

void CCalculate::calMultiProb()
{
	int i,h0,h1;
	int regionlen = 0;
	if (m_info->region->m_regiontype==0) {
		regionlen = m_info->region->m_regionlen+1;
	} else 
		regionlen = m_info->region->m_regionlen*m_info->region->m_sequences+1;

	for(i=0; i<regionlen; i++) {
		for(h0=0; h0<m_info->motif->motifsetshittimes[0]+1; h0++) {//for each hits
			for(h1=0; h1<m_info->motif->motifsetshittimes[1]+1; h1++) {
				//use different prefixset according to different h1 and h2
				if (h0*h1==0 && (h0+h1)>0) {
					twoHitsCase1(i,h0,h1);
				} else {
					twoHitsCase0(i,h0,h1);
				}
			}
		}
	}
	
}

/************************************************************************/
/* use multiprefixset[2]                                                */
/************************************************************************/
void CCalculate::twoHitsCase0(int regionindex, int hits0, int hits1)
{
	int j,q;
	for(j=m_info->motif->m_multiprefixtotalsize[2]-1; j>-1; j--) {
		if (hits0==0 && hits1==0) {		
			for(q=0; q<ALPHA_NUM; q++)
				m_info->result->multiprob[regionindex%(m_info->motif->m_motiflen+1)][0][0][j][q] = 1;
		} else {//hits0>0, hits1>0
			if (m_info->motif->multiprefixset[2][j]->length>0) {
				m_info->result->multiprob[regionindex%(m_info->motif->m_motiflen+1)][hits0][hits1][j][0] = 
					calSingleCase1(regionindex,hits0,hits1,j,2);
				//m_info->motif->multiprefixset[2][j]->printprefix();
				//printf("  regionindex %d, m0 hits %d, m1 hits %d:  %.40f \n",regionindex,hits0,hits1, m_info->result->multiprob[regionindex%(m_info->motif->m_motiflen+1)][hits0][hits1][j][0]);
			} else {
				for(q=0; q<ALPHA_NUM; q++){
					m_info->result->multiprob[regionindex%(m_info->motif->m_motiflen+1)][hits0][hits1][j][q] = 
					calSingleCase0(regionindex,hits0,hits1,j,2,q);
					//m_info->motif->multiprefixset[2][j]->printprefix();
					//printf("  regionindex %d, m0 hits %d, m1 hits %d:  %.40f \n",regionindex,hits0,hits1, m_info->result->multiprob[regionindex%(m_info->motif->m_motiflen+1)][hits0][hits1][j][q]);
				}
			}

		}
	}
}


void CCalculate::twoHitsCase1(int regionindex, int hits0, int hits1)
{
	int j,q, setindex;
	if (hits0==0) {//calculate hits1
		setindex = 1;
	} else //calculate hits0
		setindex = 0;
	for(j=m_info->motif->m_multiprefixtotalsize[setindex]-1; j>-1; j--) {
		if (m_info->motif->multiprefixset[setindex][j]->length>0) {
			m_info->result->multiprob[regionindex%(m_info->motif->m_motiflen+1)][hits0][hits1][j][0] = 
				calSingleCase1(regionindex, hits0, hits1, j, setindex);
			//m_info->motif->multiprefixset[setindex][j]->printprefix();
			//printf("  regionindex %d, m0 hits %d, m1 hits %d:  %.40f \n",regionindex,hits0,hits1, m_info->result->multiprob[regionindex%(m_info->motif->m_motiflen+1)][hits0][hits1][j][0]);
		} else {
			for(q=0; q<ALPHA_NUM; q++) {
				m_info->result->multiprob[regionindex%(m_info->motif->m_motiflen+1)][hits0][hits1][j][q] = 
				calSingleCase0(regionindex, hits0, hits1, j, setindex, q);
				//m_info->motif->multiprefixset[setindex][j]->printprefix();
				//printf("  regionindex %d, m0 hits %d, m1 hits %d:  %.40f \n",regionindex,hits0,hits1, m_info->result->multiprob[regionindex%(m_info->motif->m_motiflen+1)][hits0][hits1][j][q]);
			}
		}
	}
}


double CCalculate::calSingleCase1(int regionindex, int hits0, int hits1, int prefixindex, int setindex)
{
	//compare motif length and region length
	if ((setindex<2 && m_info->motif->m_motifsetslen[setindex]>regionindex)
		||(setindex==2 && m_info->motif->m_motiflen > regionindex)) {//it's impossiable to hit
		return 0;
	}
	int delta_hits0, delta_hits1,newhits0,newhits1;
	int newregion;
	int newprefix;
	int newsetindex;
	double retvalue;
	int boundary;
	int diff;
	if(regionindex>m_info->region->m_regionlen 
		&& regionindex%(m_info->region->m_regionlen)>0 
		&& !(regionindex%(m_info->region->m_regionlen)>m_info->motif->multiprefixset[setindex][prefixindex]->length)) {
		boundary = regionindex%(m_info->region->m_regionlen);
	} else 
		boundary = 0;
	if (m_info->motif->m_fullhit[setindex][prefixindex]>-1) {//hits
		/////////////////////////decide newhits0 and newhits1////////////////
		delta_hits0 = m_info->motif->m_cutnum[setindex][0][0][prefixindex][boundary];
		delta_hits1 = m_info->motif->m_cutnum[setindex][0][1][prefixindex][boundary];
		////////////////////decide new hits///////////////////////////////
		if (hits0==0 && hits1>0) {
			newhits1 = ((hits1-delta_hits1)>0)?(hits1-delta_hits1):0;
			newhits0 = 0;
		} else if (hits0>0 && hits1==0) {
			newhits0 = ((hits0-delta_hits0)>0)?(hits0-delta_hits0):0;
			newhits1 = 0;
		} else {
			newhits1 = ((hits1-delta_hits1)>0)?(hits1-delta_hits1):0;
			newhits0 = ((hits0-delta_hits0)>0)?(hits0-delta_hits0):0;
		}
		////////////////////decide newprefix and newsetindex////////////
		if (hits0*hits1==0 || (newhits1>0 && newhits0>0)) {// skip int the orignal prefixset
			newprefix = m_info->motif->multiskiplist[setindex][0][prefixindex];
			newsetindex = setindex;
		} else if (newhits1==0) {//jump to prefix set of motif0
			newprefix = m_info->motif->multiskiplist[setindex][0][prefixindex];//first skip
			diff = m_info->motif->multiprefixset[setindex][prefixindex]->length-
				m_info->motif->multiprefixset[setindex][newprefix]->length;
			boundary = ((boundary-diff)>0)?(boundary-diff):0;
			newsetindex = 0;
			delta_hits0 = m_info->motif->m_degcutnum[newsetindex][newprefix][boundary];
			newprefix = m_info->motif->degskiplist[newsetindex][newprefix]; 
			newhits0 = ((newhits0-delta_hits0)>0)?(newhits0-delta_hits0):0;
		} else {//jump to prefix set of motif1
			newprefix = m_info->motif->multiskiplist[setindex][0][prefixindex];
			diff = m_info->motif->multiprefixset[setindex][prefixindex]->length-
				m_info->motif->multiprefixset[setindex][newprefix]->length;
			boundary = ((boundary-diff)>0)?(boundary-diff):0;
			newsetindex = 1;
			delta_hits1 = m_info->motif->m_degcutnum[newsetindex][newprefix][boundary];
			newprefix = m_info->motif->degskiplist[newsetindex][newprefix];
			newhits1 = ((newhits1-delta_hits1)>0)?(newhits1-delta_hits1):0;
		}
		///////////////////decide newregion/////////////////////////////////
		newregion = regionindex-(m_info->motif->multiprefixset[setindex][prefixindex]->length)
			+ (m_info->motif->multiprefixset[newsetindex][newprefix]->length);
		///////////////////decide return value//////////////////////////////
		if (m_info->motif->multiprefixset[newsetindex][newprefix]->length>0) {
			retvalue =
				m_info->result->multiprob[newregion%(m_info->motif->m_motiflen+1)][newhits0][newhits1][newprefix][0];
		} else {
			retvalue =
				m_info->result->multiprob[newregion%(m_info->motif->m_motiflen+1)][newhits0][newhits1][newprefix][m_info->motif->multiprefixset[setindex][prefixindex]->getLast()];
		}
	} else {//not hits
		int k;
		for(k=0; k<ALPHA_NUM; k++) {
			/////////////////////////decide newhits0 and newhits1////////////////
			delta_hits0 = m_info->motif->m_cutnum[setindex][k][0][prefixindex][0];
			delta_hits1 = m_info->motif->m_cutnum[setindex][k][1][prefixindex][0];
			////////////////////decide new hits///////////////////////////////
			if (hits0==0 && hits1>0) {
				newhits1 = ((hits1-delta_hits1)>0)?(hits1-delta_hits1):0;
				newhits0 = 0;
			} else if (hits0>0 && hits1==0) {
				newhits0 = ((hits0-delta_hits0)>0)?(hits0-delta_hits0):0;
				newhits1 = 0;
			} else {
				newhits1 = ((hits1-delta_hits1)>0)?(hits1-delta_hits1):0;
				newhits0 = ((hits0-delta_hits0)>0)?(hits0-delta_hits0):0;
			}
			////////////////////decide newprefix and newsetindex////////////
			if (hits0*hits1==0 || (newhits1>0 && newhits0>0)) {// skip int the orignal prefixset
				newprefix = m_info->motif->multiskiplist[setindex][k][prefixindex];
				newsetindex = setindex;
			} else if (newhits1==0) {//jump to prefix set of motif0
				//(oldprefix + addchar) first jumps to the prefix in prefixset
				//find that hits1 downs to zero, the the new prefix jumps to the prefixset of motif0
				newprefix = m_info->motif->multiskiplist[setindex][k][prefixindex];
				diff = m_info->motif->multiprefixset[setindex][prefixindex]->length+1-
					m_info->motif->multiprefixset[setindex][newprefix]->length;
				boundary = ((boundary-diff)>0)?(boundary-diff):0;
				newsetindex = 0;
				delta_hits0 = m_info->motif->m_degcutnum[newsetindex][newprefix][boundary];
				newprefix = m_info->motif->degskiplist[newsetindex][newprefix];
				newhits0 = ((newhits0-delta_hits0)>0)?(newhits0-delta_hits0):0;
			} else {//jump to prefix set of motif1
				newprefix = m_info->motif->multiskiplist[setindex][k][prefixindex];
				diff = m_info->motif->multiprefixset[setindex][prefixindex]->length+1-
				m_info->motif->multiprefixset[setindex][newprefix]->length;
				boundary = ((boundary-diff)>0)?(boundary-diff):0;
				newsetindex = 1;
				delta_hits1 = m_info->motif->m_degcutnum[newsetindex][newprefix][boundary];
				newprefix = m_info->motif->degskiplist[newsetindex][prefixindex];
				newhits1 = ((newhits1-delta_hits1)>0)?(newhits1-delta_hits1):0;
			}
			///////////////////decide newregion/////////////////////////////////
			newregion = regionindex-(m_info->motif->multiprefixset[setindex][prefixindex]->length+1)
				+(m_info->motif->multiprefixset[newsetindex][newprefix]->length);
			if (m_info->motif->multiprefixset[newsetindex][newprefix]->length>0) {
				prob_c[k] =
					m_info->result->multiprob[newregion%(m_info->motif->m_motiflen+1)][newhits0][newhits1][newprefix][0];
			} else {
				prob_c[k] =
					m_info->result->multiprob[newregion%(m_info->motif->m_motiflen+1)][newhits0][newhits1][newprefix][k];
			}
			if(regionindex>m_info->region->m_regionlen && 
				(regionindex-m_info->motif->multiprefixset[setindex][prefixindex]->length)%(m_info->region->m_regionlen)==0) {
				prob_trans[k] = m_info->region->initprobability[k];
			} else
				prob_trans[k] = m_info->region->transitionmatrix[m_info->motif->multiprefixset[setindex][prefixindex]->getLast()][k];
		}
		retvalue = 0;
		for(k=0; k<ALPHA_NUM; k++) {
			retvalue += prob_c[k]*prob_trans[k];
		}
	}
	return retvalue;
}

double CCalculate::calSingleCase0(int regionindex, int hits0, int hits1, int prefixindex, int setindex, int q)
{
	if ((setindex<2 && m_info->motif->m_motifsetslen[setindex]>regionindex)
		||(setindex==2 && m_info->motif->m_motiflen > regionindex)) {//it's impossiable to hit
		return 0;
	}
	int k;
	int delta_hits0, delta_hits1,newhits0,newhits1;
	int newregion;
	int newprefix;
	int newsetindex;
	double retvalue;
	for(k=0; k<ALPHA_NUM; k++) {
		/////////////////////////decide newhits0 and newhits1////////////////
		delta_hits0 = m_info->motif->m_cutnum[setindex][k][0][prefixindex][0];
		delta_hits1 = m_info->motif->m_cutnum[setindex][k][1][prefixindex][0];
		////////////////////decide new hits///////////////////////////////
		if (hits0==0 && hits1>0) {
			newhits1 = ((hits1-delta_hits1)>0)?(hits1-delta_hits1):0;
			newhits0 = 0;
		} else if (hits0>0 && hits1==0) {
			newhits0 = ((hits0-delta_hits0)>0)?(hits0-delta_hits0):0;
			newhits1 = 0;
		} else {
			newhits1 = ((hits1-delta_hits1)>0)?(hits1-delta_hits1):0;
			newhits0 = ((hits0-delta_hits0)>0)?(hits0-delta_hits0):0;
		}
		////////////////////decide newprefix and newsetindex////////////
		if (hits0*hits1==0 || (newhits1>0 && newhits0>0)) {// skip int the orignal prefixset
			newprefix = m_info->motif->multiskiplist[setindex][k][prefixindex];
			newsetindex = setindex;
		} else if (newhits1==0) {//jump to prefix set of motif0
			//(oldprefix + addchar) first jumps to the prefix in prefixset
			//find that hits1 downs to zero, the the new prefix jumps to 
			newprefix = m_info->motif->multiskiplist[setindex][k][prefixindex];
			newprefix = m_info->motif->degskiplist[0][newprefix]; 
			delta_hits0 = m_info->motif->m_degcutnum[0][newprefix][0];
			newhits0 = ((newhits0-delta_hits0)>0)?(newhits0-delta_hits0):0;
			newsetindex = 0;
		} else {//jump to prefix set of motif1
			newprefix = m_info->motif->multiskiplist[setindex][k][prefixindex];
			newprefix = m_info->motif->degskiplist[1][prefixindex];
			delta_hits1 = m_info->motif->m_degcutnum[1][newprefix][0];
			newhits1 = ((newhits1-delta_hits1)>0)?(newhits1-delta_hits1):0;
			newsetindex = 1;
		}
		///////////////////decide newregion/////////////////////////////////
		newregion = regionindex-(m_info->motif->multiprefixset[setindex][prefixindex]->length+1)
			+(m_info->motif->multiprefixset[newsetindex][newprefix]->length);
		if (m_info->motif->multiprefixset[newsetindex][newprefix]->length>0) {
			prob_c[k] =
				m_info->result->multiprob[newregion%(m_info->motif->m_motiflen+1)][newhits0][newhits1][newprefix][0];
		} else {
			prob_c[k] =
				m_info->result->multiprob[newregion%(m_info->motif->m_motiflen+1)][newhits0][newhits1][newprefix][k];
		}
		if (regionindex>0 && regionindex%m_info->region->m_regionlen==0) {
			prob_trans[k] = m_info->region->initprobability[k];
		} else
			prob_trans[k] = m_info->region->transitionmatrix[q][k];
	}
	retvalue = 0;
	for(k=0; k<ALPHA_NUM; k++) {
		retvalue += prob_c[k]*prob_trans[k];
	}
	return retvalue;
}

