               // MotifInfo.cpp: implementation of the CMotifInfo class.
//
//////////////////////////////////////////////////////////////////////

#include "MotifInfo.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMotifInfo::CMotifInfo()
{
	m_motiflen = 0;
	fdo = true;
	m_motif = NULL;
	m_multimotifs = NULL;
	prefixset = NULL;
	m_index = 0;
	p_alpha = new CAlphabet();
	hittimes = 1;
	for(int i=0; i<ALPHA_NUM; i++) {
		skiplist[i] = NULL;
		for(int k=0; k<MOTIFSETS+1; k++) {
			multiskiplist[k][i] = NULL;
			for(int y=0; y<MOTIFSETS; y++)
				m_cutnum[k][i][y] = NULL;
		}
	}
	m_prefixsize = NULL;
	m_prefixsizestart = NULL;
	for(int j=0; j<MOTIFSETS; j++) {
		m_motifsets[j] = NULL;
		degskiplist[j] = NULL;
		m_degcutnum[j] = NULL; 
		motifsetshittimes[j] = 1;
	}
	for(int s = 0; s<MOTIFSETS+1; s++) {
		multiprefixset[s] = NULL;
		m_fullhit[s] = NULL;
	}
	//memset(motifsetshittimes, 0, sizeof(int)*MOTIFSETS);

}

CMotifInfo::~CMotifInfo()
{
	if (p_alpha!=NULL) {
		delete p_alpha;
	}
	if (m_motif!=NULL) {
		delete m_motif;
	}
	if (prefixset!=NULL) {
		for(int i=0; i<m_prefixtotalsize; i++) {
			if (prefixset[i]!=NULL) {
				delete prefixset[i];
			}
		}
		free(prefixset);
	}
	for(int j=0; j<ALPHA_NUM; j++) {
		if (skiplist[j]!=NULL) {
			delete skiplist[j];
		}
		for(int p=0; p<MOTIFSETS+1; p++) {
			if (multiskiplist[p][j]!=NULL) {
				delete multiskiplist[p][j];
			}
			for(int y =0; y<MOTIFSETS; y++) {
				if (m_cutnum[p][j][y]!=NULL) {
					for(int h=0; h<m_multiprefixtotalsize[p]; h++) {
						delete m_cutnum[p][j][y][h];
					}
					free(m_cutnum[p][j][y]);
				}
			}
		}
	}
	if(m_prefixsize!=NULL)
		delete m_prefixsize;
	if (m_prefixsizestart!=NULL) {
		delete m_prefixsizestart;
	}
	if (m_multimotifs!=NULL) {
		for(int k=0; k<m_multinum; k++)
			delete m_multimotifs[k];
		free(m_multimotifs);
	}
	for(int s=0; s<MOTIFSETS;s++) {
		if (m_motifsets[s]!=NULL) {
			for(int p=0; p<m_motifsetsnum[s];p++)
				delete m_motifsets[s][p];
			free(m_motifsets[s]);
		}
		if (degskiplist[s]!=NULL) {
			delete degskiplist[s];
		}
		if (m_degcutnum[s]!=NULL) {
			for(int t=0; t<m_multiprefixtotalsize[MOTIFSETS]; t++) {
				delete m_degcutnum[s][t];
			}
			free(m_degcutnum[s]);
		}
	}

	for(int q=0; q<MOTIFSETS+1; q++) {
		if (multiprefixset[q]!=NULL) {
			for(int x=0; x<m_multiprefixtotalsize[q]; x++)
				delete multiprefixset[q][x];
			free(multiprefixset[q]);
		}
		if (m_fullhit[q]!=NULL) {
			delete m_fullhit[q];
		}
	}
}

void CMotifInfo::initMotif()
{
	if (m_motiflen>0 && m_motif==NULL) {
		m_motif = new CAlphabet::ALPHABET[m_motiflen];
		m_index = 0;
		
	}
	if (m_prefixsize!=NULL) {
		delete m_prefixsize;
		m_prefixsize = NULL;
	}
	if (m_prefixsizestart!=NULL) {
		delete m_prefixsizestart;
		m_prefixsizestart = NULL;
	}
	m_prefixsize = new int[m_motiflen+1];
	m_prefixsizestart = new int[m_motiflen+1];
}

bool CMotifInfo::insert(char c,int motifindex,int setindex)
{
	CAlphabet::ALPHABET a;
	switch(c) {
	case 'A':
	case 'a':
		a = CAlphabet::A;
		break;
	case 'C':
	case 'c':
		a = CAlphabet::C;
		break;
	case 'G':
	case 'g':
		a = CAlphabet::G;
		break;
	case 'T':
	case 't':
		a = CAlphabet::T;
		break;
	case 'R':
	case 'r':
		a = CAlphabet::R;
		break;
	case 'W':
	case 'w':
		a = CAlphabet::W;
		break;
	case 'S':
	case 's':
		a = CAlphabet::S;
		break;
	case 'Y':
	case 'y':
		a = CAlphabet::Y;
		break;
	case 'N':
	case 'n':
		a = CAlphabet::N;
		break;
	default:
		a = CAlphabet::ERROR;
	}
	if (a==CAlphabet::ERROR) {
		return false;
	} else {
		if (motifindex==-1) {//one motif
			if (m_index==m_motiflen) {
				return false;
			} else
				m_motif[m_index++] = a;
		} else if (setindex==-1) {//one set
			if (m_index==m_motiflen) {
				return false;
			} else
				m_multimotifs[motifindex][m_index++] = a;
		} else {
			if (m_index==m_motifsetslen[setindex]) {
				return false;
			}
			m_motifsets[setindex][motifindex][m_index++] = a;
		}
			
		return true;
	}
}


void CMotifInfo::preProcess()
{
	switch(type) {
	case 0:
		initPrefixset();
		initSkiplist();
		break;
	case 1:
		preProcessMulti();
		break;
	case 2:
		preProcessMotifSets();
		break;
	default:
		return;
	}
}

void CMotifInfo::initPrefixset()
{
	//calculate the size of prefixset
	initPrefixsize();
	prefixset = (CPrefix**) malloc(m_prefixtotalsize*sizeof(CPrefix*));
	//init prefix of length 0
	prefixset[0] = new CPrefix(0);
	int i,j,k;
	int startindex = 0;
	int endindex = m_prefixsize[0];
	int indexofpre = m_prefixsize[0];
	for(i = 1; i<m_motiflen+1; i++) {
		//init prefix of length i
		//the number of prefix of length i is m_prefixsize[i]
		for(j = 0; j<p_alpha->extendnum[m_motif[i-1]]; j++) {
			//for each prefix of length i-1, append relation[m_motif[i]][j] to the end
			for(k = startindex; k<endindex; k++) {
				prefixset[indexofpre] = new CPrefix(i);
				prefixset[indexofpre]->append(prefixset[k]->m_prefix,p_alpha->relation[m_motif[i-1]][j]);
				indexofpre++;
			}
		}
		startindex = endindex;
		endindex += m_prefixsize[i];
	}
}

void CMotifInfo::initSkiplist()
{
	int i,j;
	for(i=0; i<ALPHA_NUM; i++) {
		skiplist[i] = new int[m_prefixtotalsize];
		for(j=0; j<m_prefixtotalsize; j++) {
			if (prefixset[j]->length<m_motiflen) {
				skiplist[i][j] = findprefix(j, prefixset[j], (CAlphabet::ALPHABET)i);
			} else 
				skiplist[i][j] = findprefix(j, prefixset[j], CAlphabet::ERROR);
		}
	}
}

/************************************************************************/
/* calcluate the size of prefixsize                                     */
/************************************************************************/
void CMotifInfo::initPrefixsize()
{
	m_prefixsize[0] = 1;//the number of prefix of lenght 0 is 1
	m_prefixsizestart[0] = 0;
	m_prefixtotalsize = m_prefixsize[0];
	for(int i=1; i<m_motiflen+1; i++) {
		//calculate m_prefixsize[i], that is the number of prefix of length i
		if (m_motif[i-1]>ALPHA_NUM-1) {//need to extend
			if (m_motif[i-1]==CAlphabet::N) {
				m_prefixsize[i] = m_prefixsize[i-1]<<2;
			} else 
				m_prefixsize[i] = m_prefixsize[i-1]<<1;			
		} else
			m_prefixsize[i] = m_prefixsize[i-1];
		m_prefixsizestart[i] = m_prefixtotalsize;
		m_prefixtotalsize += m_prefixsize[i];
	}	

}

void CMotifInfo::printprefix()
{
	for(int j=0; j<MOTIFSETS+1; j++) {
		for(int i=0; i<m_multiprefixtotalsize[j]; i++) {
			printf("%d  ",i);
			multiprefixset[j][i]->printprefix();
			if (m_fullhit[j]!=NULL) {
				printf("  %d",m_fullhit[j][i]);
			}
			printf("\n");
		}
	}	
}
/************************************************************************/
/* index-- the index of old prefix in prefixset                         */
/* add  -- the char append to the old prefix                            */
/* return value: the index of new prefix after append                   */
/************************************************************************/
int CMotifInfo::findprefix(int index, CPrefix* oldprefix, CAlphabet::ALPHABET add)
{
	if (oldprefix==NULL) {
		return -1;
	}
	if (oldprefix->length<m_motiflen) {//not the whole motif
		CAlphabet::ALPHABET next = m_motif[oldprefix->length];
		if (p_alpha->compatible[next][add]) {//compatiable
			return m_prefixsize[oldprefix->length]*(p_alpha->compatibleindex[next][add]+1) + index;
		}
	}
	//generate a new prefix
	CPrefix* newprefix = NULL;
	if (add!=CAlphabet::ERROR) {
		newprefix = new CPrefix(oldprefix->length+1);
		newprefix->append(oldprefix->m_prefix,add);
	} else {
		newprefix = new CPrefix(oldprefix->length-1);
		//todo
		newprefix->removeFirst(oldprefix);
	}
	
	/*
	 * check whether all char in newprefix of region[startpos, endpos] 
	 * is compatiable to orginal prefix	
	 */
	int findlen = 0;//if prestart=newprefix->length, findlen = 0
	int prestart = 0;
	while (prestart<newprefix->length) {
		int preindex = prestart;
		int orgindex = 0;
		while (preindex<newprefix->length && orgindex<m_motiflen && p_alpha->compatible[m_motif[orgindex]][newprefix->m_prefix[preindex]]) {
			preindex++;
			orgindex++;
		}
		if (preindex==newprefix->length) {//find it
			findlen = newprefix->length-prestart;
			break;
		} else {//should cut one char in start
			prestart++;
		}
	}
	int ret = findprefixoflen(newprefix,findlen);
	delete newprefix;
	return ret;
}

int CMotifInfo::findprefixoflen(CPrefix *newprefix, int len)
{
	if (len==0) {
		return 0;
	} else {
		int index = m_prefixsizestart[len];
		while (index<m_prefixsizestart[len+1]) {//search all prefix of length len
			//compare newprefix and prefix in prefixset
			int newindex = newprefix->length-len;
			int preindex = 0;
			while (newindex<newprefix->length && preindex<len && newprefix->m_prefix[newindex]==prefixset[index]->m_prefix[preindex]) {
				newindex++;
				preindex++;
			}
			if (newindex==newprefix->length) {//find it
				return index;
			} else 
				index++;
		}
		return 0;
	}
}

void CMotifInfo::printskiplist()
{
	if (type==2) {
		for(int k=0; k<MOTIFSETS+1; k++) {
			for(int i=0; i<ALPHA_NUM; i++) {
				for(int j=0; j<m_multiprefixtotalsize[k]; j++) {
					multiprefixset[k][j]->printprefix();
					printf(" + ");
					p_alpha->printalp((CAlphabet::ALPHABET)i); 
					printf(" ==> ");
					multiprefixset[k][multiskiplist[k][i][j]]->printprefix();
					printf("\n");
					for(int s=0; s<MOTIFSETS; s++) {
						printf("cut in %d \n", s);
						for(int b=0; b<multiprefixset[k][j]->length+1; b++) {
							printf(" boundary %d: %d, ", b, m_cutnum[k][i][s][j][b]);
						}
						printf("\n");
					}
					printf("\n");
				}
			}
		}
	} else {
		for(int i=0; i<ALPHA_NUM; i++) {
			for(int j=0; j<m_prefixtotalsize; j++) {
				prefixset[j]->printprefix();
				if (prefixset[j]->length<m_motiflen) {
					printf(" + ");
					p_alpha->printalp((CAlphabet::ALPHABET)i);
				} else {
					printf(" cut first ");
				}
				printf(" ==> ");
				prefixset[skiplist[i][j]]->printprefix();
				printf("\n");
			}
		}
	}	
}

void CMotifInfo::initMultiMotifs(int index)
{
	if (type==1) {
		if (m_multinum>0 && m_multimotifs==NULL) {
			m_multimotifs = (CAlphabet::ALPHABET**)malloc(sizeof(CAlphabet::ALPHABET*)*m_multinum);
			for(int i=0; i<m_multinum; i++)// each motifs with the same length
				m_multimotifs[i] = new CAlphabet::ALPHABET[m_motiflen];
		}
	} else {//init several sets
		if (m_motifsetsnum[index]>0 && m_motifsets[index]==NULL) {
			m_motifsets[index] = (CAlphabet::ALPHABET**)malloc(sizeof(CAlphabet::ALPHABET*)*m_motifsetsnum[index]);
			for(int k=0; k<m_motifsetsnum[index]; k++) {
				m_motifsets[index][k] = new CAlphabet::ALPHABET[m_motifsetslen[index]];
			}
		}		
	}
}

void CMotifInfo::printmultimotifs()
{
	int i,j,k;
	if (type==1) {
		printf(m_motifname);
		printf("\n");
		for(i=0; i<m_multinum; i++) {
			printf("\t");
			for(j=0; j<m_motiflen; j++) {
				p_alpha->printalp(m_multimotifs[i][j]);
			}
			printf("\n");
		}
	} else {
		printf(m_motifname);
		printf("\n");
		for(k=0; k<MOTIFSETS; k++) {
			for(i=0; i<m_motifsetsnum[k]; i++) {
				printf("\t");
				for(j=0; j<m_motifsetslen[k]; j++) {
					p_alpha->printalp(m_motifsets[k][i][j]);
				}
				printf("\n");
			}
			printf("\n\n");
		}
	}
}

// m_motiflen should be set to the length of CAlphabet::ALPHABET* sequences
// m_multinum should be set to the number of sequences
// totallength is the length of prefix
// we init startpos, prefix according to sequences and return the real length of prefix
int CMotifInfo::initMultiPrefixset(int totallength, CPrefix** prefix, int* startpos, CAlphabet::ALPHABET** sequences, int* hits,int setindex)
{
	//allocate memory for prefixsetE
	prefix[0] = new CPrefix(0);//epsilon
	if (hits!=NULL) {
		hits[0] = -1;
	}
	int prefixsize = 1;
	startpos[0] = 0;
	int startindex = 0;
	for(int i=1; i<m_motiflen+1; i++) {
		startindex = prefixsize;
		startpos[i] = startindex;
		for(int j=0; j<m_multinum; j++) {
			//put the prefix of m_multimotifs[j] with length i into prefixset
			int ch = startindex;
			while (ch<prefixsize && prefix[ch]->isSame(sequences[j],0,i)==false) {
				ch++;
			}
			if (ch==prefixsize) {//should add it
				prefix[prefixsize] = new CPrefix(i);	
				prefix[prefixsize]->duplicate(sequences[j],i);
				if (hits!=NULL) {
					if (i==m_motiflen) {
						hits[prefixsize] = setindex;
					} else {
						hits[prefixsize] = -1;
					}
				}
				prefixsize++;
			}
		}
	}
	return prefixsize;
}

void CMotifInfo::initMultiSkiplist(int* skip[], CPrefix** prefix, int prefixsize)
{
	int i,j;
	for(i=0; i<ALPHA_NUM; i++) {
		for(j=0; j<prefixsize; j++) {
			if (prefix[j]->length<m_motiflen) {
				skip[i][j] = findmultiprefix(prefix, prefix[j], (CAlphabet::ALPHABET)i, prefixsize,true);
			} else 
				skip[i][j] = findmultiprefix(prefix, prefix[j], CAlphabet::ERROR, prefixsize,true);
		}
	}
}


// find the new prefix index constructed from oldprefix in prefixes
// newprefix = oldprefix+a
// prefixsize is the size of prefixes
// m_motiflen should be set to be the longest length of motifs in prefixes
// m_prefixstartsize should be set to be the start index of each length prefix int prefixes

int CMotifInfo::findmultiprefix(CPrefix** prefixes, CPrefix* oldprefix, CAlphabet::ALPHABET add, int prefixsize, bool cut)
{
	if (oldprefix==NULL) {
		return -1;
	}
	//generate a new prefix
	CPrefix* newprefix = NULL;
	if (add!=CAlphabet::ERROR) {
		newprefix = new CPrefix(oldprefix->length+1);
		newprefix->append(oldprefix->m_prefix,add);
	} else {//the full motif
		if (cut) {
			newprefix = new CPrefix(oldprefix->length-1);
			newprefix->removeFirst(oldprefix);
		} else {
			newprefix = new CPrefix(oldprefix->length);
			newprefix->duplicate(oldprefix->m_prefix,oldprefix->length);
		}
	}
	
	/*
	 * check whether all char in newprefix of region[prestart, end] 
	 * is the same to same prefix in the prefixset
	 */
	int prestart = 0;
	int indexofprefixset = 0;
	int end = 0;
	if (newprefix->length-prestart>m_motiflen) {
		prestart = newprefix->length-m_motiflen;
	}
	while (prestart<newprefix->length) {
		indexofprefixset = m_prefixsizestart[newprefix->length-prestart];
		if (newprefix->length-prestart == m_motiflen) {
			end = prefixsize;
		} else 
			end = m_prefixsizestart[newprefix->length-prestart+1];
		while (indexofprefixset<end && prefixes[indexofprefixset]->isSame(newprefix->m_prefix,prestart,newprefix->length)==false) {
			indexofprefixset++;
		}
		if (indexofprefixset==end) {// not find
			prestart++;
		} else {//find it
			delete newprefix;
			return indexofprefixset;
		}
	}
	delete newprefix;
	return 0;
}

void CMotifInfo::preProcessMulti()
{
	m_prefixtotalsize = m_motiflen*m_multinum+1;//it is a upper bound, it's possiable that total size is less than this value
	prefixset = (CPrefix**) malloc(m_prefixtotalsize*sizeof(CPrefix*));
	m_prefixsizestart = new int[m_motiflen+1];
	m_prefixtotalsize = initMultiPrefixset(m_prefixtotalsize, prefixset, m_prefixsizestart,m_multimotifs, NULL,0);
	for(int i=0; i<ALPHA_NUM; i++) {
		skiplist[i] = new int[m_prefixtotalsize];
	}
	initMultiSkiplist(skiplist, prefixset, m_prefixtotalsize);
}

void CMotifInfo::preProcessMotifSets()
{
	int i,j;
	int* prefixsizestart[MOTIFSETS];
	for(i=0; i<MOTIFSETS; i++) {
		m_motiflen = m_motifsetslen[i];
		m_multinum = m_motifsetsnum[i];
		m_multiprefixtotalsize[i] = m_motiflen*m_multinum+1;
		multiprefixset[i] = (CPrefix**) malloc(m_multiprefixtotalsize[i]*sizeof(CPrefix*));
		prefixsizestart[i] = new int[m_motiflen+1];
		m_fullhit[i] = new int[m_multiprefixtotalsize[i]];
		m_multiprefixtotalsize[i] = initMultiPrefixset(m_multiprefixtotalsize[i],multiprefixset[i],prefixsizestart[i],m_motifsets[i], m_fullhit[i], i);		
		for(j=0; j<ALPHA_NUM; j++) {
			multiskiplist[i][j] = new int[m_multiprefixtotalsize[i]];
		}
		m_prefixsizestart = prefixsizestart[i];
		initMultiSkiplist(multiskiplist[i],multiprefixset[i],m_multiprefixtotalsize[i]);
		cutNumber(i);
	}
	//combine the prefix to prefixset
	int maxlen = 0;
	m_multiprefixtotalsize[MOTIFSETS] = 0;
	for(i=0; i<MOTIFSETS; i++) {//find the longest
		if (m_motifsetslen[i]>maxlen) {
			maxlen = m_motifsetslen[i];
		}
		m_multiprefixtotalsize[MOTIFSETS]+= m_multiprefixtotalsize[i];
	}
	m_motiflen  = maxlen;
	m_prefixsizestart = new int[m_motiflen+1];
	m_fullhit[MOTIFSETS] = new int[m_multiprefixtotalsize[MOTIFSETS]];
	multiprefixset[MOTIFSETS] = (CPrefix**) malloc(m_multiprefixtotalsize[MOTIFSETS]*sizeof(CPrefix*));//it is a upper bound
	m_multiprefixtotalsize[MOTIFSETS] = Combine(m_multiprefixtotalsize, multiprefixset, prefixsizestart, multiprefixset[MOTIFSETS]);
	for(j=0; j<ALPHA_NUM; j++) {
		multiskiplist[MOTIFSETS][j] = new int[m_multiprefixtotalsize[MOTIFSETS]];
	}
	initMultiSkiplist(multiskiplist[MOTIFSETS],multiprefixset[MOTIFSETS],m_multiprefixtotalsize[MOTIFSETS]);
	cutNumber(MOTIFSETS);
	////////////////////////////init degskiplist/////////////////////////////////////////
	int* temp = m_prefixsizestart;
	m_prefixtotalsize = m_multiprefixtotalsize[MOTIFSETS];
	for(int motifindex = 0; motifindex<MOTIFSETS; motifindex++) {
		degskiplist[motifindex] = new int[m_multiprefixtotalsize[MOTIFSETS]];
		m_motiflen = m_motifsetslen[motifindex];
		m_prefixsizestart = prefixsizestart[motifindex];
		for(int prefixindex=0; prefixindex<m_multiprefixtotalsize[MOTIFSETS]; prefixindex++) {
			degskiplist[motifindex][prefixindex] = findmultiprefix(multiprefixset[motifindex],multiprefixset[MOTIFSETS][prefixindex],CAlphabet::ERROR,m_multiprefixtotalsize[motifindex],false);
		}
	}
	m_motiflen = maxlen;
	degcunNumber();
	/////////////////////////////realase memeory/////////////////////////
	for(i=0; i<MOTIFSETS; i++)
		delete prefixsizestart[i];
	delete temp;
	m_prefixsizestart = NULL;
//	printprefix();
//	printskiplist();
//	printfdegskiplist();
} 

/************************************************************************/
/* combine the all prefix into prefixset and return the real number of 
prefixset                                                               */
/************************************************************************/
int CMotifInfo::Combine(int *prefixtotalsize, CPrefix **prefix[], int *prefixsizestart[], CPrefix** comprefix)
{
	int prefixsize = 0;
	int currlen = 0;
	int end = 0;
	int ch = 0;
	while (currlen<m_motiflen+1) {
		//insert prefix with length currlen
		m_prefixsizestart[currlen] = prefixsize;
		for(int j=0; j<MOTIFSETS; j++) {
			if (currlen>m_motifsetslen[j]) {
				continue;
			}
			if (currlen+1<m_motifsetslen[j]) {
				end = prefixsizestart[j][currlen+1];
			} else
				end = prefixtotalsize[j];
			for(int k = prefixsizestart[j][currlen]; k<end; k++) {
				//add prefix[j][k] to prefixset
				ch = m_prefixsizestart[currlen];
				while (ch<prefixsize && comprefix[ch]->isSame(prefix[j][k]->m_prefix,0,currlen)==false) {
					ch++;
				}
				if (ch==prefixsize) {//should add it
					comprefix[prefixsize] = new CPrefix(currlen);
					comprefix[prefixsize]->duplicate(prefix[j][k]->m_prefix,currlen);
					if (currlen==m_motiflen) {
						m_fullhit[MOTIFSETS][prefixsize] = 1;
					} else 
						m_fullhit[MOTIFSETS][prefixsize] = -1;
					prefixsize++;
				} 	
			}
		}
		currlen++;
	}
	return prefixsize;
}

int CMotifInfo::initCutNumber(int setindex, int addchar, int multimoitfindex, int prefixindex, int boundary)
{
	//the oldprefix is
	CPrefix* oldprefix = multiprefixset[setindex][prefixindex];
	CPrefix* newprefix;
	int cutlength = 0;
	if (oldprefix->length==m_motiflen) {
		cutlength = oldprefix->length - multiprefixset[setindex][multiskiplist[setindex][addchar][prefixindex]]->length;
		newprefix = new CPrefix(oldprefix->length);
		newprefix->duplicate(oldprefix->m_prefix,oldprefix->length);
	} else {
		cutlength = oldprefix->length+1 - multiprefixset[setindex][multiskiplist[setindex][addchar][prefixindex]]->length;
		newprefix = new CPrefix(oldprefix->length+1);
		newprefix->append(oldprefix->m_prefix,(CAlphabet::ALPHABET)addchar);
	}
	
	if (cutlength==0) {
		delete newprefix;
		return 0;
	}
	int startpos = 0;
	int endpos = startpos+m_motifsetslen[multimoitfindex]; 
	int j;
	int cutnum = 0;
	while (startpos<cutlength && !(endpos>newprefix->length)) {
		if (!(startpos<boundary && endpos>boundary)) {
			//it is cross the boundary, no cut
			j=0;
			while (j<m_motifsetsnum[multimoitfindex] && newprefix->isSame(m_motifsets[multimoitfindex][j],startpos, endpos, 0, m_motifsetslen[multimoitfindex])==false ) {
				j++;
			}
			if (j<m_motifsetsnum[multimoitfindex]) {
				cutnum++;
			}			
		}
		startpos++;
		endpos = startpos+m_motifsetslen[multimoitfindex];
	}
	delete newprefix;
	return cutnum;
}

void CMotifInfo::cutNumber(int setindex)
{
	for(int addchar=0; addchar<ALPHA_NUM; addchar++) {
		for(int motifindex=0; motifindex<MOTIFSETS; motifindex++) {
			m_cutnum[setindex][addchar][motifindex] = (int**)malloc(sizeof(int*)*m_multiprefixtotalsize[setindex]);
			for(int prefixindex=0; prefixindex<m_multiprefixtotalsize[setindex]; prefixindex++) {
				m_cutnum[setindex][addchar][motifindex][prefixindex] = new int[multiprefixset[setindex][prefixindex]->length+1];
				for(int boundary=0; boundary<multiprefixset[setindex][prefixindex]->length+1; boundary++) {
					m_cutnum[setindex][addchar][motifindex][prefixindex][boundary] = 
					initCutNumber(setindex,addchar,motifindex,prefixindex,boundary);
				}
			}
		}
	}
}


void CMotifInfo::degcunNumber()
{
	int startpos, endpos,cutlength,cutnum,j;
	for(int motifindex=0; motifindex<MOTIFSETS; motifindex++) {
		m_degcutnum[motifindex] = (int**)malloc(sizeof(int*)*m_multiprefixtotalsize[MOTIFSETS]);
		for(int prefixindex = 0; prefixindex<m_multiprefixtotalsize[MOTIFSETS]; prefixindex++) {
			m_degcutnum[motifindex][prefixindex] = new int[multiprefixset[MOTIFSETS][prefixindex]->length+1];
			for(int boundary =0; boundary<multiprefixset[MOTIFSETS][prefixindex]->length+1; boundary++) {
				CPrefix* oldprefix = multiprefixset[MOTIFSETS][prefixindex];
				CPrefix* newprefix = multiprefixset[motifindex][degskiplist[motifindex][prefixindex]];
				cutlength = oldprefix->length-newprefix->length;
				if (cutlength==0) {
					m_degcutnum[motifindex][prefixindex][boundary] = 0;
				} else {
					cutnum = 0;
					startpos = 0;
					endpos = startpos+m_motifsetslen[motifindex];
					while (!(endpos>m_motiflen)&& startpos<cutlength) {
						if(!(startpos<boundary && endpos>boundary)) {
							j=0;
							while (j<m_motifsetsnum[motifindex] && oldprefix->isSame(m_motifsets[motifindex][j],startpos, endpos, 0, m_motifsetslen[motifindex])==false ) {
								j++;
							}
							if (j<m_motifsetsnum[motifindex]) {
								cutnum++;
							}
						}
						startpos++;
						endpos = startpos+m_motifsetslen[motifindex];
					}
					m_degcutnum[motifindex][prefixindex][boundary] = cutnum;
				}
			}
		}
	}
}


void CMotifInfo::printfdegskiplist()
{
	for(int i=0; i<m_multiprefixtotalsize[MOTIFSETS]; i++) {
		multiprefixset[MOTIFSETS][i]->printprefix();
		printf("\n");
		for(int j= 0; j<MOTIFSETS; j++) {
			printf(" to set %d :", j);
			multiprefixset[j][degskiplist[j][i]]->printprefix();
			for(int k=0; k<multiprefixset[MOTIFSETS][i]->length+1; k++)
				printf(" boundary %d, cutnumber %d", k, m_degcutnum[j][i][k]);
			printf("\n");
		}
		printf("\n");
	}
}
