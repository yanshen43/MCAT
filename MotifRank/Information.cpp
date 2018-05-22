                                                                             // Information.cpp: implementation of the CInformation class.
//
//////////////////////////////////////////////////////////////////////

#include "Information.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CInformation::CInformation()
{
	motifset = NULL;
	region = new CRegionInfo();
	result = new CResultInfo();
	alp = new CAlphabet();
}

CInformation::~CInformation()
{
	if (motifset!=NULL) {
		delete motifset;
	}
	if (region!=NULL) {
		delete region;
	}
	if (result!=NULL) {
		delete result;
	}
	if (alp!=NULL) {
		delete alp;
	}
}

/************************************************************************/
/* read a double number from hile, ignore blank line                    */
/* if read a double number successfully, readflag will be true          */
/* if read the end of file, endflag will be true                        */
/************************************************************************/
double CInformation::readDouble(FILE *hfile, bool* readflag, bool* endflag)
{
	char one[1];
	char number[64];
	memset(number,0,sizeof(char)*64);
	int size=0;
	while(size<64 && fread(one,sizeof(char),1,hfile)>0) {
		if ((*one == '\r') || (*one == '\n')) {//read the end
			if (size==0) {//blank line
				*readflag= false;
				*endflag = false;
				return 0;
			} else
				break;
		} else {
			number[size++] = *one;
		}	
	}
	if (size==0) {
		*endflag = true;
		return 0;
	}
	*endflag = false;
	*readflag = true;
	char* endptr;
	double ret = strtod(number, &endptr);
	if(ret==0) {
		*endflag=true;
	} 
	return ret;
}

/************************************************************************/
/* read a integer from hile, ignore blank line                          */
/* if read integer successfully, readflag will be true                  */
/* if read the end of file, endflag will be true                        */
/************************************************************************/
int CInformation::readInt(FILE *hfile, bool* readflag, bool* endflag)
{
	char one[1];
	char number[64];
	memset(number,0,sizeof(char)*64);
	int size=0;
	while(size<64) {
		fread(one,sizeof(char),1,hfile);
		if ((*one == '\r') || (*one == '\n')) {//read the end
			if (size==0) {//blank line
				*readflag = false;
				*endflag  = false;
				return 0;
			} else
				break;
		} else {
			number[size++] = *one;
		}	
	}
	if (size==0) {
		*endflag = true;
		return 0;
	}
	*endflag = false;
	*readflag = true;
	return atoi(number);
}

/************************************************************************/
/* read a char from hile, ignore blank line                             */
/* if read a char successfully, readflag will be true                   */
/* if read the end of file, endflag will be true                        */
/************************************************************************/
char CInformation::readChar(FILE *hfile, bool *readflag, bool *endflag)
{
	char one[1];
	while (true) {
		if (fread(one,sizeof(char),1,hfile)==0) {
			*endflag = true;//the end of the file
			return '0';
		} else if ( (*one>'a'-1 && *one<'z'+1) || (*one>'A'-1 && *one<'Z'+1)) {
			*readflag = true;
			*endflag = false;
			return *one;
		}
	}
}

void CInformation::readLine(FILE *hfile, bool *readflag, bool *endflag, int *len, char* oneline)
{
	char one[1];
	*len = 0;
	while (true) {
		if(fread(one,sizeof(char),1,hfile)==0) {
			*endflag = true;
			*readflag = false;
			return;
		} else if ((*one == '\r') || (*one == '\n')) {
			if (*len==0) {
				continue;
			} else {
				*readflag = true;
				*endflag = false;
				return;
			}
		} else {//if ((*one>'a'-1 && *one<'z'+1) || (*one>'A'-1 && *one<'Z'+1)) {
			if(*len>MAX_LEN-1) {
				*readflag = false;
				*endflag = true;
				return;
			}
			oneline[(*len)++] = one[0];
		} 
	}
}


/************************************************************************/
/* process the motif information                                        */
/************************************************************************/
void CInformation::preProcess()
{
	if (motifset!=NULL && motifset->numberofmotif>0) {
		int maxlen = 0;
		int maxmotif = 0;
		int maxhit = 0;
		int maxmultihit[MOTIFSETS];
		memset(maxmultihit,0,sizeof(int)*MOTIFSETS);
		for(int i=0; i<motifset->numberofmotif; i++) {
			if (motifset->motifs[i]->fdo==false) {
				continue;
			}
			motifset->motifs[i]->preProcess();
			if (motifset->motifs[i]->m_prefixtotalsize>maxlen) {
				maxlen = motifset->motifs[i]->m_prefixtotalsize;
			}
			if (motifset->motifs[i]->m_motiflen>maxmotif) {
				maxmotif = motifset->motifs[i]->m_motiflen;
			}
			if (motifset->motifs[i]->hittimes>maxhit) {
				maxhit = motifset->motifs[i]->hittimes;
			}
			for(int j=0; j<MOTIFSETS; j++) {
				if (motifset->motifs[i]->motifsetshittimes[j]>maxmultihit[j]) {
					maxmultihit[j] = motifset->motifs[i]->motifsetshittimes[j];
				}
			}
		}
		//changed
		if (motifset->type<2) {
			result->inittotal(maxmotif, maxlen, maxhit);
		} else 
			result->inittotal(maxmotif, maxlen, maxmultihit);
	
	}
}

/************************************************************************/
/* init the promoter sequences                                          */
/*
 *	@param hfile: the FILE handle to the fasta file
 *  @return int:
 *  0	read sequences succussfully 
 *  -1	otherwise
 */
/************************************************************************/
int CInformation::initSequences(FILE *hfile)
{
	if (hfile==NULL) {
		return -1;
	}
	char line[MAXREGION];
	memset(line, 0, sizeof(char)*MAXREGION);
	int num = 0;
	while (readLine(hfile, line, MAXREGION, '\n',0) != -1) {//the decription of promoter 
		memset(line, 0, sizeof(char)*MAXREGION);
		num = readLine(hfile,line,MAXREGION, '>',1);
		while (num==0) {//read until the end or read something
			num = readLine(hfile,line,MAXREGION, '>',1);
		}
		if (num==-1) {//the end and nothing read
			return -1;
		} else if (region->m_regionlen==0) {//init the region length
			region->m_regionlen = num;
		} else if (region->m_regionlen!=num) {//if the sequences are not the same length
			return -1;
		}
		region->sequences[region->m_sequences] = new char[region->m_regionlen+1];
		memset(region->sequences[region->m_sequences],0,sizeof(char)*(region->m_regionlen+1));
		memcpy(region->sequences[region->m_sequences], line,sizeof(char)*region->m_regionlen);//copy from the read line to sequences
		region->m_sequences++;
		if (region->m_sequences==MAXSEQUENCE) {//exceed the max sequnces num
			return -1;
		}
	}
	if (region->m_sequences<1) {
		return -1;
	//} else if (region->m_sequences>1) {
	//	return -2;
	} else {
		if (rc==true) {
			region->initRCpromoter();
		}
		return 0;
	}
}
/************************************************************************/
/* init motifs                                                          */
/*
 *	@param hfile: the FILE handle to the motif file
 ***************************file format***********************************
 line 1: type of motifs
 line 2: the number of motifs
 
 TYPE 0:
 @ID MCB
 ACWYSNNCA
 
 TYPE 1
 @ID MCB
 @SIZE 1
 ATGCTAC

 TYPE 2
 @ID MOTIF1
 @SIZE 1 1
 ATGCTAC
 CAAATTA
 
 TYPE 3
 @ID MCB
 A	0.833333	0.000000	0.000000	0.000000	0.000000	0.250000
 C	0.083333	1.000000	0.000000	1.000000	0.000000	0.000000
 G	0.000000	0.000000	1.000000	0.000000	1.000000	0.000000
 T	0.083333	0.000000	0.000000	0.000000	0.000000	0.750000

 TYPE 4
 @ID MOTIF3
 A	0.833333	0.000000	0.000000	0.000000	0.000000	0.250000
 C	0.083333	1.000000	0.000000	1.000000	0.000000	0.000000
 G	0.000000	0.000000	1.000000	0.000000	1.000000	0.000000
 T	0.083333	0.000000	0.000000	0.000000	0.000000	0.750000
 
 A	0.833333	0.000000	0.000000	0.000000	0.000000	0.250000  0.500000
 C	0.083333	1.000000	0.000000	1.000000	0.000000	0.000000  0.500000
 G	0.000000	0.000000	1.000000	0.000000	1.000000	0.000000  0.000000
 T	0.083333	0.000000	0.000000	0.000000	0.000000	0.750000  0.000000
 
 *************************************************************************
 *  @return int:
 *  0	read motifs succussfully 
 *  -1	otherwise
 */
/************************************************************************/
int CInformation::initMotif(FILE *hfile)
{
	if (hfile==NULL) {
		return -1;
	}
	/////////read the number of motifs to rank and init motifset///////////
	bool readflag = false;
	bool endflag = true;
	int size;
	while (readflag==false) {
		size = readInt(hfile,&readflag,&endflag);
		if (endflag==true) {//the end of file, invalid file format
			return -1;
		}
	}
	if (size<1) {//number of motifs should greater than 0
		return -1;
	} else {
		motifset  = new CMotifSet(size);
	}
	/////////////read the type of motifs///////////////////////////////////
	readflag = false;
	endflag = true;
	while (readflag == false) {
		motifset->type = readInt(hfile,&readflag,&endflag);
		if (readflag==false && endflag == true) {
			return false;
		}
	}
	///////////////read motifs according to the type///////////////////////
	switch(motifset->type) {
	case 0://a motif word on IUPAC
		if (readMotifCase0(hfile)==-1) {
			return -1;
		}
		break;
	case 1://a set of motif words
		if (readMotifCase1(hfile)==-1) {
			return -1;
		}
		break;
	case 2:
		if (readMotifCase2(hfile)==-1) {
			return -1;
		}
		break;
	case 3:
		motifset->type = 1;
		if (readMotifCase3(hfile)==-1) {
			return -1;
		}
		break;
	case 4:
		motifset->type = 2;
		if (readMotifCase4(hfile)==-1) {
			return -1;
		}
		break;
	default:
		return -1;
	}
	return 0;
}

int CInformation::initMatrix(FILE *hfile)
{
	////////////////////4.read the init probability//////////////////////
	int i,j;
	bool readflag, endflag;
	for(i=0; i<ALPHA_NUM; i++) {
		readflag = false;
		endflag = true;
		while (readflag==false) {
			region->initprobability[i] = readDouble(hfile,&readflag, &endflag);//init initprobability[i][j]
			if (endflag==true) {
				return -1;
			}
		}
		if (region->initprobability[i]<0 || region->initprobability[i]==0) {//check validation
			return -1;
		}
	}
	//////////////////5.read the transition probability///////////////////
	for(i=0;i<ALPHA_NUM; i++) {
		for(j=0; j<ALPHA_NUM;j++) {
			readflag = false;
			endflag = true;
			while (readflag==false) {
				region->transitionmatrix[i][j] = readDouble(hfile, &readflag, &endflag);//init transitionmatrix[i][j][k]
				if (endflag==true) {
					return -1;
				}
			}
			if (region->transitionmatrix[i][j]<0 || region->transitionmatrix[i]==0) {//check validation
				return -1;
			}
		}
	}
	return 0;
}

/************************************************************************/
/* read from hfile until the endtoken                                   */
/*
 *	@param hfile: the file handle
 *  @param line: store the string read
 *  @param maxlinelen: the max size of line
 *  @param endtoken: read from hfile unitl endtoken or EOF
 *  @param readtype: the type of read
                     -1: filter all
                     0: read all
					 1: check {A,C,G,T}
					 2: check {A,C,G,T,R,W,S,Y,N}
					 3: leave number and alphabet
					 4: leave number
 *  @return int: the size of char read
 *   -1 file end and read nothing or there are some error
 *    0 read nothing
 */
/************************************************************************/
int CInformation::readLine(FILE *hfile, char line[], int maxlinelen, char endtoken, int readtype)
{
	if (hfile==0 || (readtype!=-1 && line==NULL)) {
		return -1;
	}
	char one[1];
	int numread = 0;
	if (readtype==-1) {
		while (fread(one,sizeof(char),1,hfile)!=0) {
			if (*one==endtoken) {//the end, not add the char into the buf
				return 1;
			}
		}
	} else if (readtype==0) {
		while (fread(one,sizeof(char),1,hfile)!=0) {
			if (*one==endtoken) {//the end, not add the char into the buf
				return numread;
			} else  {//add all
				line[numread++] = *one;
				if (numread>maxlinelen) {
					return -1;
				}
			}
		}
	} else if (readtype==1 || readtype==2) {
		while (fread(one,sizeof(char),1,hfile)!=0) {
			if (*one==endtoken) {//the end, not add the char into the buf
				return numread;
			} else if(alp->checkValid(readtype,one[0])) {//{A,C,G,T} or {A,C,G,T,R,W,Y,S,N}
				line[numread++] = *one;
				if (numread>maxlinelen) {
					return -1;
				}
			}
		}
	} else if(readtype==3) {
		while (fread(one,sizeof(char),1,hfile)!=0) {
			if (*one==endtoken) {//the end, not add the char into the buf
				return numread;
			} else if((*one>'A'-1 && *one<'Z'+1) || (*one>'a'-1 && *one<'z'+1) || (*one>'0'-1 && *one<'9'+1)) {//{A,C,G,T}
				line[numread++] = *one;
				if (numread>maxlinelen) {
					return -1;
				}
			}
		}
	} else {
		while (fread(one,sizeof(char),1,hfile)!=0) {
			if (*one==endtoken) {//the end, not add the char into the buf
				return numread;
			} else if( *one>'0'-1 && *one<'9'+1) {//should be number
				line[numread++] = *one;
				if (numread>maxlinelen) {
					return -1;
				}
			}
		}
	}
	if (numread==0) {
		return -1;
	} else 
		return numread;
}



/************************************************************************/
/* TYPE=0, all motifs are deterministic motifs based on based alphabet  */
/* or IUPAC alphabet                             						*/
/************************************************************************/
int CInformation::readMotifCase0(FILE* hfile)
{
	//read until @
	if (readLine(hfile,NULL,0,'@',-1)==-1) {
		return -1;
	}
	int index = 0;
	while(index<motifset->numberofmotif) {
		//first readname
		if (readName(hfile,index)==-1) {
			return -1;
		}
		//then read motifs
		char oneline[MAX_LEN];
		memset(oneline,0,sizeof(char)*MAX_LEN);
		if ((motifset->motifs[index]->m_motiflen = readLine(hfile,oneline,MAX_LEN,'@',2))<1) {//read the motif unitl next @
			return -1;
		}
		motifset->motifs[index]->initMotif();//allocate memory for motif with index
		motifset->motifs[index]->type = 0;
		for(int j=0; j<motifset->motifs[index]->m_motiflen;j++){
			motifset->motifs[index]->insert(oneline[j],-1,-1); 
		}
		motifset->motifs[index]->hittimes = calHit(motifset->motifs[index]->m_motif,motifset->motifs[index]->m_motiflen);//calculate the hit times
		index++;
	}
	return 0;
}

/************************************************************************/
/* TYPE=1, the moitif is a set of words on basic alphabet               */
/************************************************************************/
int CInformation::readMotifCase1(FILE *hfile)
{
	int index = 0;
	int ret = 0;
	int length = 0;
	char buffer[MAX_LEN];
	while (index<motifset->numberofmotif) {
		//read until @
		if (readLine(hfile,NULL,0,'@',-1)==-1) {
			return -1;
		}
		//first read name
		if (readName(hfile,index)==-1) {
			return -1;
		}
		//then read size of set
		if (readSize(hfile,index,0)==-1) {
			return -1;
		}
		length = 0;
		for(int i=0; i<motifset->motifs[index]->m_multinum; i++) {
			memset(buffer,0,sizeof(char)*MAX_LEN);
			ret = readLine(hfile,buffer,MAX_LEN,'\n',1);
			while (ret==0) {
				ret = readLine(hfile,buffer,MAX_LEN,'\n',1);
			}
			if (ret==-1) {
				return -1;
			} else {
				if (length==0) {
					length = ret;
					motifset->motifs[index]->m_motiflen = ret;
					motifset->motifs[index]->type = 1;
					motifset->motifs[index]->initMultiMotifs(-1);
				} else if (length!=ret) {//different length, wrong
					return -1;
				} 	
				motifset->motifs[index]->m_index = 0;
				for(int k=0; k<motifset->motifs[index]->m_motiflen;k++) {//copy from buffer to motif
					motifset->motifs[index]->insert(buffer[k],i,-1);
				}
			}
		}
		motifset->motifs[index]->hittimes = calHit(motifset->motifs[index]->m_multimotifs,motifset->motifs[index]->m_multinum,motifset->motifs[index]->m_motiflen);
		index++;
	}
	for(int j=0; j<motifset->numberofmotif; j++)
		motifset->motifs[j]->printmultimotifs();
	return 0;
}
/************************************************************************/
/* TYPE=2, the moitif is two sets of words on basic alphabet            */
/************************************************************************/
int CInformation::readMotifCase2(FILE *hfile)
{
	int index = 0;
	int ret = 0;
	int length = 0;
	char buffer[MAX_LEN];
	while (index<motifset->numberofmotif) {
		//read until @
		if (readLine(hfile,NULL,0,'@',-1)==-1) {
			return -1;
		}
		//first read name
		if (readName(hfile,index)==-1) {
			return -1;
		}
		for(int s=0; s<MOTIFSETS; s++) {
			if (readSize(hfile,index,s)==-1) {//then read size of set
				return -1;
			}
			length = 0;
			for(int i=0; i<motifset->motifs[index]->m_motifsetsnum[s]; i++) {//read motifs
				memset(buffer,0,sizeof(char)*MAX_LEN);
				ret = readLine(hfile,buffer,MAX_LEN,'\n',1);
				while (ret==0) {
					ret = readLine(hfile,buffer,MAX_LEN,'\n',1);
				}
				if (ret==-1) {
					return -1;
				} else {
					if (length==0) {//the first
						length = ret;
						motifset->motifs[index]->m_motifsetslen[s] = ret;
						motifset->motifs[index]->type = 2;
						motifset->motifs[index]->initMultiMotifs(s);
					} else if (length!=ret) {//different length, wrong
						return -1;
					} 	
					motifset->motifs[index]->m_index = 0;
					for(int k=0; k<motifset->motifs[index]->m_motifsetslen[s];k++) {//copy from buffer to motif
						motifset->motifs[index]->insert(buffer[k],i,s);
					}
				}
			}
			motifset->motifs[index]->motifsetshittimes[s] =calHit(motifset->motifs[index]->m_motifsets[s],motifset->motifs[index]->m_motifsetsnum[s],motifset->motifs[index]->m_motifsetslen[s]);
		}
		index++;
	}
	for(int j=0; j<motifset->numberofmotif; j++)
		motifset->motifs[j]->printmultimotifs();
	return 0;
}

int CInformation::readMotifCase3(FILE *hfile)
{
	int index = 0;
	int ret = 0;
	int length = 0;
	int num = 0;
	char buffer[LINEMAX];
	int results[MAX_LEN];
	int alphbetindex = 0;
	CConvertor* convertor = new CConvertor(motifset->numberofmotif);
	while (index<motifset->numberofmotif) {
		//read until @
		if (readLine(hfile,NULL,0,'@',-1)==-1) {
			return -1;
		}
		//first read name
		if (readName(hfile,index)==-1) {
			return -1;
		}
		//then read the motif
		length = 0;
		for(int i=0; i<ALPHA_NUM; i++) {
			memset(buffer,0,sizeof(char)*LINEMAX);
			ret = readLine(hfile,buffer,LINEMAX,'\n',0);
			while (ret==0) {
				ret = readLine(hfile,buffer,LINEMAX,'\n',0);
			}
			if (ret==-1) {
				delete convertor;
				return -1;
			} 
			memset(results,0,sizeof(int)*MAX_LEN);
			num = parserLine(results,MAX_LEN,buffer,ret,&alphbetindex);
			if (num<1) {
				delete convertor;
				return -1;
			} else if (length==0) {
				length = num;
				//convertor->motiflen[index] = length;
				convertor->initAll(index,length);
			} else if (num!=length) {//different length
				delete convertor;
				return -1;
			}
			convertor->duplicate(index,results,alphbetindex);
		}
		motifset->motifs[index]->type = 1;
		index++;
	}
	convertor->transformation(motifset, region, -1, rc);
	delete convertor;
	/*for(int j=0; j<motifset->numberofmotif; j++) {
		motifset->motifs[j]->printmultimotifs();
		printf("%d\n",motifset->motifs[j]->hittimes);
	}*/
	return 0;
}

/************************************************************************/
/* read the name of motif with index in the motifset                    */
/************************************************************************/
int CInformation::readName(FILE *hfile, int index)
{
	if(hfile==NULL) {
		return -1;
	}
	//the following two char shound be "ID"
	char name[3];
	memset(name,0,sizeof(char)*3);
	if (fread(name,sizeof(char),2,hfile)!=2) {
		return -1;
	}
	if (strcmp(name,"ID")!=0) {
		return -1;
	}
	char line[MAXNAME];
	memset(line,0,sizeof(char)*MAXNAME);
	if (readLine(hfile,line,MAXNAME,'\n',3)>0) {
		memcpy(motifset->motifs[index]->m_motifname,line, MAXNAME*sizeof(char));
	} else {
		//read nothing
		char buffer[32];
		//buffer[0] = 'M';
		//buffer[1] = 'O';
		//buffer[2] = 'T';
		//buffer[3] = 'I';
		//buffer[4] = 'F';
		//buffer[5] = '_';
		//_itoa(index,buffer+sizeof(char)*6,10);
		sprintf(buffer,"MOTIF_%d",index);
		memcpy(motifset->motifs[index]->m_motifname,buffer,sizeof(char)*32);
	}
	return 0;
}

int CInformation::readSize(FILE *hfile, int index, int setindex)
{
	if(hfile==NULL) {
		return -1;
	}
	//read until @
	if (readLine(hfile,NULL,0,'@',-1)==-1) {
		return -1;
	}
	//the following four chars should be "SIZE"
	char name[5];
	memset(name,0,sizeof(char)*5);
	if (fread(name,sizeof(char),4,hfile)!=4) {
		return -1;
	}
	if (strcmp(name,"SIZE")!=0) {
		return -1;
	}
	char line[MAXNAME];
	memset(line,0,sizeof(char)*MAXNAME);
	int size = 0;
	if (readLine(hfile,line,MAXNAME,'\n',4)>0) {
		size = atoi(line);
		if (size>0) {
			if (motifset->type==1) {
				motifset->motifs[index]->m_multinum = size;
			} else if (motifset->type==2) {
				motifset->motifs[index]->m_motifsetsnum[setindex] = size;
			} else {
				return -1;
			}
			return 0;
		} else {
			return -1;
		}
	} else {
		return -1;
	}
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
int CInformation::parserLine(int *results, int maxresult, char *buffer, int maxlen, int* alphbetindex)
{
	if (result==NULL || buffer==NULL) {
		return -1;
	}
	//find the first char in A,C,G,T
	int pivot = 0;
	while (pivot<maxlen && alp->checkValid(1,buffer[pivot])==false) {
		pivot++;
	}
	if (pivot==maxlen) {//can't find A,C,G,T
		return -1;
	} 
	*alphbetindex = alp->getIndex(buffer[pivot]);
	int numofdouble = 0;
	int indexofdouble = 0;
	pivot++;
	char value[64];
	char* endptr;
	while (numofdouble<maxresult) {
		indexofdouble =0;
		memset(value,0,sizeof(char)*64);
		while (pivot<maxlen) {
			if ((buffer[pivot]>'0'-1 && buffer[pivot]<'9'+1) || buffer[pivot]=='.') {//check the current 
				value[indexofdouble++] = buffer[pivot];
			} else {//the end of a double value
				if (indexofdouble>0) {//read something,then conver it
					results[numofdouble++] = strtol(value,&endptr,10);
				} 
				indexofdouble=0;
				pivot++;
				break;
			}
			pivot++;
		}
		if (pivot==maxlen) {//read the end and convert the last one
			if (indexofdouble>0) {
				results[numofdouble++] = strtol(value, &endptr,10);
			}
			break;
		}
	}
	return numofdouble;
}

int CInformation::readMotifCase4(FILE *hfile)
{
	int index = 0;
	int ret = 0;
	int length = 0;
	int num = 0;
	char buffer[LINEMAX];
	int results[MAX_LEN];
	int alphbetindex = 0;
	CConvertor* convertor[MOTIFSETS];
	int c;
	for(c=0; c<MOTIFSETS; c++) {
		convertor[c]= new CConvertor(motifset->numberofmotif);
	}
	while (index<motifset->numberofmotif) {
		//read until @
		if (readLine(hfile,NULL,0,'@',-1)==-1) {
			return -1;
		}
		//first read name
		if (readName(hfile,index)==-1) {
			return -1;
		}
		//then read the motif
		for(int setindex=0; setindex<MOTIFSETS; setindex++) {
			length = 0;
			for(int i=0; i<ALPHA_NUM; i++) {
				memset(buffer,0,sizeof(char)*LINEMAX);
				ret = readLine(hfile,buffer,LINEMAX,'\n',0);
				while (ret==0) {
					ret = readLine(hfile,buffer,LINEMAX,'\n',0);
				}
				if (ret==-1) {
					for(c=0; c<MOTIFSETS; c++) {
						delete convertor[c];
					}
					return -1;
				} 
				memset(results,0,sizeof(int)*MAX_LEN);
				num = parserLine(results,MAX_LEN,buffer,ret,&alphbetindex);
				if (num<1) {
					for(c=0; c<MOTIFSETS; c++) {
						delete convertor[c];
					}
					return -1;
				} else if (length==0) {
					length = num;
					//convertor->motiflen[index] = length;
					convertor[setindex]->initAll(index,length);
				} else if (num!=length) {//different length
					for(c=0; c<MOTIFSETS; c++) {
						delete convertor[c];
					}
					return -1;
				}
				convertor[setindex]->duplicate(index,results,alphbetindex);
			}
		}
		motifset->motifs[index]->type = 2;
		index++;
	}
	for(c=0; c<MOTIFSETS; c++) {
		convertor[c]->transformation(motifset,region,c,rc);
		delete convertor[c];
	}
	/*for(int j=0; j<motifset->numberofmotif; j++) {
		motifset->motifs[j]->printmultimotifs();
	//	printf("%d\n",motifset->motifs[j]->hittimes);
	}*/
	return 0;	
}

int CInformation::calHit(CAlphabet::ALPHABET* m, int mlen)
{
	if (m==NULL) {
		return 1;
	}
	if (rc==false) {//not consider the reverse and compliment sequences
		return calHitonSeq(region->sequences[0],region->m_regionlen, m, mlen, NULL);			
	} else {
		int hit[2];
		int hitnum = 0;
		CList* hitpos[2];
		hitpos[0] = new CList();
		hitpos[1] = new CList();
		hit[0] = calHitonSeq(region->sequences[0],region->m_regionlen,m,mlen,hitpos[0]);
		hit[1] = calHitonSeq(region->rcsequences[0],region->m_regionlen,m,mlen,hitpos[1]);
		hitnum = hit[0];
		CNode* pivot1 = hitpos[1]->head;
		for(int j=0; j<hit[1]; j++) {//check whether add hittims[1]
			int index = 0;
			CNode* pivot0 = hitpos[0]->head;
			while (pivot0!=NULL && index<hit[0] && (pivot0->position+pivot1->position+mlen != region->m_regionlen)) {
				index++;
				pivot0 = pivot0->next;
			}
			if (index==hit[0]) {//add it
				hitnum++;
			}
			pivot1 = pivot1->next;
		}
		delete hitpos[0];
		delete hitpos[1];
		if (hitnum<1) {
			return 1;
		} else 
			return hitnum;
	}
}


int CInformation::calHitonSeq(char* seq, int seqlen, CAlphabet::ALPHABET* m, int mlen, CList* hitpostion)
{
	int index;
	int hitnum = 0;
	for(int i=0; i<seqlen-mlen; i++) {
		//check wether seq[i]~seq[i+mlen-1] is the same to m
		index = i;
		while (index<i+mlen && alp->compatible[m[index-i]][alp->getIndex(seq[index])]==true) {
			index++;
		}
		if (index==i+mlen) {//find a hit
			hitnum++;
			if (hitpostion!=NULL) {
				hitpostion->insert(new CNode(i));
			}
		}
	}
	if (hitnum<1 && hitpostion==NULL) {
		return 1;
	} else 
		return hitnum;
}

int CInformation::calHit(CAlphabet::ALPHABET** m, int mnum, int mlen)
{
	if (m==NULL) {
		return 1;
	}
	if (rc==false) {//not consider the reverse and compliment sequences
		return calHitonSeq(region->sequences[0],region->m_regionlen, m, mnum, mlen, NULL);			
	} else {
		int hit[2];
		int hitnum = 0;
		CList* hitpos[2];
		hitpos[0] = new CList();
		hitpos[1] = new CList();
		hit[0] = calHitonSeq(region->sequences[0],region->m_regionlen,m, mnum, mlen,hitpos[0]);
		hit[1] = calHitonSeq(region->rcsequences[0],region->m_regionlen,m, mnum, mlen,hitpos[1]);
		hitnum = hit[0];
		CNode* pivot1 = hitpos[1]->head;
		for(int j=0; j<hit[1]; j++) {//check whether add hittims[1]
			int index = 0;
			CNode* pivot0 = hitpos[0]->head;
			while (pivot0!=NULL && index<hit[0] && (pivot0->position+pivot1->position+mlen != region->m_regionlen)) {
				index++;
				pivot0 = pivot0->next;
			}
			if (index==hit[0]) {//add it
				hitnum++;
			}
			pivot1 = pivot1->next;
		}
		delete hitpos[0];
		delete hitpos[1];
		if (hitnum<1) {
			return 1;
		} else 
			return hitnum;
	}
}

int CInformation::calHitonSeq(char* seq, int seqlen, CAlphabet::ALPHABET** m, int mnum, int mlen,CList* hitpostion)
{
	int index;
	int hitnum = 0;
	for(int i=0; i<seqlen-mlen; i++) {
		//check wether seq[i]~seq[i+mlen-1] is the same to m
		for(int j=0; j<mnum; j++) {
			index = i;
			while (index<i+mlen && alp->compatible[m[j][index-i]][alp->getIndex(seq[index])]==true) {
				index++;
			}
			if (index==i+mlen) {//find a hit
				hitnum++;
				if (hitpostion!=NULL) {
					hitpostion->insert(new CNode(i));
				}
				break;
			}
		}
	}
	if (hitnum<1 && hitpostion==NULL) {
		return 1;
	} else 
		return hitnum;
}
