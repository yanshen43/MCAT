// MotifInfo.h: interface for the CMotifInfo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MOTIFINFO_H__4DB9849E_3EEC_4820_B7B9_5AE004281B9B__INCLUDED_)
#define AFX_MOTIFINFO_H__4DB9849E_3EEC_4820_B7B9_5AE004281B9B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Alphabet.h"
#include "Prefix.h"
#include "stdlib.h"
#include "stdio.h"
#define MAX_LEN 64
#define MOTIFSETS 2
#define MAXNAME 64

class CMotifInfo  
{
public:	
	void preProcess();
	bool insert(char c, int motifindex, int setindex);
	void initMotif();	
	void initMultiMotifs(int index);
	void printmultimotifs();
	CMotifInfo();
	virtual ~CMotifInfo();
	int m_alphabet;		// the size of alphabet for motif
	int m_prefixtotalsize;    //total number of prefix set
	CPrefix** prefixset;	  //the list of all prefix,from the shortest to the longest
	int* skiplist[ALPHA_NUM]; //the skip list, skiplist[i][j] is the index of Sm(prefixset[j]+i) in prefixset
	int hittimes;			  //times the motif hits the region
	int type;				 
	
	//the type of motifs, 0 means deterministic motif,1 for multiply motifs, 2 for two motifs
	char m_motifname[MAXNAME];		//the name of the motif
	////////////////TYPE 0: member for deterministic motif on basic alphabet or IUPAC//////
	CAlphabet::ALPHABET* m_motif;	//the information of the motif
	int m_motiflen;					// the input length of motif
	////////////////TYPE 1: member for multiply motifs//////////////////////////////////
	CAlphabet::ALPHABET** m_multimotifs;//size of m_multimotifs is m_multinum*m_multilen
	int m_multinum;		//number of motif words in the set	
	bool fdo;
	////////////////TYPE 2: member for two motif sets////////////////////////////////////	
	CAlphabet::ALPHABET** m_motifsets[MOTIFSETS];
	int m_motifsetslen[MOTIFSETS];	//motif length of each motif set
	int m_motifsetsnum[MOTIFSETS];   //number of motifs in each set
	int motifsetshittimes[MOTIFSETS];		//hittimes of motifs in each set
	int* multiskiplist[MOTIFSETS+1][ALPHA_NUM];
	CPrefix** multiprefixset[MOTIFSETS+1];
	int m_multiprefixtotalsize[MOTIFSETS+1];
	int* m_fullhit[MOTIFSETS+1];
	int** m_cutnum[MOTIFSETS+1][ALPHA_NUM][MOTIFSETS];
	int* degskiplist[MOTIFSETS];
	int** m_degcutnum[MOTIFSETS];
	///////////////////////////////////////////////////////////////////////////////
	int m_index;		// the length
	CAlphabet* p_alpha;
	int* m_prefixsize;
	int* m_prefixsizestart;
	// in skiplist[i][j], j is the index of current prefix, i is the next char alphabet
	// skiplist[i][j] is the index of prefix after skip
private:
	void degcunNumber();
	void printfdegskiplist();
	void cutNumber(int setindex);
	int initCutNumber(int setindex, int addchar, int multimoitfindex, int prefixindex,int boundary);
	int Combine(int* prefixtotalsize, CPrefix** prefix[MOTIFSETS], int* prefixsizestart[MOTIFSETS], CPrefix** comprefix);
	void preProcessMotifSets();
	void preProcessMulti();
	void initMultiSkiplist(int* skip[ALPHA_NUM], CPrefix** prefix, int prefixsize);
	int initMultiPrefixset(int totallength, CPrefix** prefix, int* startpos, CAlphabet::ALPHABET** sequences, int* hits, int setindex);
	void printskiplist();
	void findstartend(int* start, int* end);
	int findprefixoflen(CPrefix* newprefix, int len);
	void printprefix();
	void initSkiplist();
	void initPrefixset();
	void initPrefixsize();
	int findprefix(int index, CPrefix* oldprefix, CAlphabet::ALPHABET add);
	int findmultiprefix(CPrefix** prefixes, CPrefix* oldprefix, CAlphabet::ALPHABET add, int prefixsize, bool cut);
		
	void releasePrefix();
};

#endif // !defined(AFX_MOTIFINFO_H__4DB9849E_3EEC_4820_B7B9_5AE004281B9B__INCLUDED_)
