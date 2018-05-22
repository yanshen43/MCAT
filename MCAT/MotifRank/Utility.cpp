// Utility.cpp: implementation of the CUtility class.
//
//////////////////////////////////////////////////////////////////////

#include "Utility.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CUtility::CUtility()
{

}

CUtility::~CUtility()
{

}

/************************************************************************/
/*  Init the information of sequences, motif and matrix                 */
/*  @param info: the point to CInformation, shouldn't be a NULL point
 *	@param argc: the size of argv
 *  @param argv: the filename of sequences, motif and matrix
 *  @return:
 *  0	if the information is correct
 *  -1	otherwise
 */
/************************************************************************/
int CUtility::readFile(CInformation* info,int argc, char* argv[]) 
{
	if (info==NULL) {
		return -1;
	}
	int ret;
	/*init the rc*/
	info->rc = (atoi(argv[4])==1);
	/*init sequence file*/
	FILE* seqfile = fopen(argv[1],"r");
	if (seqfile==NULL) {
		printf("can't open sequences file.\n");
		return -1;
	} else {
		ret = info->initSequences(seqfile);
		fclose(seqfile);
		if (ret==-1) {
			printf("the file format of sequences file is wrong.\n");
			return -1;
		} else if (ret==-2)	{
			printf("can't support multiply sequences now.\n");
			return -1;
		}
	}
	/*init matrix file*/
	FILE* matrixfile = fopen(argv[3],"r");
	if (matrixfile==NULL) {
		printf("can't oepn matrix file.\n");
		return -1;
	} else  {
		ret = info->initMatrix(matrixfile);
		fclose(matrixfile);
		if (ret==-1) {
			printf("the file format of matrix file is wrong.\n");
			return -1;
		}
	}
	/*init motif file*/
	FILE* motiffile = fopen(argv[2],"r");
	if (motiffile==NULL) {
		printf("can't open motif file.\n");
		return -1;
	} else {
		ret = info->initMotif(motiffile);
        strncpy(info->m_filename, argv[2], 64);
		fclose(seqfile);
		if (ret==-1) {
			printf("the file format of motif file is wrong.\n");
			return -1;
		}
	}
	return 0;
}
