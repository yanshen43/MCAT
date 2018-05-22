#include "Information.h"
#include "Utility.h"
#include "Rank.h"
#define ARGCNUM 5

int main(int argc, char *argv[])
{
	if (argc!=ARGCNUM) {
		printf("Invocation:\n%s sequence_file motif_file matrix_file <0, 1 double standed?>\n",
                argv[0]);
		exit(-1);
	} else {
		CRank* rank = new CRank();
		CUtility* cu = new CUtility();
		CInformation* info =new CInformation();
		int ret = cu->readFile(info, argc, argv);
		if (ret==0) {
			rank->doRanking(info);
		} 
		delete info;
		delete rank;
		delete cu;
	}
	return 0;
}
