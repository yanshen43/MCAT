#include "backgroundDistribution.h"
#include "Globals.h"
#include "LogTable.h"
#include "StartModels.h"
#include "utils.h"
#include "em/EM.h"
#include "em/hoNullModel.h"
#include "em/hoUtils.h"
#include "nucleotides/conservationScores.h"
#include "refinementPhase/MotifContainer.h"
#include "refinementPhase/Refine_Motifs.h"

int main(int argc, char *argv[])
{
	/* wrapper method with runtime measurement */
	fprintf(stderr, "\n");
	fprintf(stderr, "====================\n");
	fprintf(stderr, "==  run XXmotif   ==\n");
	fprintf(stderr, "====================\n");

    long time1 = time( NULL );

	Global G( argc, argv ); /* Initialize global variables */

	LogTable L; /* Precalculate logs */

	printf("posSet - entries: %d, maxMultSeq: %d, minLength: %d, maxLength: %d, total: %d\n",
			Global::posSet->nent, Global::posSet->max_MultSeq, Global::posSet->min_leng, Global::posSet->max_leng, Global::posSet->total[Global::posSet->nent]);
	if(Global::negSet != NULL)printf("NegSet - entries: %d, maxMultSeq: %d, minLength: %d, maxLength: %d, total: %d\n",
			Global::negSet->nent, Global::negSet->max_MultSeq, Global::negSet->min_leng, Global::negSet->max_leng, Global::negSet->total[Global::negSet->nent]);

	if( !Global::aa ){
		setConservationProbs();
//		L.setMotifNbCorrection();
	}

	MotifContainer startModels;
	StartModels seeds;

	Motif::initHoMotif();

	if( Global::aa )
		MadonaPro::initialize( 6, 3, 7 );

	if( Global::bindingSiteFile == NULL && Global::profFile == NULL && Global::startMotif == NULL ){
		seeds.findInitialMotifs( startModels, 1e-3 );
	}
	else{
		seeds.initStartMotif( startModels );
	}

	if( Global::em ){
		EM em( startModels );
		em.go();
		if( Global::cv )
			if( Global::testSet )
				em.score();
	}
	else if( !( Global::cv ) ){

		MotifRefinement refinement( startModels );
		refinement.start();
	}
	
	printf("\n");
	fprintf(stderr,"\n------------time: %ld seconds (%0.2f minutes)---------------\n",
        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);

	return 0;	
}
