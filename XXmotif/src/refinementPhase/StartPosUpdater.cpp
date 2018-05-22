#include "StartPosUpdater-inl.h"
#include "empiricalPvalCalibration.h"


StartPosUpdater::StartPosUpdater() : pValCalculator(PVal_Calculator::getInstance()) {

	/* correction factor: minimum motif length is 5 for na and 2 for aa */
	sitesAlloc1 		= Global::posSet->total[Global::posSet->nent] / (Global::aa ? 2 : 5);
	sortedSites  		= (sorted_sites*)malloc( sitesAlloc1 * sizeof(sorted_sites));
	sortedSites_elong  	= (sorted_sites*)malloc( sitesAlloc1 * sizeof(sorted_sites));
	sitesAlloc2 		= Global::posSet->total[Global::posSet->nent];
	bestStartPos 		= (sorted_sites*)malloc( sitesAlloc2 * sizeof(sorted_sites));

	motifsPerSequenceCount = (int*)malloc((Global::posSet->nent+1) * sizeof(int));

	for(int i=0; i<sitesAlloc1; i++){
		sortedSites[i]  	  = (sorted_sites)malloc(sizeof(sorted_sites_type));
		sortedSites_elong[i]  = (sorted_sites)malloc(sizeof(sorted_sites_type));
	}
	for(int i=0; i<sitesAlloc2; i++){
		bestStartPos[i] = (sorted_sites)malloc(sizeof(sorted_sites_type));
	}

	pvalCons = (double*)malloc(sitesAlloc1* sizeof(double));

	kmer = (unsigned char*) calloc(PWM_LENGTH, sizeof(unsigned char));

	log_posSetSites = (double*)malloc((Global::posSet->nent+1) * sizeof(double));
	for(int i=1; i<=Global::posSet->nent; i++){
		log_posSetSites[i] = LogTable::LOG_i[Global::posSet->total[i]];
	}
}

StartPosUpdater::~StartPosUpdater(){
	for(int i=0; i<sitesAlloc1; i++){
		free(sortedSites[i]);
		free(sortedSites_elong[i]);
	}
	for(int i=0; i<sitesAlloc2; i++){
		free(bestStartPos[i]);
	}
	free(sortedSites);
	free(sortedSites_elong);
	free(bestStartPos);

	free(pvalCons);

	free(kmer);
	free(log_posSetSites);
}

void StartPosUpdater::fill_startPosList_with_bestSubset(Motif* motif, double& bestPval, int& bestMotifNb){
	if(nbSites == 0) return;
	sortSitesByScore(sortedSites, nbSites);
	if(sortedSites[0]->score == 1) return;
	int bestCounter = 1; int bestMotifNbPlus = 1; int bestCounterPlus = 1;
	getBestSubset(motif->getPWM(), motif->getMotifColumns(), sortedSites, nbSites, motif->getPosSetSize(), bestPval,
			bestMotifNb, bestCounter, bestMotifNbPlus, bestCounterPlus, motif->getTotalSites());

	sortSitesByStartPos(sortedSites, bestCounter);

	fill_startPosList_with_sites(motif->getStartPosList(), motif->getFirstMotifColumn(), bestMotifNb);
}

void StartPosUpdater::set_BindingSites(Motif* motif, const ThresholdChecker &pValThreshold){
	if(nbSites == 0) return;
	sortSitesByScore(sortedSites, nbSites);
	if(sortedSites[0]->score == 1) return;

	int bestMotifNb = 0;
	while(bestMotifNb < nbSites && pValThreshold.satisfies(sortedSites[bestMotifNb]->score)) bestMotifNb++;

	sortSitesByStartPos(sortedSites, bestMotifNb);

	fill_startPosList_with_matches(motif->getBindingSites(), motif->getFirstMotifColumn(), bestMotifNb);
}

//void StartPosUpdater::update_pwm_with_bestSubset_PlusAll(Motif* motif, double pseudo){
//	if(nbBestSites == 0) return;
//	sortSitesByScore(bestStartPos, nbBestSites);
//	if(bestStartPos[0]->score == 1) return;
////	for(int i=0; i<nbBestSites; i++){
////		for(int j=bestStartPos[i]->startPos; j<bestStartPos[i]->startPos+2; j++){
////			fprintf(stderr, "%c", AlphaChar(Global::posSet->entity[bestStartPos[i]->sequence]->S[0][j], Global::A));
////		}
////		fprintf(stderr, "\t%d/%d\t%f\n", bestStartPos[i]->sequence, bestStartPos[i]->startPos, bestStartPos[i]->score);
////	}
//
//	int bestMotifNb_plus = getBestSubset_plus(motif->getPWM(), motif->getMotifColumns(), bestStartPos, nbBestSites, motif->getPosSetSize(), motif->getTotalSites());
//
//	assert(bestMotifNb_plus > 0);
//	//cerr << "bestMotifNb_plus: " << bestMotifNb_plus << endl;
//	/* update PWM */
//	motif->updatePWM_OOL(bestStartPos, motif->getPWM(), motif->getMotifColumns(), bestMotifNb_plus, pseudo);
//}
//
//void StartPosUpdater::update_pwm_with_bestSubset_All(Motif* motif, double pseudo){
//	if(nbSites == 0) return;
//	sortSitesByScore(sortedSites, nbSites);
//	if(sortedSites[0]->score == 1) return;
//
//	int bestMotifNb = 0;
//	double bestPval;
//	int bestCounter; int bestMotifNbPlus; int bestCounterPlus;
//	getBestSubset(motif->getPWM(), motif->getMotifColumns(), sortedSites, nbSites, motif->getPosSetSize(),
//			bestPval, bestMotifNb, bestCounter, bestMotifNbPlus, bestCounterPlus, motif->getTotalSites());
//
//	assert(bestMotifNb > 0);
//	/* update PWM */
//	motif->updatePWM_OOL(sortedSites, motif->getPWM(), motif->getMotifColumns(), bestMotifNb, pseudo);
//}
//
//void StartPosUpdater::update_pwm_with_bestSubset_Plus(Motif* motif, double pseudo){
//	if(nbSites == 0) return;
//	sortSitesByScore(sortedSites, nbSites);
//	if(sortedSites[0]->score == 1) return;
//
//	int bestMotifNb_plus = getBestSubset_plus(motif->getPWM(), motif->getMotifColumns(), sortedSites, nbSites, motif->getPosSetSize(), motif->getTotalSites());
//
//	assert(bestMotifNb_plus > 0);
//	/* update PWM */
//	motif->updatePWM_OOL(sortedSites, motif->getPWM(), motif->getMotifColumns(), bestMotifNb_plus, pseudo);
//	//cerr << *motif;
//}
//
//
//void StartPosUpdater::update_elongPWM_with_bestSubset_Plus(Motif* motif, double **elongPWM, motif_columns_type& tmpColumns,	int& bestMotifNb){
//	if(nbSites == 0) return;
//	sortSitesByScore(sortedSites_elong, nbSites);
//	if(sortedSites_elong[0]->score == 1)return;
//
//	updateElongPWM(motif->getPWM(), elongPWM, tmpColumns);
//
//	int oldSites = motif->getTotalSites(); if(oldSites == 0) oldSites = maxSites;
//
//	//cerr << "getBestSubset_plus" << endl;
//	bestMotifNb = getBestSubset_plus(elongPWM, tmpColumns, sortedSites_elong, nbSites, motif->getPosSetSize(), motif->getTotalSites());
//
//	//cerr << "nbSites: " << nbSites << "\tbestMotifNb: " << bestMotifNb << endl;
//	if(bestMotifNb != 0) motif->updatePWM_OOL(sortedSites_elong, elongPWM, tmpColumns, bestMotifNb, Global::pseudo);
//}

void StartPosUpdater::update_PWM_with_bestSubset(Motif* motif, int& bestMotifNb, double& bestPval, bool plus, double pseudo){
	if(nbSites == 0) return;
	sortSitesByScore(sortedSites, nbSites);
	if(sortedSites[0]->score == 1)return;

	int bestCounter = 1; int bestMotifNbPlus = 1; int bestCounterPlus = 1;
	getBestSubset(motif->getPWM(), motif->getMotifColumns(), sortedSites, nbSites, motif->getPosSetSize(),
			bestPval, bestMotifNb, bestCounter, bestMotifNbPlus, bestCounterPlus, motif->getTotalSites());
	//fprintf(stderr, "bestMotifNb: %d, bestPval: %e\n", bestMotifNb, exp(bestPval));

	assert(bestMotifNb > 0);

	int motifNb = plus ? bestMotifNbPlus : bestMotifNb;

	motif->updatePWM_OOL(sortedSites, motif->getPWM(), motif->getMotifColumns(), motifNb, pseudo);
}


void StartPosUpdater::update_elongPWM_with_bestSubset(Motif* motif, double **elongPWM, motif_columns_type& tmpColumns,
		int& bestMotifNb, double& bestPval, bool plus, double pseudo){
	if(nbSites == 0) return;
	sortSitesByScore(sortedSites_elong, nbSites);
	if(sortedSites_elong[0]->score == 1)return;

	int bestCounter = 1; int bestMotifNbPlus = 1; int bestCounterPlus = 1;
	getBestSubset(elongPWM, tmpColumns, sortedSites_elong, nbSites, motif->getPosSetSize(), bestPval,
			bestMotifNb, bestCounter, bestMotifNbPlus, bestCounterPlus, motif->getTotalSites());
	//fprintf(stderr, "bestMotifNb: %d, bestPval: %e\n", bestMotifNb, exp(bestPval));

	assert(bestMotifNb > 0);

	int motifNb = plus ? bestMotifNbPlus : bestMotifNb;

	motif->updatePWM_OOL(sortedSites_elong, elongPWM, tmpColumns, motifNb, pseudo);
}

double StartPosUpdater::fill_startPosList_with_elongSites(StartPosContainer& realStartPosList, int bestMotifNb, motif_columns_type& tmpColumns){
	if (bestMotifNb == 0) return 1;

	/* store significant oops start positions ( for new matrix ) */
	realStartPosList.clear();
	int firstMotifColumn = tmpColumns.front();

	/* reset motifs per sequence counter */
	memset(motifsPerSequenceCount, 0, (Global::posSet->nent+1)* sizeof(int));

	for(int i=0, motifCounter = 0;;i++){
		if(++motifCounter > bestMotifNb) break;

		motifsPerSequenceCount[sortedSites[i]->sequence]++;
		if(motifsPerSequenceCount[sortedSites[i]->sequence] > Global::maxMotifsPerSequence) continue;

		realStartPosList.push_back(StartPos(sortedSites_elong[i]->sequence,
			static_cast<int32_t>(sortedSites_elong[i]->startPos-firstMotifColumn + 1)));
	}
	return sortedSites_elong[bestMotifNb-1]->score;
}

void StartPosUpdater::update_sites_with_empirical_Pvals(Motif* motif){
	double finalPval, pValScore;
	sortSitesByScore(bestStartPos, nbBestSites);

	int region = motif->getEnrichment().endRegion - motif->getEnrichment().startRegion + 1;
	int length = motif->getMotifLength();
	int size = sizeof(unsigned char)*length;

	memcpy(kmer, Global::posSet->entity[bestStartPos[nbBestSites-1]->sequence]->S[0] + bestStartPos[nbBestSites-1]->startPos, size);

	//fprintf(stderr, "worstKmer: ");
	//for(int i=0; i<motif->getMotifLength(); i++)fprintf(stderr, "%c", AlphaChar(kmer[i], Global::A));
	//fprintf(stderr, " \n");

	/* update pVals of bestStartPos list */
	empiricalPvalCalibration &kmerHash = empiricalPvalCalibration::getInstance();
	kmerHash.reinitialize(motif, kmer);

	sortSitesByStartPos(bestStartPos, nbBestSites);

	//motif_columns_type sortedColumns = pValCalculator.getSortedColumns(motif->getPWM(), motif->getMotifColumns());

	int sites_counter = 0;
	for(int i=0; i<nbBestSites; i++){
		/* initialize sites */
		const int32_t seq = bestStartPos[i]->sequence;
		const int32_t pos = bestStartPos[i]->startPos;

		if(pos < 1 || pos + length - 1 > Global::posSet->entity[seq]->n) continue;

		int len = Global::posSet->entity[seq]->n - length + 1;
//		if(motif->getType() == PALINDROME && Global::revcomp) len = Global::posSet->entity[seq]->n/2 - length + 1;
		int delta_0 = std::min(len, region);

		finalPval = pValCalculator.calculatePval(kmerHash, kmer, motif->getMotifColumns(), \
				seq, pos, len, delta_0, size, length, Global::posSet->avgLength[motif->getPosSetSize()], \
				pValScore, Global::multipleOccurrence);

		//if(seq == 1) fprintf(stderr, "before: %f, after: %f\n ", bestStartPos[i]->score, finalPval);
		//if(seq < 10)fprintf(stderr, "before: %f, after: %f\n ", bestStartPos[i]->score, finalPval);

		//finalPval = std::max(bestStartPos[i]->score, finalPval);
		//fprintf(stderr, "emp: %f\n", finalPval);
		if(finalPval == -1) finalPval = bestStartPos[i]->score;
		if(finalPval == 1) continue;

		//fprintf(stderr, "final Pval: %f\n", finalPval);

		if(Global::multipleOccurrence){
			if(!updateMultOcc(sortedSites, sites_counter, seq, pos, length, finalPval, Global::posSet)){
				continue;
			}
		}else if(finalPval < sortedSites[seq-1]->score){
			//cerr << i << "\t" << j << "\t" << pValScore << "\t" << Global::posSet->avgLength-length+1 << "\t" << finalPval << endl;
			sortedSites[seq-1]->sequence = seq;
			sortedSites[seq-1]->startPos = pos;
			sortedSites[seq-1]->score = static_cast<float>(finalPval);
		}
		sites_counter++;
	}
	if(Global::multipleOccurrence){
		sortedSites[sites_counter]->score = 1;

		nbSites = sites_counter;
		maxSites = Global::posSet->total[motif->getPosSetSize()] - motif->getPosSetSize()*(length-1);
		if(Global::revcomp){
			maxSites = Global::posSet->total[motif->getPosSetSize()] - 2*motif->getPosSetSize()*(length-1);
		}
	}else{
		nbSites = motif->getPosSetSize();
		maxSites = motif->getPosSetSize();
	}
}

StartPosUpdater& StartPosUpdater::getInstance()
{
  static StartPosUpdater instance;
  return instance;
}

