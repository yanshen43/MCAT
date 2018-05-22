#include "Refine_Motifs.h"
#include "../output.h"
#include "../aminoacids/madonaPro.h"

/*find best Motifs by iterating till convergence */
void MotifRefinement::start(){
	int oldMergedMotifsNb = _models.getMotifNb();

	printf("Number of motifs: %d\n\n", _models.getMotifNb());

	if(!Global::noRefinementPhase){
		printf("\n\tITERATIVELY REFINE PWMS\n\n");

		bool lastIteration = false;
		int round = 1;
		do{
			if(round != 1 && _models.getMotifNb() == oldMergedMotifsNb) lastIteration = true;
			//if(_round != 1 && _models.getMotifNb() == oldMergedMotifsNb){
			//	if(!_maximizeMotifLength) _maximizeMotifLength = true;
			//	else _lastIteration = true;
			//}
			/* find best motif by iterating them */
			if (Global::maxIterations > 0){
				_models.iterate(Global::maximizeMotifLength, lastIteration);
			}
			if(_models.getMotifNb() == 0 || Global::maxIterations<=1) lastIteration = true;

			/* sort found motifs*/
			_models.sort_and_filter(100);

			if(!lastIteration){
				oldMergedMotifsNb = _models.getMotifNb();
				/* merge motifs and remove overlapping motifs */
				//cerr << endl << "merge models" << endl;
				_models.merge(true);
				//cerr << endl << "merged" << endl;
				//_models.merge(false);
			}
			round++;
		}while(!lastIteration);
	}

	/* Rescale E-values by taking square root.
	 * (Improved correspondance between theoretical and reported
	 * E-values on random sequences)
	 */
	const list<Motif*> &motifs = _models.getMotifs();
	for (list<Motif*>::const_iterator mit=motifs.begin(); mit!=motifs.end(); ++mit) {
	  (*mit)->setPval(0.5 * (*mit)->getPval());
	}

	int maxFilter = std::numeric_limits<int>::max();
	if(Global::noRefinementPhase) maxFilter = 10;
	_models.filter(5, maxFilter, Global::finalFilterThreshold);

//	_models.updatePWM(1e-10);

	_models.setBindingSites(Global::instanceThreshold);

	if (Global::removeHomology) {
		_models.removeHomologousMotifs();
	}

//  /* filter motifs differing too much from their reverse complement */
//  if (!Global::aa) {
//    std::cout << std::endl;
//    _models.revcompFilter();
//  }

	/* write results into files */
	std::cout << std::endl << "********** Final results **********" << std::endl << std::endl;

	Output::printOutput(_models, 10);
}
