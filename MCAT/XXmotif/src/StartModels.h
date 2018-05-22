#ifndef OVERREPPATTERNSCOMPOSITE_H_
#define OVERREPPATTERNSCOMPOSITE_H_

#include "refinementPhase/MotifContainer.h"
#include "elongationPhase/elongationCandidates.h"
#include "aminoacids/madonaPro.h"

class StartModels {

public:
	void findInitialMotifs(MotifContainer&, double pValThreshold);
	void initStartMotif(MotifContainer&);

private:
	void addElongatedToStartModels(elongList&, MotifContainer&);

	int removeRedundantMotifs(elongList&);

	void printOverrepPatternsIUPAC(elongList&, int = 10);
};

#endif /*OVERREPPATTERNSCOMPOSITE_H_*/
