#ifndef EM_H
#define EM_H

#include "../refinementPhase/MotifContainer.h"

#include "hoUtils.h"

class EM{
public:

	EM( MotifContainer& startModels ) : _models( startModels ){}

	MotifContainer& getModels(){
		return _models;
	}

	void cv( int folds );

	void go( int N, int sequenceNr[] );
	void go();

	void score();

private:

	MotifContainer& _models;
};

#endif /* EM_H */
