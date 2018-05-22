#ifndef _DISOCONS_H
#define _DISOCONS_H

#include "MProGlobal.h"
#include "AbstractSupplementaryInformationProvider.h"
//#include "../refinementPhase/Motif.h"

class AbstractKmer;

#include <vector>

class DisoCons : public AbstractSupplementaryInformationProvider {

private:
	typedef std::vector<std::vector<double> > vecvec_t;
	static const int dcMatrixRes = 100;
	/** precomputed values for disorder and conservation */
	vecvec_t diso;
	vecvec_t cons;
	/** P values are stored in Pdc[diso][cons] */
	vecvec_t Pdc;
	void readDCFile(const std::string &filename,
		const std::vector<Sequence> &sset, vecvec_t &diso, vecvec_t &cons);

	static DisoCons* instance;

	/**
	 * Initialize data structures.
	 */
	DisoCons(const MProGlobal::SeqInfo_t &S);
	~DisoCons();

public:

	static AbstractSupplementaryInformationProvider const* getInstance();

	/**
	 * Compute disorder/conservation P value for the positions defined
	 * by the true match positions of the given kmer in sequence seq,
	 * starting at position startPos.
	 */
	double getPvalue(AbstractKmer const* const kmer, const int &seq,
			const int &startPos) const;

	double getPvalue(motif_columns_type const * const cols, const int &seq,
			const int &startPos) const;

};

#endif
