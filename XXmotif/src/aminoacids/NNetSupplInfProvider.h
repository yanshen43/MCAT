#ifndef _NNET_SUPPL_INF_PROVIDER_H
#define _NNET_SUPPL_INF_PROVIDER_H

#include <vector>

#include "MProGlobal.h"
#include "AbstractSupplementaryInformationProvider.h"

class NNet;

class NNetSupplInfProvider: public AbstractSupplementaryInformationProvider {

  public:
    /* implemented as a singleton */
    static AbstractSupplementaryInformationProvider const* getInstance();

    double getPvalue(AbstractKmer const* const kmer, const int &seq,
        const int &startPos) const;

    double getPvalue(motif_columns_type const * const cols, const int &seq,
        const int &startPos) const;

  private:
    static NNetSupplInfProvider* instance;
    static NNet* nn;

    /* Externally calculated P-values are sampled in 1000 bins,
     * i.e. with a resolution of 1e-3.
     * TODO read RES from P-value file
     */
    static const int RES = 1000;
    std::vector<double> ptable;
    /* Precalculated table of log(P) for efficency. */
    std::vector<double> log_ptable;

    /* sequence properties for NN prediction,
     * currently five properties per residue are used
     * TODO parameter instead of constant
     */
    static const int NUM_PROPS = 5;
    std::vector<std::vector<std::vector<double> > > props;

    /**
     * Constructor and destructor
     */
    NNetSupplInfProvider(const MProGlobal::SeqInfo_t &S);
    ~NNetSupplInfProvider();

};

#endif
