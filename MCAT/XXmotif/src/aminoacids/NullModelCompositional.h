#ifndef NULLMODELCOMPOSITIONAL_H_
#define NULLMODELCOMPOSITIONAL_H_

#include "../UngappedKmer.h"

class NullModelCompositional {

  public:
  static double getProbability(unsigned char const * const kmer, const int length);
  static double getProbability(char const * const kmer, const int length = -1);
  static double getProbability(const UngappedKmer &kmer);

  static void initialize();

  static double getConditional(const uint64_t &index) {
    return bgCompositionalCond[index];
  }

//  static float const* getProbabilities() {
//    return bgCompositionalProb;
//  }
//
//  static float const* getConditionals() {
//    return bgCompositionalCond;
//  }

  static void destruct() {
    if (bgCompositionalProb != NULL) free(bgCompositionalProb);
    if (bgCompositionalCond != NULL) free(bgCompositionalCond);
  }

private:

  typedef float REAL;
  static void countCompositionProbs();
  static void compositionCount_rec(const int startPos, const int maxPos, const int order, unsigned char *kmer, double *counts, const unsigned char * const seq);

  static int      modelOrder;   /** order of background model */
  static REAL*   bgCompositionalProb;  /** negset composition-dependent conditional probabilites */
  static REAL*   bgCompositionalCond;  /** negset composition-dependent conditional probabilites */
  static REAL    alpha;      /** average number of pseudocounts per k-mer */
};

#endif /* NULLMODELCOMPOSITIONAL_H_ */
