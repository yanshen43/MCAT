#include <cstdlib>

#include "NullModelCompositional.h"
#include "MProGlobal.h"
#include "../NullModel.h"
#include "../Globals.h"
#include "../UngappedKmer.h"

int NullModelCompositional::modelOrder;
NullModelCompositional::REAL* NullModelCompositional::bgCompositionalProb;
NullModelCompositional::REAL* NullModelCompositional::bgCompositionalCond;
NullModelCompositional::REAL NullModelCompositional::alpha;

void NullModelCompositional::initialize() {
  modelOrder = NullModel::getOrder();
  alpha = NullModel::getAlpha();
  countCompositionProbs();
}

void NullModelCompositional::compositionCount_rec(const int startPos,
    const int maxPos, const int order, unsigned char *kmer, double *counts,
    const unsigned char * const seq) {
  for (int currentPos = startPos; currentPos <= maxPos; ++currentPos) {
    kmer[order] = seq[currentPos];
    ++bgCompositionalProb[UngappedKmer(kmer, order + 1)];
    ++counts[order];
    if (order < modelOrder && currentPos < maxPos) {
      compositionCount_rec(currentPos + 1, maxPos, order + 1, kmer, counts, seq);
    }
  }
}

void NullModelCompositional::countCompositionProbs() {
  const int window_size = 11; /* size of sliding window */
  const int max_k = modelOrder + 1;
  unsigned char* const kmer_chars = new unsigned char[max_k];
  double* const counts = new double[max_k];

  ss_type set = Global::negSet != NULL ? Global::negSet : Global::posSet;

  UngappedKmer maxKmer(UngappedKmer::greatestWithLength(modelOrder + 1));
  bgCompositionalProb = static_cast<REAL*>(calloc(maxKmer + 1, sizeof(REAL)));
  bgCompositionalCond = static_cast<REAL*>(calloc(maxKmer + 1, sizeof(REAL)));

  /* initialize counts */
  for (int order = 0; order <= modelOrder; ++order) {
    counts[order] = 0;
  }

  /* For each position in pos/negSet, ... */
  for (int seq = 1; seq <= set->nent; ++seq) {
    for (int startPos = 1; startPos <= set->entity[seq]->n; ++startPos) {
      int maxPos = std::min(startPos + window_size - 1,
          (int) set->entity[seq]->n);
      for (int pos = startPos; pos <= maxPos; ++pos) {
        if (set->entity[seq]->S[0][pos] == 0) {
          maxPos = pos - 1;
          break;
        }
      }
      kmer_chars[0] = set->entity[seq]->S[0][startPos];
      if (modelOrder > 0) {
        compositionCount_rec(startPos + 1, maxPos, 1, kmer_chars, counts,
            set->entity[seq]->S[0]);
      }
    }
  }

  /* copy probabilities for 0th order */
  double sum_prob = 0;
  double sum_cond = 0;
  for (unsigned char i = 1; i <= nAlpha(Global::A); ++i) {
    bgCompositionalProb[UngappedKmer(&i, 1)]
        = NullModel::getProbs()[NullModel::sub2ind(&i, 0)];
    bgCompositionalCond[UngappedKmer(&i, 1)]
        = NullModel::getConds()[NullModel::sub2ind(&i, 0)];
    sum_prob += bgCompositionalProb[UngappedKmer(&i, 1)];
    sum_cond += bgCompositionalCond[UngappedKmer(&i, 1)];
    assert(sum_prob == sum_cond); /* must be identical for 0th order */
  }
  bgCompositionalProb[UngappedKmer(std::string(1, MProGlobal::getS().alphabet.wildcardChar()))] = static_cast<REAL>(sum_prob);
  bgCompositionalCond[UngappedKmer(std::string(1, MProGlobal::getS().alphabet.wildcardChar()))] = static_cast<REAL>(sum_cond);

  /* now, normalize and compute conditional probs */
  for (int order = 1; order <= modelOrder; ++order) {

    const UngappedKmer lastK = UngappedKmer::greatestWithLength(order + 1);

    /* fudge dollar */
    for (UngappedKmer kmer = UngappedKmer::smallestWithLength(order + 1); kmer
        <= lastK; ++kmer) {
      if (UngappedKmer::containsDollar(kmer)) {
        counts[order] += (Global::dollarScaleFactor - 1)
            * bgCompositionalProb[kmer];
        bgCompositionalProb[kmer] *= Global::dollarScaleFactor;
      }
    }

    /* marginalize true counts */
    for (UngappedKmer kmer = UngappedKmer::smallestWithLength(order + 1); kmer
        <= lastK; ++kmer) {
      if (kmer.lastWildcardPosition() != UngappedKmer::undef) {
        UngappedKmer::ungapped_kmer_list_t expansions(kmer.expansions());
        bgCompositionalProb[kmer] = 0;
        for (UngappedKmer::ungapped_kmer_list_t::const_iterator it =
            expansions.begin(); it != expansions.end(); ++it) {
          bgCompositionalProb[kmer] += bgCompositionalProb[*it];
        }
      }
    }

    /* add pseudocounts from lower order and normalize */
    const double numberOfPseudocounts = pow(nAlpha(Global::A), order + 1)
        * alpha;
    double sum_prob = 0;
    for (UngappedKmer kmer = UngappedKmer::smallestWithLength(order + 1); kmer
        <= lastK; ++kmer) {
      double prob = 1;
      for (int stop = 0; stop < kmer.length(); ++stop) {
        const int start = std::max(0, stop - order + 1);
        const UngappedKmer sub = kmer.subKmer(start, stop);
        prob *= bgCompositionalCond[sub];
      }
      bgCompositionalProb[kmer]
          = static_cast<REAL> ((bgCompositionalProb[kmer]
              + numberOfPseudocounts * prob) / (counts[order]
              + numberOfPseudocounts));
      if (kmer.lastWildcardPosition() == UngappedKmer::undef) sum_prob += bgCompositionalProb[kmer];
    }

    /* marginalize terminal X and compute conditional probabilities */
    const UngappedKmer lastPrefix(UngappedKmer::greatestWithLength(order));
    double sum_marg = 0;
    for (UngappedKmer prefix = UngappedKmer::smallestWithLength(order); prefix
        <= lastPrefix; ++prefix) {
      UngappedKmer kX(prefix.mutatedCopy(prefix.length(), 0));
      if (prefix.lastWildcardPosition() == UngappedKmer::undef) {
        sum_marg += bgCompositionalProb[kX];
      }
      double sum_cond = 0;
      for (uint8_t res = 0; res <= nAlpha(Global::A); ++res) {
        UngappedKmer curr(kX.mutatedCopy(kX.length() - 1, res));
        if (bgCompositionalProb[kX] == 0) {
          bgCompositionalProb[curr] = 0;
          sum_cond = 1;
        } else {
          bgCompositionalCond[curr] = bgCompositionalProb[curr]
              / bgCompositionalProb[kX];
          if (res > 0) sum_cond += bgCompositionalCond[curr];
        }
      }
      if (!(std::abs(sum_cond - 1) < 1e-3)) {
        fprintf(stderr, "P_comp_cond[%s]-1 == %g >= 1e-3\t which is BAD!\n",
            kX.toString().c_str(), sum_cond - 1);
      }
      assert(std::abs(sum_cond-1)<1e-3);
    }
    if (!(std::abs(sum_marg - 1) < 1e-3)) {
      fprintf(stderr, "sum_marg compProb[%d] = %f\n", order, sum_marg);
    }
    assert(std::abs(sum_marg-1)<1e-3);
  }

  delete[] kmer_chars;
  delete[] counts;
}

double NullModelCompositional::getProbability(char const * const kmer,
    const int length) {
  static unsigned char trans[100];
  const int len = (length == -1) ? static_cast<int> (strlen(kmer)) : length;
  for (int i = 0; i < len; ++i) {
    trans[i] = AlphaCode(static_cast<unsigned char>(kmer[i]), Global::A);
  }
  return getProbability(trans, len);
}

double NullModelCompositional::getProbability(unsigned char const * const kmer,
    const int length) {
  double p;
  if (length <= modelOrder + 1) {
    p = bgCompositionalProb[UngappedKmer(kmer, length)];
  } else {
    p = bgCompositionalProb[UngappedKmer(kmer, modelOrder + 1)];
    for (int stop = modelOrder + 1; stop < length; ++stop) {
      p *= bgCompositionalCond[UngappedKmer(kmer + stop - modelOrder,
          modelOrder + 1)];
    }
  }
  return p;
}

double NullModelCompositional::getProbability(const UngappedKmer &kmer) {
  return getProbability(kmer.toString().c_str());
}
