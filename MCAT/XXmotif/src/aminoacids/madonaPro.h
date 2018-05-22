#ifndef MADONNAPRO_H_
#define MADONNAPRO_H_

#include "MProGlobal.h"
#include "sequence.h"
#include "../AbstractKmer.h"
#include "../Globals.h"
#include "../SmallKmer.h"
#include "../elongationPhase/Kmer.h"
#include "../elongationPhase/Match.h"
#include "../elongationPhase/elongationCandidates.h"
#include "../memoryPool/pool_alloc.h"

#ifdef GTEST_CONFIGURATION
#include <gtest/gtest.h>
#endif

#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

/** Container for results of one sequence for order statistics */
class SequenceResult {
public:
	SequenceResult() :
		p_pos(1), p_pos_seq(1), p_diso(1), p_set(1) {
	}
	/** p value of best match */
	double p_pos;
	/** p value for whole sequence */
	double p_pos_seq;
	/** p value for disorder/conservation (per instance) */
	double p_diso;
	/** p value for whole set if sequences included up to this one */
	double p_set;
	/** comparison by sequence p values */
	bool operator<(const SequenceResult &other) const {
		if (p_pos_seq < other.p_pos_seq) {
			return true;
		} else {
			return (p_pos_seq == other.p_pos_seq) && (p_diso < other.p_diso);
		}
	}
};

class MadonaPro {
#ifdef GTEST_CONFIGURATION
	FRIEND_TEST(MadonaProTest, aaTrimerStringTest);
	FRIEND_TEST(MadonaProTest, allDimerAndTrimerSubstringsTest);
#endif

public:
	static void findInitialMotifs(elongList &el);
	static void initialize(const int _charBits, const int _gapBits,
			const int _seedSize);
	/**
	 * Calibrate the given kmer results using order statistics. Recalculation of
	 * match scores can be turned off (for initial seeds which already contain
	 * correct match scores).
	 */
	static void
	calibrateMotif(Kmer &kres, bool updateMatchScores = true);

private:
	static bool initialized;

	static long totalNumSeeds;

	/** max size (matches plus wildcards) of initial seed motifs */
	static int seedSize;

	/** contains all match positions in positive set in 3d array [gap1][gap2][index] */
	static std::vector<std::vector<std::vector<MatchContainer> > >
			trimerPositions;

	/** set of initial kmers to be calibrated */
	static std::set<SmallKmer*, KmerPtrOrderCmp> initialKmerSet;

	static std::string aaTrimerString(const int gap1, const int gap2,
			const int c1, const int c2, const int c3);
	static int aaTrimerId(const int c1, const int c2, const int c3);

	/** collect dimer/trimer instances match positions */
	static void gatherAADimerAndTrimerPositions();

	/**
	 * Generates all dimer and trimer submotifs of the given (very restricted)
	 * regular expression (only square brackets are allowed to specify alternative
	 * characters at a single position). The implementation is brute-force and ugly,
	 * but it works.
	 * Useful for tracking a set of seeds for a known motif.
	 */
	static std::set<std::string> allDimerAndTrimerSubstrings(
			const std::unordered_set<std::string> &str);
	static std::set<std::string> allDimerAndTrimerSubstrings_simple(
			const std::string &str);

	/**
	 * Calculate (background corrected) score of a match of the given kmer and sequence.
	 */
	struct KmerScores_t {
		KmerScores_t() :
			left(0), right(0) {
		}
		KmerScores_t(const double sl, const double sr) :
			left(sl), right(sr) {
		}
		double left;
		double right;
	};
	static KmerScores_t scoreKmerMatch(const AbstractKmer &kmer,
			const Sequence &seq, const int &startPos, const int &leftLastIndex);

	struct MatchData {
		MatchData() : p_pos(std::numeric_limits<double>::max()), p_dc(std::numeric_limits<double>::max()),
				p_comb(std::numeric_limits<double>::max()), best_p_dc(std::numeric_limits<double>::max()) {}
		MatchContainer::iterator match; /* points to the seed in the current Kmer */
		MadonaPro::KmerScores_t scores; /* sequence scores for left and right part */
		double p_pos; /* p-value for sequence match */
		double p_dc; /* p-value for disorder/conservation */
		double p_comb; /* combined p-value for sequence score and d/c score */
		double best_p_dc; /* best p_dc for all following elements in list of possible matches */
		bool operator>(const MatchData &other) const {
			return scores.left+scores.right > other.scores.left+other.scores.right;
		}
		struct cmpPtrGreaterScore {
			bool operator()(const std::shared_ptr<MatchData> &s1,
				const std::shared_ptr<MatchData> &s2) {
			return *s1 > *s2;
			}
		};
		struct cmpPtrLessPcomb {
			bool operator()(const std::shared_ptr<MatchData> &s1,
				const std::shared_ptr<MatchData> &s2) {
				return s1->p_comb < s2->p_comb;
			}
		};
	};

	typedef std::list<std::shared_ptr<MatchData>, Pool_alloc<std::shared_ptr<MatchData> > > MatchDataContainer;

	struct TempMatchData {
		TempMatchData() :
			match(0), p_pos(1) {
		}
		Match* match; /* pointer to the "real" match */
		KmerScores_t scores; /* scores for left and right part */
		double p_pos;
		bool operator>(const TempMatchData &other) const {
			return match->score > other.match->score;
		}
	};

	/**
	 * Calibrates an initial kmer by searching for seed matches in order of descending
	 * score until a seed has been found for each sequence and a subsequent call to
	 * calibrateInitialMotif.
	 */
	static void update_min_p_dc(MatchDataContainer &cont, MatchDataContainer::reverse_iterator start);
	static std::shared_ptr<Kmer> calibrateInitialMotif(SmallKmer *stateKmer);

	static double getPPos(const AbstractKmer &kmer, const TempMatchData &match,
			const KmerScores_t &best, const int leftPart);

};

#endif /* MADONNAPRO_H_ */
