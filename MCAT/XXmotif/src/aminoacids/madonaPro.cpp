#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <math.h>
#include <set>
#include <stdint.h>
#include <string>
#include <sstream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ios_base;

#include "madonaPro.h"
#include "alphabet.h"
#include "columnState.h"
#include "IntervalSet.h"
#include "sequence.h"
#include "utils.h"
#include "DisoCons.h"
#include "NNetSupplInfProvider.h"
#include "../branch_and_bound-inl.h"
#include "../AbstractKmer.h"
#include "../alphabet.h"
#include "../Globals.h"
#include "../LogTable.h"
#include "../SmallKmer.h"
#include "../UniversalKmer.h"
#include "../utils.h"
#include "../elongationPhase/Extension.h"
#include "../elongationPhase/ExtensionTable.h"
#include "../memoryPool/pool_alloc.h"
#include "../ProgressBar.h"
#include "../pValCalculation.h"

int MadonaPro::seedSize;
std::vector<std::vector<std::vector<MatchContainer> > > MadonaPro::trimerPositions;
std::set<SmallKmer*, KmerPtrOrderCmp> MadonaPro::initialKmerSet;
long MadonaPro::totalNumSeeds = 0;

bool MadonaPro::initialized = false;

void MadonaPro::initialize(const int _charBits, const int _gapBits, const int _seedSize) {
  assert(!initialized);
	seedSize = _seedSize;
	SmallKmer::init(_charBits, _gapBits);
	if (Global::suppInfMode == SUPP_DISOCONS) {
		Global::supplInf = DisoCons::getInstance();
	} else if (Global::suppInfMode == SUPP_NNET) {
		Global::supplInf = NNetSupplInfProvider::getInstance();
	} else if (Global::suppInfMode == SUPP_NO) {
		// ok
	} else {
		std::cerr << "Unknown supplementary information mode: " << Global::suppInfMode << std::endl;
		exit(1);
	}
	Extension::initialize();
  initialized = true;
}

void MadonaPro::findInitialMotifs(elongList &el) {

  assert(initialized);
//	std::stringstream s;
//	s << "/home/eckhart/temp/" << StateType::toString(Global::type_of_states) << ".tex";
//	S.states.writeTikz(21,s.str());
//	exit(1);

	cout << endl << "collecting dimer and trimer positions... " << std::flush;
	gatherAADimerAndTrimerPositions();
	cout << "done" << endl;
	cout << endl;

	{
    std::set<std::string> trackedMotifs;
    if (Global::trackedOnly) {
      trackedMotifs.insert(Global::trackedMotifs.begin(), Global::trackedMotifs.end());
    } else {
      trackedMotifs = allDimerAndTrimerSubstrings(Global::trackedMotifs);
    }
    Global::trackedMotifs.clear();
    Global::trackedMotifs.insert(trackedMotifs.begin(), trackedMotifs.end());
	}

	typedef std::list<std::shared_ptr<Kmer> > result_list_t;
	result_list_t kmerResults;
	{
		int count = 0;
		if (Global::trackedOnly) {
			for (std::set<SmallKmer*>::const_iterator k_it = initialKmerSet.begin(); k_it != initialKmerSet.end(); ++k_it) {
				if (Global::isTracked((*k_it)->bestAAString())) {
					++count;
				}
			}
		}
		const int total = Global::trackedOnly ? count : static_cast<int>(initialKmerSet.size());
		const int every = std::max(1, total / 1000);
		const double log_ecut = log(Global::extensionECut);
		count = 0;
		cout << "Calibrating " << total << " motifs... " << std::flush;

		typedef std::vector<std::list<int> > aalist_t;
		aalist_t sigAA;
		sigAA.resize(MProGlobal::getS().states.size());
		for (int s = 0; s < static_cast<int>(MProGlobal::getS().states.size()); ++s) {
			for (int aa = 0; aa < MProGlobal::getS().alphabet.size(); ++aa) {
				for (StateLib::pair_list_type::const_iterator sig_it = MProGlobal::getS().states.sigStates[aa].begin(); sig_it
						!= MProGlobal::getS().states.sigStates[aa].end(); ++sig_it) {
					if ((int)s==sig_it->second) {
						std::list<int>::iterator it = sigAA[s].begin();
						while (it!=sigAA[s].end() && MProGlobal::getS().states[s][aa] < MProGlobal::getS().states[s][*it]) ++it;
						sigAA[s].insert(it,aa);
					}
				}
			}
		}

		for (std::set<SmallKmer*>::const_iterator k_it = initialKmerSet.begin(); k_it != initialKmerSet.end(); ++k_it) {
			if (!Global::trackedOnly || Global::isTracked((*k_it)->bestAAString())) {

				std::shared_ptr<Kmer> kres;
        if (Global::supplInf != NULL) {
					kres = std::shared_ptr<Kmer>(new Kmer(*k_it));
					/* gather state seeds from signigicantly matching exact aa seeds */
					const int gap1 = (*k_it)->gapsAfter(0);
					const int gap2 = (*k_it)->numMatches() == 2 ? 0 : (*k_it)->gapsAfter(1);
					const int s1 = (*k_it)->charAt(0);
					const int s2 = (*k_it)->charAt(1);
					for (std::list<int>::const_iterator aa1 = sigAA[s1].begin(); aa1 != sigAA[s1].end(); ++aa1) {
						for (std::list<int>::const_iterator aa2 = sigAA[s2].begin(); aa2 != sigAA[s2].end(); ++aa2) {
							std::list<int> third;
							if ((*k_it)->numMatches() == 2) {
								third.push_back(0);
							} else {
								third.insert(third.begin(), sigAA[(*k_it)->charAt(2)].begin(),
										sigAA[(*k_it)->charAt(2)].end());
							}
							for (std::list<int>::const_iterator aa3 = third.begin(); aa3 != third.end(); ++aa3) {
								const int id = aaTrimerId(*aa1, *aa2, *aa3);
								kres->seeds.insert(kres->seeds.end(), trimerPositions[gap1][gap2][id].begin(),
										trimerPositions[gap1][gap2][id].end());
							}
						}
					}
					calibrateMotif(*kres, true);
				} else {
					kres = calibrateInitialMotif(*k_it);
				}

				/* sort seeds by (seq,pos) offset to make motifs comparable by == */
				kres->seeds.sort(Match::cmpMatches);

				/* Check if result has any chance to be extended. If not,
				 * forget it immediately (needed to save memory for large sequence sets).
				 */
				if (kres->seeds.size() > 0 && kres->setSize >= Global::minCoverage && ((int)kmerResults.size() <= Global::extensionMinCut || kres->p_set <= kmerResults.back()->p_set || kres->p_set <= log_ecut)) {
					result_list_t::iterator it = kmerResults.begin();
					bool present = false;
					for (; it != kmerResults.end() && kres->p_set > (*it)->p_set; ++it) {
						if (AbstractKmer::haveIdenticalGapPattern(*(kres->getKmer()), *((*it)->getKmer()))
								&& (kres->seeds == (*it)->seeds)) {
							present = true;
							break;
						}
					}
					if (present) {
						continue;
					}

					/* We are now at the first element which is at least as bad as the new one.
					 * So, insert new element here.
					 */
					it = kmerResults.insert(it, kres);

					/* Now check the rest of the list contains an equivalent state kmer with
					 * the same seeds (and worse p-val). If so, remove it.
					 */
					++it;
					while (it != kmerResults.end()) {
						if (AbstractKmer::haveIdenticalGapPattern(*(kres->getKmer()), *((*it)->getKmer()))
								&& (kres->seeds == (*it)->seeds)) {
							it = kmerResults.erase(it);
						} else {
							++it;
						}
					}

					/* Check if list is longer than minimum number of extended seeds,
					 * and if the trailing motifs are worse then the extension cutoff.
					 * In this case, remove them.
					 */
					double last_p = kmerResults.back()->p_set;
					if (last_p > log_ecut && (int) kmerResults.size() > Global::extensionMinCut) {
						int num_remove = 0;
						for (result_list_t::reverse_iterator back = kmerResults.rbegin(); back != kmerResults.rend()
								&& (*back)->p_set == last_p; ++back, ++num_remove)
							;
						if ((int) kmerResults.size() - num_remove > Global::extensionMinCut) {
							for (int i = 0; i < num_remove; ++i) {
								kmerResults.pop_back();
							}
						}
					}
				}
				++count;
				if (!Global::batch_mode && (count % every == 0 || count == total)) {
					cout << "\rCalibrating " << total << " motifs... " << progressBar(count, total, 40) << std::flush;
				}
			} else {
				delete *k_it;
			}
		}

		/* check that really sorted by p-val */
		for (result_list_t::const_iterator it = kmerResults.begin(); it != kmerResults.end(); ++it) {
			result_list_t::const_iterator next = it;
			++next;
			if (next != kmerResults.end()) {
				assert((*it)->p_set <= (*next)->p_set);
			}
		}

		cout << kmerResults.size() << " done" << endl << endl;
	}

	cout << "***" << endl << "*** total number of seeds of all calibrated motifs: " << totalNumSeeds << endl << "***" << endl;

	/*
	 * show results
	 */
	if (Global::DEBUG) {
		int rank = 0;
		int with_current_rank = 0;
		double last_p = std::numeric_limits<double>::max();
		for (result_list_t::const_iterator it = kmerResults.begin(); it != kmerResults.end(); ++it) {
			++with_current_rank;
			if (last_p != (*it)-> p_set) {
				rank += with_current_rank;
				last_p = (*it)->p_set;
				with_current_rank = 0;
			}
			if (rank <= Global::extensionMinCut || Global::isTracked((*it)->getKmer()->bestAAString())) {
				cout << "Rank " << std::setw(6) << rank << (Global::isTracked(
            (*it)->getKmer()->bestAAString()) ? " (tracked)" : "") << ": "
            << std::setw(MProGlobal::getS().D_MAX + 2)
            << (*it)->getKmer()->bestAAString() << std::setw(5)
            << (Global::isTracked((*it)->getKmer()->bestAAString()) ? " <<< "
                : "") << (*it)->p_set << endl;
			}
		}
		long sumOfSeeds = 0;
		for (result_list_t::const_iterator it = kmerResults.begin(); it != kmerResults.end(); ++it) {
			sumOfSeeds += (*it)->seeds.size();
		}
		cout << endl << endl;
		cout << std::flush;
		int resCount = 0;
		for (result_list_t::const_iterator it = kmerResults.begin(); it != kmerResults.end(); ++it) {
			++resCount;
			const std::string &ks = (*it)->getKmer()->bestAAString();
			if (Global::isTracked(ks) && resCount <= Global::extensionMinCut) {
        std::cerr << "Rank " << std::setw(6) << resCount << (Global::isTracked(
            ks) ? " (tracked)" : "") << ": " << **it << endl << endl;
      }
		}

	}
	/////////////////////////////////////////////////////////////////////////////// end DEBUG

	cout << endl << endl;

	/* remove too bad / too much */
	{
		result_list_t::iterator first_bad = kmerResults.begin();
		int numExt = 1;
		double last_p = std::numeric_limits<double>::max();
		const double log_ecut = log(Global::extensionECut);
		while (first_bad != kmerResults.end()) {
			if (last_p != (*first_bad)->p_set) {
				if (((*first_bad)->p_set > log_ecut && numExt > Global::extensionMinCut) || numExt
						> Global::extensionMaxCut) {
					break;
				}
			}
			last_p = (*first_bad)->p_set;
			++first_bad;
			++numExt;
		}
		kmerResults.erase(first_bad, kmerResults.end());
	}

	/* extend the rest */
	Extension::graphVizOptions *vizopts = new Extension::graphVizOptions;
	vizopts->showEdgesToLosers = true;
	vizopts->showEdgesToLookups = true;
	std::string grfilename(Global::outputDirectory);
	grfilename += "/extGraph.dot";
	std::ofstream grout(grfilename.c_str());
	vizopts->stream = &grout;
	*(vizopts->stream) << "digraph G {" << endl;
	int numExt = 1;
	cout << "Extending initial motifs... " << std::flush;
	for (result_list_t::const_iterator it = kmerResults.begin(); it != kmerResults.end(); ++it) {
		vizopts->id_suffix = "___" + (*it)->getKmer()->bestAAString();
		std::shared_ptr<Kmer> ptr(*it);
		ptr->enrichment.startRegion = 0;
		ptr->enrichment.endRegion = Global::posSet->max_leng;
		elongList extList;
		extList.push_back(ptr);
		if (Global::DEBUG) {
			cout << "... " << (*it)->getKmer()->bestAAString() << " "
          << (*it)->getKmer()->toString() << (Global::isTracked(
          (*it)->getKmer()->bestAAString()) ? " (tracked)" : "") << " => "
          << std::flush;
		}
		extList = Extension::getBestExtensions(extList, 0, vizopts, false);
		el.insert(el.end(), extList.begin(), extList.end());
		if (Global::DEBUG > 0) {
			int i = 1;
			for (result_list_t::const_iterator e_it = extList.begin(); e_it != extList.end(); ++e_it) {
				if (i > 1) cout << "; ";
				cout << (*e_it)->getKmer()->bestAAString() << " " << (*e_it)->getKmer()->toString();
				++i;
			}
			cout << endl;
		} else {
			if (!Global::batch_mode) {
				cout << "\rExtending initial motifs... " << progressBar(numExt, static_cast<int>(kmerResults.size()), 40) << std::flush;
			}
		}
		++numExt;
	}
	cout << endl;
	*(vizopts->stream) << "}" << endl;
	grout.close();
	delete vizopts;

	/* keep only best for identical motifs with the same seed */
	el.sort(Kmer::lessPtr);
	for (elongList::iterator it = el.begin(); it != el.end(); ++it) {
		elongList::iterator next = it;
		++next;
		while (next != el.end()) {
			if (AbstractKmer::haveIdenticalGapPattern(*((*next)->getKmer()), *((*it)->getKmer())) && ((*next)->seeds
					== (*it)->seeds)) {
				next = el.erase(next);
			} else {
				++next;
			}
		}
	}

	cout << endl << "Extended " << (numExt - 1) << " motifs to " << el.size() << " unique motifs" << endl;
	int outCount = 1;
	for (std::list<std::shared_ptr<Kmer> >::const_iterator it = el.begin(); it != el.end(); ++it, ++outCount) {
		if (outCount <= (Global::DEBUG > 0 ? Global::extensionMinCut : 5)) {
      cout << endl << "Rank " << outCount << (Global::isTracked(
          (*it)->getKmer()->bestAAString()) ? " (tracked)" : "") << endl
          << **it << endl;
    }
	}

	trimerPositions.clear();
	initialKmerSet.clear();
	ExtensionTable::freeTable();

}

void MadonaPro::calibrateMotif(Kmer &kres, bool updateMatchScores) {
	if (updateMatchScores) {

		/* calculate kmer split points for branch & bound */
		const int kmerMatches = kres.getKmer()->numMatches();
		const int leftPart = kmerMatches > 5 ? 5 : kmerMatches; // number of characters in left part of kmer split
		const int rightPart = kmerMatches > 5 ? kmerMatches - leftPart : 0;
		assert(leftPart + rightPart == kmerMatches);

		/* calculate best possible match scores */
		KmerScores_t bestScores;
		for (int i = 1; i <= kmerMatches; ++i) {
			double &score = (i <= leftPart) ? bestScores.left : bestScores.right;
			score += MProGlobal::getS().states[kres.getKmer()->charAt(i - 1)].nthSignificantScore(0)
							- Global::negBg_log[MProGlobal::getS().alphabet.ctoi(
							    MProGlobal::getS().states[kres.getKmer()->charAt(i - 1)].nthSignificantChar(
										0))];
		}

		/*
		 * find worst match
		 */
		double worstScore_left = std::numeric_limits<double>::max();
		double worstScore_right = std::numeric_limits<double>::max();
		for (MatchContainer::iterator mit = kres.seeds.begin(); mit != kres.seeds.end();) {
			const KmerScores_t s = scoreKmerMatch(*(kres.getKmer()), MProGlobal::getS()._posSet[mit->seq], mit->pos, leftPart - 1); // -1 because scoreKmerMatch character indices start at 0
			if (s.left < -500 || s.right < -500) {
				mit = kres.seeds.erase(mit);
				continue;
			}
			if (s.left < worstScore_left) {
				worstScore_left = s.left;
			}
			if (s.right < worstScore_right) {
				worstScore_right = s.right;
			}
			++mit;
		}

		/*
		 * calculate P-values of matches
		 */
		Kmer_Generator_Reloaded<aa_hash_t> &kgen = Kmer_Generator_Reloaded<aa_hash_t>::getInstance();
		kgen.reinitialize(*(kres.getKmer()), leftPart, worstScore_left - (bestScores.right
				- worstScore_right), worstScore_right - (bestScores.left - worstScore_left));
		for (MatchContainer::iterator mit = kres.seeds.begin(); mit
        != kres.seeds.end(); ++mit) {
      double p = kgen.get_kmer_probability_matched_only(
          kres.getKmer()->matchedCharsStartingAt(
              MProGlobal::getS()._posSet[mit->seq], mit->pos));
      p *= PVal_Calculator::getInstance().getOverrepCorrection(
          kres.getKmer()->numMatches());
      mit->score = static_cast<float> (p);
    }

	}

	/*
	 * find best matches in each sequence
	 */
	std::vector<SequenceResult> seqResults(MProGlobal::getS()._posSet.size());
	for (MatchContainer::const_iterator mit = kres.seeds.begin(); mit != kres.seeds.end(); ++mit) {
		if (mit->score < seqResults[mit->seq].p_pos) {
			seqResults[mit->seq].p_pos = mit->score;
			if (Global::supplInf != NULL) {
				seqResults[mit->seq].p_diso = Global::supplInf->getPvalue(kres.getKmer(), mit->seq, mit->pos);
			}
		}
	}
	/*
	 * calculate per-sequence P-values
	 */
	for (size_t i = 0; i < seqResults.size(); ++i) {
		seqResults[i].p_pos_seq = 1 - pow((1-seqResults[i].p_pos), Global::posSet->avgLength[Global::posSet->nent] - kres.getKmer()->length() + 1);
	}
	/*
	 * order statistic
	 */
	std::sort(seqResults.begin(), seqResults.end());
//	for (size_t i = 0; i < seqResults.size(); ++i) {
//		fprintf(stderr, "[%2lu]p_pos: %g\tp_pos_seq: %g\tp_diso_seq: %g\t p_set: %g\n", i, seqResults[i].p_pos,
//				seqResults[i].p_pos_seq, seqResults[i].p_diso_seq, seqResults[i].p_set);
//	}
	kres.p_set = std::numeric_limits<double>::max();
	double sum_log_diso_seq = 0;
	const double log_bonferroni = calculate_log_bonferonni(kres.getKmer()->getMotifColumns(), LogTable::LOG_i[Global::neff_discrete]);
	for (int i = 0; i < static_cast<int>(seqResults.size()); ++i) {
		double log_p = calculateOrderStatisticsPvalue(i + 1,
				static_cast<int>(seqResults.size()),
				seqResults[i].p_pos_seq);
		if (Global::supplInf != NULL) {
			sum_log_diso_seq += log(seqResults[i].p_diso);
			const double log_p_diso = PVal_Calculator::getCombinedPval_log(
					sum_log_diso_seq, i+1);
			log_p = PVal_Calculator::getWeightedPval_log(log_p,
					log_p_diso, Global::disoconsWeight);
		}
		log_p += log_bonferroni;
		seqResults[i].p_set = log_p;
		if (seqResults[i].p_set <= kres.p_set) {
			kres.p_pos = seqResults[i].p_pos;
			kres.p_set = log_p;
			kres.setSize = i + 1;
		}
	}
	/*
	 * remove seeds which are no longer significant
	 */
	MatchContainer::iterator m_it = kres.seeds.begin();
	const double thresh = kres.p_pos;
	while (m_it != kres.seeds.end()) {
		if (m_it->score > thresh) {
			m_it = kres.seeds.erase(m_it);
		} else {
			++m_it;
		}
	}
}

inline std::shared_ptr<Kmer> MadonaPro::calibrateInitialMotif(SmallKmer *stateKmer) {
	std::shared_ptr<Kmer> result(new Kmer(stateKmer));

	static unsigned char kmer[3];

	/*
	 * Generate all aa kmers for current state kmer with scores and p-values.
	 */
	double max_score = -std::numeric_limits<double>::max();
	double thresh = 0;
	for (int i = 0; i < stateKmer->numMatches(); ++i) {
		const char chr = stateKmer->charAt(i);
		const double score = MProGlobal::getS().states[chr][MProGlobal::getS().alphabet.ctoi(MProGlobal::getS().states[chr].nthSignificantChar(0))]
		                             - Global::negBg_log[MProGlobal::getS().alphabet.ctoi(MProGlobal::getS().states[chr].nthSignificantChar(0))];
		thresh += score;
		if (score > max_score) {
			max_score = score;
		}
	}
	thresh -= max_score;

	Kmer_Generator_Reloaded<aa_hash_t> &kgen = Kmer_Generator_Reloaded<aa_hash_t>::getInstance();
	kgen.reinitialize(*stateKmer, 5, thresh, 0);
	const Kmer_Generator_Reloaded<aa_hash_t>::mapType scoreMap = kgen.getScoreMap();

	/*
	 * Iterate over kmers in order of descending score.
	 * Break if a match has been found in all sequences, also
	 * in case of multiple-occurrence-per-sequence model.
	 * Reason in case of mops: it is highly unlikely that all instances
	 * found until this point will be deemed significant, so the rest
	 * can be skipped safely.
	 */
	IntervalSet matchedSeqs(0, static_cast<int>(MProGlobal::getS()._posSet.size()) - 1);
	const int gap1 = stateKmer->gapsAfter(0);
	const int gap2 = stateKmer->numMatches() == 3 ? stateKmer->gapsAfter(1) : 0;
	for (Kmer_Generator_Reloaded<aa_hash_t>::mapType::const_reverse_iterator mapIt = scoreMap.rbegin(); (!matchedSeqs.containsAll())
			&& (mapIt != scoreMap.rend()); ++mapIt) {
		/*
		 * For each score, process the associated list of kmers.
		 */
		for (Kmer_Generator_Reloaded<aa_hash_t>::id_list_type::const_iterator klist_it = mapIt-> second.begin(); klist_it
				!= mapIt->second.end(); ++klist_it) {
			/*
			 * If the kmer has any hits in the _posSet, add it to the list of seeds.
			 */
			uint64_t id = *klist_it;
			kmer[0] = static_cast<uint8_t>((id & kgen.getChrMask()) + 1);
			id >>= kgen.getIdShift();
			kmer[1] = static_cast<uint8_t>((id & kgen.getChrMask()) + 1);
			id >>= kgen.getIdShift();
			kmer[2] = static_cast<uint8_t>(stateKmer->numMatches() == 3 ? (id & kgen.getChrMask()) + 1 : 0);
			const int index = aaTrimerId(kmer[0], kmer[1], kmer[2]);
			const double pv = kgen.get_kmer_probability_matched_only(kmer);
			/*
			 * Add all sequences matching the kmer to the list of matched seqs.
			 */
			for (MatchContainer::const_iterator match = trimerPositions[gap1][gap2][index].begin(); match
					!= trimerPositions[gap1][gap2][index].end(); ++match) {
				matchedSeqs.add(match->seq);
				/*
				 * Remember match as a possible seed. The seeds are
				 * added in descending P-value order, so the seeds will
				 * initially be sorted (if disocons is not in use).
				 */
				Match m(*match);
				m.score = static_cast<float>(pv);
				result->seeds.push_back(m);
			}
		}
	}
	calibrateMotif(*result, false);
	return result;
}

/*
 * Helper functions
 *
 *
 *
 *
 */

int MadonaPro::aaTrimerId(const int c1, const int c2, const int c3) {
	return (c1 * MProGlobal::getS().alphabet.size() + c2) * MProGlobal::getS().alphabet.size() + c3;
}

std::string MadonaPro::aaTrimerString(const int gap1, const int gap2, const int c1, const int c2, const int c3) {
	std::stringstream str;
	str << MProGlobal::getS().alphabet.itoc(c1);
	for (int i = 0; i < gap1; ++i) {
		str << MProGlobal::getS().alphabet.wildcardChar();
	}
	str << MProGlobal::getS().alphabet.itoc(c2);
	if (c3 != MProGlobal::getS().alphabet.wildcardIndex()) {
		for (int i = 0; i < gap2; ++i) {
			str << MProGlobal::getS().alphabet.wildcardChar();
		}
		str << MProGlobal::getS().alphabet.itoc(c3);
	}
	return str.str();
}

/*
 * store all trimer positions in array [gap1][gaps][index]
 */
void MadonaPro::gatherAADimerAndTrimerPositions() {
	trimerPositions.resize(seedSize - 1);
	const size_t num = static_cast<size_t>(pow(MProGlobal::getS().alphabet.size(), 3));
	for (int i = 0; i < static_cast<int>(trimerPositions.size()); ++i) {
		trimerPositions[i].resize(seedSize - 2);
		for (int j = 0; j < static_cast<int>(trimerPositions[i].size()); ++j) {
			trimerPositions[i][j].resize(num);
		}
	}

	const std::vector<Sequence> &set = MProGlobal::getS()._posSet;
	for (uint32_t seq = 0; seq < set.size(); ++seq) {
		for (int32_t pos = 0; pos < static_cast<int>(set[seq].size()); ++pos) {
			if (set[seq][pos] == MProGlobal::getS().alphabet.wildcardIndex()) {
				continue;
			}
			// dimers
			for (int gap1 = 0; gap1 <= MProGlobal::getS().D_MAX && gap1 + 2 <= seedSize && pos + gap1 + 1 < static_cast<int>(set[seq].size()); ++gap1) {
				bool masked = false;
				for (int p = pos + 1; p <= pos + gap1 + 1; ++p) {
					if (set[seq][p] == MProGlobal::getS().alphabet.wildcardIndex()) {
						masked = true;
					}
				}
				if (masked) {
					continue;
				}
				assert(!(set[seq][pos]==0||set[seq][pos+gap1+1]==0));
				trimerPositions[gap1][0][aaTrimerId(set[seq][pos], set[seq][pos + gap1 + 1], 0)].push_back(Match(seq,
						pos, 1));
				if (trimerPositions[gap1][0][aaTrimerId(set[seq][pos], set[seq][pos + gap1 + 1], 0)].size() == 1) {
					//	calibrateMotif all state kmers with score >= thresh at all positions
					for (StateLib::pair_list_type::const_iterator it1 = MProGlobal::getS().states.sigStates[set[seq][pos]].begin(); it1
							!= MProGlobal::getS().states.sigStates[set[seq][pos]].end(); ++it1) {
						for (StateLib::pair_list_type::const_iterator it2 =
						    MProGlobal::getS().states.sigStates[set[seq][pos + gap1 + 1]].begin(); it2
								!= MProGlobal::getS().states.sigStates[set[seq][pos + gap1 + 1]].end(); ++it2) {
							SmallKmer *k = new SmallKmer(gap1, it1->second, it2->second);
							if (initialKmerSet.find(k) == initialKmerSet.end()) {
								initialKmerSet.insert(k);
							} else {
								delete k;
							}
						}
					}
				}
			}
			// trimers
			for (int gap1 = 0; gap1 <= MProGlobal::getS().D_MAX && gap1 + 3 <= seedSize; ++gap1) {
				for (int gap2 = 0; gap2 <= MProGlobal::getS().D_MAX && gap2 + gap1 + 3 <= seedSize && pos + gap1 + gap2 + 2
						< static_cast<int>(set[seq].size()); ++gap2) {
					bool masked = false;
					for (int p = pos + 1; p <= pos + gap1 + gap2 + 2; ++p) {
						if (set[seq][p] == MProGlobal::getS().alphabet.wildcardIndex()) {
							masked = true;
						}
					}
					if (masked) {
						continue;
					}
					assert(!(set[seq][pos]==0||set[seq][pos+gap1+1]==0||set[seq][pos+gap1+gap2+2]==0));
					int index = aaTrimerId(set[seq][pos], set[seq][pos + gap1 + 1], set[seq][pos + gap1 + gap2 + 2]);
					trimerPositions[gap1][gap2][index].push_back(Match(seq, pos, 1));
					if (trimerPositions[gap1][gap2][index].size() == 1) {
						for (StateLib::pair_list_type::const_iterator it1 = MProGlobal::getS().states.sigStates[set[seq][pos]].begin(); it1
								!= MProGlobal::getS().states.sigStates[set[seq][pos]].end(); ++it1) {
							for (StateLib::pair_list_type::const_iterator it2 = MProGlobal::getS().states.sigStates[set[seq][pos + gap1
									+ 1]].begin(); it2 != MProGlobal::getS().states.sigStates[set[seq][pos + gap1 + 1]].end(); ++it2) {
								for (StateLib::pair_list_type::const_iterator it3 = MProGlobal::getS().states.sigStates[set[seq][pos
										+ gap1 + gap2 + 2]].begin(); it3 != MProGlobal::getS().states.sigStates[set[seq][pos + gap1
										+ gap2 + 2]].end(); ++it3) {
									SmallKmer *k = new SmallKmer(gap1, gap2, it1->second, it2->second, it3->second);
									if (initialKmerSet.find(k) == initialKmerSet.end()) {
										initialKmerSet.insert(k);
									} else {
										delete k;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

std::set<std::string> MadonaPro::allDimerAndTrimerSubstrings_simple(const std::string &str) {
	std::set<std::string> substrings;
	for (size_t pos = 0; pos + 1 < str.length(); ++pos) {
		if (str[pos] != MProGlobal::getS().alphabet.wildcardChar()) {
			std::string( s1);
			s1 += str[pos];
			for (size_t gap1 = 0; pos + gap1 + 1 < str.length(); ++gap1) {
				if (str[pos + gap1 + 1] != MProGlobal::getS().alphabet.wildcardChar()) {
					std::string s2(s1);
					for (size_t i = 0; i < gap1; ++i) {
						s2 += MProGlobal::getS().alphabet.wildcardChar();
					}
					s2 += str[pos + gap1 + 1];
					substrings.insert(s2);
					for (size_t gap2 = 0; pos + gap1 + gap2 + 2 < str.length(); ++gap2) {
						if (str[pos + gap1 + gap2 + 2] != MProGlobal::getS().alphabet.wildcardChar()) {
							std::string s3(s2);
							for (size_t i = 0; i < gap2; ++i) {
								s3 += MProGlobal::getS().alphabet.wildcardChar();
							}
							s3 += str[pos + gap1 + gap2 + 2];
							substrings.insert(s3);
						}
					}
				}
			}
		}
	}
	return substrings;
}

std::set<std::string> MadonaPro::allDimerAndTrimerSubstrings(const std::unordered_set<std::string> &strs) {
  std::set<std::string> ret;
  for (auto s = strs.begin(); s != strs.end(); ++s) {
    const std::string& str = *s;
    /*
     * collect possible characters for each position in the result string
     */
    std::list<std::list<char> > charsAtPos;
    for (size_t strpos = 0; strpos < str.length(); ++strpos) {
      std::list<char> currentChars;
      if (str[strpos] == '[') {
        while (str[++strpos] != ']') {
          currentChars.push_back(str[strpos]);
          if (strpos >= str.length()) {
            cerr << "Syntax error in regex " << str << endl;
            exit(1);
          }
        }
      } else {
        currentChars.push_back(str[strpos]);
      }
      charsAtPos.push_back(currentChars);
    }
    /*
     * build all combinations by iterating over the positions and appending
     * every possible next char at each string built so far
     */
    std::list<std::string> stringsSoFar;
    stringsSoFar.push_back(std::string(""));
    for (std::list<std::list<char> >::const_iterator cpit = charsAtPos.begin(); cpit != charsAtPos.end(); ++cpit) {
      std::list<std::string> newStrings;
      for (std::list<std::string>::const_iterator sit = stringsSoFar.begin(); sit != stringsSoFar.end(); ++sit) {
        for (std::list<char>::const_iterator cit = cpit->begin(); cit != cpit->end(); ++cit) {
          std::stringstream s;
          s << *sit << *cit;
          newStrings.push_back(s.str());
        }
      }
      stringsSoFar = newStrings;
    }
    /*
     * now build all substrings with wildcards for each string matching the initial pattern
     */
    for (std::list<std::string>::const_iterator sit = stringsSoFar.begin(); sit != stringsSoFar.end(); ++sit) {
      std::set<string> subset = allDimerAndTrimerSubstrings_simple(*sit);
      ret.insert(subset.begin(), subset.end());
    }
  }
	return ret;
}

MadonaPro::KmerScores_t MadonaPro::scoreKmerMatch(const AbstractKmer &kmer, const Sequence &seq, const int &startPos,
		const int &lastLeftIndex) {
	double score_left = 0;
	double score_right = 0;
	int pos = startPos;
	//	for (SmallKmer::KmerConstIterator it = kmer.begin(); it != kmer.end(); ++it) {
	//	score += states[it->chr][seq[pos]];
	int i = 0;
	while (i < kmer.numMatches() - 1) {
		double &score = i <= lastLeftIndex ? score_left : score_right;
		// no need to check for wildcard since in that case
		// -infinity is the log score, and the sum of all
		// scores will also be -infinity
		//			if (s[pos] == wildcardIndex ) {
		//			}

		score += MProGlobal::getS().states[kmer.charAt(i)][seq[pos]] - Global::negBg_log[seq[pos]];
		pos += 1 + kmer.gapsAfter(i);
		++i;
	}
	double &score = i <= lastLeftIndex ? score_left : score_right;
	score += MProGlobal::getS().states[kmer.charAt(i)][seq[pos]] - Global::negBg_log[seq[pos]];
	return KmerScores_t(score_left, score_right);
}
