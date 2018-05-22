#include "Kmer.h"
#include "elongationCandidates.h"
#include "../NullModel.h"
#include <cstddef>

std::ostream& operator<<(std::ostream &os, const Kmer &res) {
	if (Global::aa) {
		const Alphabet &alphabet = MProGlobal::getS().alphabet;
		os << res.getKmer()->toString() << endl;
		std::string ks = res.getKmer()->bestAAString();
		os << "\tbest aa representation: " << ks << endl;
//		chr_last = ks[0];
//		os << "\tP(" << chr_last << ") = "
//				<< MProGlobal::getS().negBgProbs[MProGlobal::getS().alphabet.ctoi(
//						chr_last)] << "\t";
//		int num_gaps = 0;
//		for (size_t i = 1; i < ks.length(); ++i) {
//			chr_curr = ks[i];
//			if (chr_curr != alphabet.wildcardChar()) {
//				os << "P(" << chr_curr << "|" << chr_last << "," << num_gaps
//						<< ") = "
//						<< MProGlobal::getS().condProbs[num_gaps][alphabet.ctoi(
//								chr_last)][alphabet.ctoi(chr_curr)] << "   ";
//				chr_last = chr_curr;
//				num_gaps = 0;
//			} else {
//				++num_gaps;
//			}
//		}
		os << "\tP_pos: " << res.p_pos << endl << "\tP_set: " << exp(res.p_set)
				<< endl << "\tSet size: " << res.setSize << endl;
		if (!res.debugString.empty()) {
			os << res.debugString << endl;
		}
		os << "\tSeeds: ";
		size_t count = 0;
		for (MatchContainer::const_iterator m_it = res.seeds.begin(); m_it
				!= res.seeds.end(); ++m_it) {
			os << (count > 0 ? "\t       " : "") << *m_it << " (";
			for (int i = 0; i < static_cast<int>(ks.length()); ++i) {
				if (ks[i] != alphabet.wildcardChar()) {
					os
							<< alphabet.itoc(
									MProGlobal::getS()._posSet[m_it->seq][m_it->pos
											+ i]);
				} else {
					os << ".";
				}
			}
			os << ")";
			if (count < res.seeds.size() - 1) {
				os << endl;
			}
			++count;
		}
	} else {
		//os << "operator<< undefined for Kmer in DNA mode" << endl;
		std::string ks = res.getKmer()->toString(0, ElongCandidates::IUPAC_CHARS);
		os << endl << "IUPAC string: " << ks << endl;
		os << "\tSeeds: ";
		size_t count = 0;
		for (MatchContainer::const_iterator m_it = res.seeds.begin(); m_it
				!= res.seeds.end(); ++m_it) {
			os << (count > 0 ? "\t       " : "") << *m_it << " (";
			for (size_t i = 0; i < ks.length(); ++i) {
				if (ks[i] != 'N') {
					os << AlphaChar(Global::posSet->entity[m_it->seq]->S[0][m_it->pos+i], Global::A);
				} else {
					os << ".";
				}
			}
			os << ")";
			if (count < res.seeds.size() - 1) {
				os << endl;
			}
			++count;
		}
		os << endl;
		os << "nbSeeds: " << res.seeds.size() << endl;
	}
	return os;
}
