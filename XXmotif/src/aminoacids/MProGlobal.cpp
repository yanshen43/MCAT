#include "MProGlobal.h"

#include "../cs/substitution_matrix-inl.h"
#include "../cs/gonnet_matrix.h"
#include "../cs/blosum_matrix.h"
#include "../cs/profile_library-inl.h"
#include "../cs/utils-inl.h"

#include <iomanip>
#include <limits>
#include <string>
#include <vector>
#include <fstream>

MProGlobal* MProGlobal::instance;

MProGlobal::MProGlobal() {
	S.D_MAX = 7;
	S.alphabet = Alphabet(Global::A);
	S.states = StateLib(S.alphabet);
	/*
	 * convert ugly C structures to something more usable
	 */
	if (Global::negSet != NULL) {
		for (int i = 1; i <= Global::negSet->nent; i++) {
			const e_type &ent = Global::negSet->entity[i];
			std::vector<uint8_t> v((uint8_t*) (ent->S[0] + 1),
					(uint8_t*) (ent->S[0] + ent->n + 1));
			S._negSet.push_back(Sequence(v));
		}
	}
	for (int i = 1; i <= Global::posSet->nent; i++) {
		const e_type &ent = Global::posSet->entity[i];
		std::vector<uint8_t> v((uint8_t*) (ent->S[0] + 1),
				(uint8_t*) (ent->S[0] + ent->n + 1));
		S._posSet.push_back(Sequence(v));
	}
	initStates(Global::type_of_states);
}

MProGlobal::~MProGlobal() {
  delete instance;
}

std::string MProGlobal::getFullResourceName(const std::string &name) {
	std::stringstream s;
	size_t pos = Global::argv0.find("Release/XXmotif");
	if ( pos == std::string::npos ) {
		pos = Global::argv0.find("Debug/XXmotif");
	}
	if ( pos == std::string::npos ) {
		pos = Global::argv0.find("GTest/XXmotif");
	}
	if ( pos != std::string::npos ) {
		std::string path(Global::argv0);
		path.erase(pos);
		s << path << "resources/mpro/" << name;
	} else {
		cerr << "ERROR: MPro resources not found, check your installation" << endl
				<< "argv0 is '" << Global::argv0 << "'" << endl;
		exit(1);
	}
	return s.str();
}

void MProGlobal::initStates(const StateType &type) {
	S.states.setType(type);
	if (type == CS_63) {
		/* build a more comfortable wrapper around Andreas' context profiles */
		std::string cslib(getFullResourceName("elm_cs.lib"));
		FILE * fin = fopen(cslib.c_str(), "r");
		cs::ProfileLibrary<cs::AminoAcid> lib(fin);
		fclose(fin);
		int index = 0;
		for (int i = 0; i < lib.size(); i++) {
			if (i != 8 && i != 58) {
				// remove two "background states":
				// - one because state index 63 is disallowed (unique ids)
				// - one for the start/end state
				S.states.push_back(ColumnState(lib[i], index, S.alphabet));
				index++;
			}
		}
		/* add state for start/end of sequence */
		S.states.push_back(ColumnState(index, S.alphabet,
				Alphabet::startStopChar));
	} else if (StateType::toString(type).substr(0,6) == "BLOSUM") {
		cs::BlosumMatrix
				blo(
						type == BLOSUM45_21 ? cs::BlosumMatrix::BLOSUM45
								: (type == BLOSUM62_21 ? cs::BlosumMatrix::BLOSUM62
										: cs::BlosumMatrix::BLOSUM80));
		for (int i = 0; i < cs::AminoAcid::instance().size(); ++i) {
			cs::ContextProfile<cs::AminoAcid> prof(0, 1);
			prof.TransformToLogSpace();
			for (int j = 0; j < cs::AminoAcid::instance().size(); ++j) {
				prof[0][j] = static_cast<float>(log2(blo.r(j, i)));
			}
			S.states.push_back(ColumnState(prof, i, S.alphabet));
		}
		/* add state for start/end of sequence */
		S.states.push_back(ColumnState(static_cast<int>(S.states.size()), S.alphabet,
				Alphabet::startStopChar));
	} else if (type == GONNET_21) {
		cs::GonnetMatrix go;
		for (int i = 0; i < cs::AminoAcid::instance().size(); ++i) {
			cs::ContextProfile<cs::AminoAcid> prof(0, 1);
			prof.TransformToLogSpace();
			for (int j = 0; j < cs::AminoAcid::instance().size(); ++j) {
				prof[0][j] = static_cast<float>(log2(go.r(j, i)));
			}
			S.states.push_back(ColumnState(prof, i, S.alphabet));
		}
		/* add state for start/end of sequence */
		S.states.push_back(ColumnState(static_cast<int>(S.states.size()), S.alphabet,
				Alphabet::startStopChar));
	} else {
		cerr << "Unknown state type: " << StateType::toString(type) << endl;
		exit(1);
	}

	S.states.sortByEntropy();

	for (int i = 0; i < static_cast<int>(S.states.size()); ++i) {
		double sum = 0;
		for (int j = 0; j < static_cast<int>(S.alphabet.size()); ++j) {
			if ((int) j == S.alphabet.wildcardIndex())
				continue;
			sum += exp(S.states[i][j]);
		}
		if (!(std::abs(1-sum)<1e-3)) {
			fprintf(stderr, "ERROR: |1-sum| %g for state %d (%s)\n", std::abs(1-sum), i, S.states[i].toString().c_str());
		}
		assert(std::abs(1-sum)<1e-3);
	}

	for (int i = 0; i < static_cast<int>(S.states.size()); ++i) {
		cerr << S.states[i].toString(static_cast<int>(S.alphabet.size())) << endl;
	}

	S.states.computeSigStates(Global::aaStateSigThresh);
	for (int i = 0; i < static_cast<int>(S.alphabet.size()); ++i) {
		cout << S.alphabet.itoc(i) << "(" << std::setw(2) << i << "): ";
		for (StateLib::pair_list_type::const_iterator it =
				S.states.sigStates[i].begin(); it
				!= S.states.sigStates[i].end(); ++it) {
			cout << "  (" << (it->first) << " " << S.states[it->second] << "),";
		}
		cout << endl;
	}
}
