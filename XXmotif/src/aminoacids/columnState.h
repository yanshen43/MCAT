/*
 * columnState.h
 *
 *  Created on: Oct 12, 2009
 *      Author: eckhart
 */

#ifndef MPRO_PROFILE_H_
#define MPRO_PROFILE_H_

#include "../cs/amino_acid.h"
#include "../cs/context_profile.h"

#include "alphabet.h"
#include <vector>

/**
 * Class representing a column profile. It is constructed from
 * Andreas' ContextProfile for a ContextProfile with just one column.
 * Furthermore, indexing starts at 1 for compatibility with Madona,
 * whereas ContextProfile indexing starts at 0.
 */
class ColumnState {

	friend std::ostream& operator<<(std::ostream&, const ColumnState&);

public:
	ColumnState() : index(-1), entropy(0) {
	}

	/**
	 * Construct from cs profile.
	 */
	ColumnState(const cs::ContextProfile<cs::AminoAcid> &prof, const int idx, const Alphabet &A);

	/**
	 * Construct a column state that matches only the given character.
	 */
	 ColumnState(const int idx, const Alphabet &A, const char matchChar);

	/**
	 *  Construct a column state with given nucleotide frequencies
	 */
	ColumnState(const double A, const double C, const double G, const double T);

	size_t alphabetSize() const {
		return p.size();
	}

	const double& operator[](const int i) const {
		return p[i];
	}

	double& operator[](const int i) {
		return p[i];
	}

	double nthSignificantScore(const int n) const {
		return p_sorted_desc[n].first;
	}

	char nthSignificantChar(const int n) const {
		return p_sorted_desc[n].second;
	}

	std::string toString(const int numBestCols=3) const;

	void setIndex(const int i);

	double getEntropy() const {
		return entropy;
	}

	bool operator<(const ColumnState &other) const {
		return entropy < other.entropy;
	}

private:
	int index; // index of state in Globals::states
	std::valarray<double> p;
	double entropy;
	typedef std::pair<double,char> scoreCharPair_type;
	std::vector<scoreCharPair_type> p_sorted_desc;
	void calculateEntropy(const Alphabet &alphabet);
};

#endif /* COLUMNSTATE_H_ */
