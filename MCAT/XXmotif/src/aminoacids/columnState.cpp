/*
 * columnState.cpp
 *
 *  Created on: Oct 23, 2009
 *      Author: eckhart
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "columnState.h"
#include "alphabet.h"
#include "utils.h"

ColumnState::ColumnState(const double A, const double C, const double G, const double T):
	p(5), p_sorted_desc(5){
	assert(!Global::aa);

	p[0] = -1e100;							p_sorted_desc[0] = scoreCharPair_type(p[0], 'N');
	p[1] = ( A == 0 ? -1e100 : log(A) );	p_sorted_desc[1] = scoreCharPair_type(p[1], 'A');
	p[2] = ( C == 0 ? -1e100 : log(C) );	p_sorted_desc[2] = scoreCharPair_type(p[2], 'C');
	p[3] = ( G == 0 ? -1e100 : log(G) );	p_sorted_desc[3] = scoreCharPair_type(p[3], 'G');
	p[4] = ( T == 0 ? -1e100 : log(T) );	p_sorted_desc[4] = scoreCharPair_type(p[4], 'T');

	std::sort(p_sorted_desc.begin(), p_sorted_desc.end(), pairCmp<
			scoreCharPair_type::first_type, scoreCharPair_type::second_type,
			std::greater<scoreCharPair_type::first_type> > ());
}

ColumnState::ColumnState(const cs::ContextProfile<cs::AminoAcid> &prof,
		const int idx, const Alphabet &A) :
	index(idx), p(A.size()), p_sorted_desc(A.size()) {
	assert(prof.logspace());
	// copy corresponding profile entries, translating
	// between the (possibly different) indices of amino acids
	// in alphabet
	for (int i = 0; i < prof.alphabet_size() + 1; ++i) {
		const char &aa = cs::AminoAcid::instance().itoc(i);
		p[A.ctoi(aa)] = prof[0][i] * log(2); // convert from log_2 (Andreas) to ln
		p_sorted_desc[A.ctoi(aa)] = scoreCharPair_type(p[A.ctoi(aa)], aa);
	}
	p[A.ctoi(A.wildcardChar())] = log(std::numeric_limits<double>::min());
	p_sorted_desc[A.ctoi(A.wildcardChar())] = scoreCharPair_type(p[A.ctoi(
			A.wildcardChar())], A.wildcardChar());
	p[A.ctoi(Alphabet::startStopChar)]
			= log(std::numeric_limits<double>::min());
	p_sorted_desc[A.ctoi(Alphabet::startStopChar)] = scoreCharPair_type(
			p[A.ctoi(Alphabet::startStopChar)], Alphabet::startStopChar);
	//			for (size_t i = 0; i < alphabetSize()-1; ++i) {
	//				std::cout << i << ": '" << A.itoc(i) << "' == " << p[i];
	//				std::cout << "  <=====  ";
	//				int oldindex = cs::AminoAcid::instance().ctoi(A.itoc(i));
	//				std::cout << oldindex << ": '" << cs::AminoAcid::instance().itoc(
	//						oldindex) << "' == " << prof[0][oldindex];
	//				std::cout << std::endl;
	//			}
	//			std::cout << A.size()-1 << ": '" << A.itoc(A.size()-1) << "' == " << p[A.size()-1] << std::endl;
	std::sort(p_sorted_desc.begin(), p_sorted_desc.end(), pairCmp<
			scoreCharPair_type::first_type, scoreCharPair_type::second_type,
			std::greater<scoreCharPair_type::first_type> > ());
	calculateEntropy(A);
}

ColumnState::ColumnState(const int idx, const Alphabet &A, const char matchChar) :
	index(idx), p(A.size()), p_sorted_desc(A.size()) {
	for (int i = 0; i < A.size(); ++i) {
		p[i] = log(std::numeric_limits<double>::min());
		p_sorted_desc[i] = scoreCharPair_type(p[i], A.itoc(i));
	}
	p[A.ctoi(matchChar)] = log(1 -(A.size()-1)*std::numeric_limits<double>::min());
	p_sorted_desc[A.ctoi(matchChar)] = scoreCharPair_type(
			p[A.ctoi(matchChar)], matchChar);
	std::sort(p_sorted_desc.begin(), p_sorted_desc.end(), pairCmp<
			scoreCharPair_type::first_type, scoreCharPair_type::second_type,
			std::greater<scoreCharPair_type::first_type> > ());
	calculateEntropy(A);
}

std::string ColumnState::toString(const int numBestCols) const {
	std::stringstream s;
	s << std::setfill('0');
	s << "[" << std::setw(2) << index;
	if (numBestCols > 0) {
		s << ":";
	}
	for (int j = 0; j < numBestCols; ++j) {
		s << nthSignificantChar(j) << std::setw(2) << (int) (exp(
				nthSignificantScore(j)) * 100 + 0.5);
		if (j < numBestCols - 1) {
			s << "|";
		}
	}
	s << "]";
	return s.str();
}

void ColumnState::setIndex(const int i) {
	index = i;
}


void ColumnState::calculateEntropy(const Alphabet &alphabet) {
	entropy = 0;
	for (int a=0; a<alphabet.size(); ++a) {
		if (a==alphabet.wildcardIndex()) continue;
		entropy -= (p[a]/log(2))*exp(p[a]);
	}
}

std::ostream& operator<<(std::ostream &os, const ColumnState &state) {
	os << state.toString();
	return os;
}
