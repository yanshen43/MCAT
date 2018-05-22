/*
 * IntervalSet.h
 *
 *  Created on: Oct 28, 2009
 *      Author: eckhart
 */

#ifndef INTERVALSET_H_
#define INTERVALSET_H_

#include <cassert>
#include <iostream>
#include <stdlib.h>
#include <string.h>

class IntervalSet {

public:
	IntervalSet(const int min, const int max) :
		minValue(min), maxValue(max), size(0) {
		isMember = (bool*)calloc(maxValue - minValue + 1, sizeof(bool));
	}

	~IntervalSet() {
		free(isMember);
	}

	void add(const int m) {
		assert(m>=minValue && m<=maxValue);
		if (!isMember[m - minValue]) {
//			std::cout << "adding " << m << " ";
			++size;
			isMember[m - minValue] = true;
		}
	}

	bool remove(const int m) {
		assert(m>=minValue && m<=maxValue);
		if (isMember[m - minValue]) {
			isMember[m - minValue] = false;
			--size;
			return true;
		}
		return false;
	}

	bool contains(const int m) {
		assert(m>=minValue && m<=maxValue);
		return isMember[m - minValue];
	}

	bool containsAll() const {
		return size == maxValue - minValue + 1;
	}

	int getSize() const {
		return size;
	}

	void clear() {
		memset(isMember, 0, (maxValue - minValue + 1) * sizeof(bool));
		size = 0;
	}

private:
	const int minValue;
	const int maxValue;
	int size;
	bool* isMember;

};

#endif /* INTERVALSET_H_ */
