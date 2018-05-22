#include "DisoCons.h"
#include "../Globals.h"
#include "../AbstractKmer.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using std::cerr;
using std::cout;
using std::endl;

DisoCons* DisoCons::instance;

AbstractSupplementaryInformationProvider const* DisoCons::getInstance() {
	if (instance == NULL) {
		instance = new DisoCons(MProGlobal::getS());
	}
	return instance;
}

DisoCons::~DisoCons() {
  /* delete statement commented out due to OS X compile error/warning */
//	delete instance;
}

DisoCons::DisoCons(const MProGlobal::SeqInfo_t &S) {
	std::string pvfilename(MProGlobal::getFullResourceName("p_disocons.dat"));
	std::ifstream pvfile(pvfilename.c_str());
	if (!pvfile.is_open()) {
		cerr << "Can't read disocons P values from " << pvfilename << endl;
		exit(1);
	} else {
		Pdc.resize(dcMatrixRes + 1);
		for (int i = 0; i <= dcMatrixRes; ++i) {
			Pdc[i].resize(dcMatrixRes + 1);
		}
		while (!pvfile.eof()) {
			std::string line;
			getline(pvfile, line);
			std::stringstream s(line);
			int i;
			int j;
			s >> i;
			s >> j;
			s >> Pdc[i][j];
		}
	}
	for (int i = 0; i <= dcMatrixRes; ++i) {
		for (int j = 0; j <= dcMatrixRes; ++j) {

		}
	}
	std::string dcfilename(Global::name);
	dcfilename = getBasename(dcfilename.c_str());
	dcfilename += ".dc";
	readDCFile(dcfilename, S._posSet, diso, cons);
}

void DisoCons::readDCFile(const std::string &dcfilename, const std::vector<Sequence> &sset, vecvec_t &dis, vecvec_t &con) {
	std::ifstream dcfile(dcfilename.c_str());
	if (!dcfile.is_open()) {
		cerr << "Can't open disorder/conservation file " << dcfilename << endl;
		exit(1);
	} else {
		dis.resize(sset.size());
		con.resize(sset.size());
		int seq = -1;
		while (!dcfile.eof()) {
			/*
			 * format of discons file: repetition of
			 * 1st line: > id
			 * 2nd line: (seqlen) disorder values, separated by a space
			 * 3rd line: (seqlen) conservation values
			 */
			std::string line;
			getline(dcfile, line);
			if (line.length() == 0 || line[0] == '>') {
				++seq;
				continue;
			} else {
				const int termOffset = Global::termMode==BOTH || Global::termMode==POS ? 1 : 0;
				dis[seq].resize(sset[seq].size());
				con[seq].resize(sset[seq].size());

				/* read disorder */
				std::stringstream s(line);
				double v;
				int i = termOffset;
				while (s >> v) dis[seq][i++] = v;
				if (!(i==(int)sset[seq].size() - termOffset)) {
					fprintf(stderr, "Index is %d, should be %d-%d==%d\n", i, (int)sset[seq].size(), termOffset, (int)sset[seq].size() - termOffset);
					assert(i==(int)sset[seq].size() - termOffset);
				}
				
				/* read conservation */
				getline(dcfile, line);
				std::stringstream t(line);
				i = termOffset;
				while (t >> v) con[seq][i++] = v;
				assert(i==(int)sset[seq].size() - termOffset);
			}
		}
		dcfile.close();
	}
}


double DisoCons::getPvalue(AbstractKmer const* const kmer, const int &seq,
		const int &startPos) const {
	const std::vector<Sequence> &sset = MProGlobal::getS()._posSet;
	double d = 0;
	double c = 0;
	int pos = startPos;
	int ignored = 0;
	for (int i=0; i < kmer->numMatches()-1; ++i) {
		if (MProGlobal::getS().alphabet[sset[seq][pos]] == MProGlobal::getS().alphabet.startStopChar) {
			/* ignore sequence ends (no values for $) */
			++ignored;
		} else {
			d += diso[seq][pos];
			c += cons[seq][pos];
		}
		pos += 1 + kmer->gapsAfter(i);
	}
	if (MProGlobal::getS().alphabet[sset[seq][pos]] == MProGlobal::getS().alphabet.startStopChar) {
		/* ignore sequence ends (no values for $) */
		++ignored;
	} else {
		d += diso[seq][pos];
		c += cons[seq][pos];
	}
	d /= kmer->numMatches() - ignored;
	c /= kmer->numMatches() - ignored;
	double P = Pdc[(int) (d * 100 + 0.5)][(int) (c * 100 + 0.5)];
	return P;
}

double DisoCons::getPvalue(motif_columns_type const* const cols, const int &seq,
		const int &startPos) const {
	const std::vector<Sequence> &sset = MProGlobal::getS()._posSet;
	const int seq_pro = seq - 1;
	double d = 0;
	double c = 0;
	int ignored = 0; // number of $ characters in matched seq
	const int offset = -cols->front();
	for (motif_columns_type::const_iterator it = cols->begin(); it != cols->end(); ++it) {
		const int pos_pro = startPos-1+*it+offset;  // Holger's indices start at 1
		if  (seq_pro>=(int)diso.size() || pos_pro >= (int)diso[seq_pro].size()) {
			fprintf(
					stderr,
					"ERROR: index: diso[%d][%d]=%f\tsize(dis)=%d\tsize(dis[%d])=%d\tid=%s\t|%s[%d]|=%d\n",
					seq_pro, pos_pro, diso[seq_pro][pos_pro],
					static_cast<int>(diso.size()), seq_pro,
					static_cast<int>(diso[seq_pro].size()),
					Global::negSet->entity[seq]->info[0],
					"posSet", seq,
					Global::posSet->entity[seq]->n);
		}
		if (MProGlobal::getS().alphabet[sset[seq_pro][pos_pro]] == MProGlobal::getS().alphabet.startStopChar) {
			/* ignore sequence ends (no values for $) */
			++ignored;
		} else {
			d += diso[seq_pro][pos_pro];
			c += cons[seq_pro][pos_pro];
		}
	}
	d /= static_cast<int>(cols->size()) - ignored;
	c /= static_cast<int>(cols->size()) - ignored;
	double P = Pdc[(int) (d * 100 + 0.5)][(int) (c * 100 + 0.5)];
	return P;
}
