#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "NNetSupplInfProvider.h"
#include "NNet.h"
#include "../Globals.h"
#include "../AbstractKmer.h"
#include "../pValCalculation.h"

NNetSupplInfProvider* NNetSupplInfProvider::instance;
NNet* NNetSupplInfProvider::nn;

AbstractSupplementaryInformationProvider const* NNetSupplInfProvider::getInstance() {
  if (instance == NULL) {
    instance = new NNetSupplInfProvider(MProGlobal::getS());
  }
  return instance;
}

NNetSupplInfProvider::~NNetSupplInfProvider() {
  /* delete statements commented out due to OS X compile error/warning */
//  delete instance;
//  delete nn;
}

/* initialize the provider instance using the given sequence set */
NNetSupplInfProvider::NNetSupplInfProvider(const MProGlobal::SeqInfo_t &S) {
  /* initialize neural network */
  nn = new NNet(Global::nnetFilename);

  /* read precalculated P-value table */
  std::string pfilename(Global::nnetFilename);
  std::ifstream pfile(pfilename.c_str());
  if (!pfile.is_open()) {
    std::cerr << "Can't open " << pfilename << std::endl;
    exit(1);
  } else {
    ptable.resize(RES + 1);
    log_ptable.resize(RES + 1);
    std::string line;
    /* skip all lines until P-value table */
    while (!pfile.eof()) {
      getline(pfile, line);
      if (line.find("##### P-values") != std::string::npos) break;
    }
    while (!pfile.eof()) {
      getline(pfile, line);
      if (line.empty()) continue;
      std::stringstream s(line);
      double i;
      s >> i;
      int idx = static_cast<int> (i * RES);
      s >> ptable[idx];
      log_ptable[idx] = log(ptable[idx]);
    }
    pfile.close();
  }

  /* read property file */
  const std::vector<Sequence> &ps = MProGlobal::getS()._posSet;
  std::string propfilename(Global::name);
  propfilename = getBasename(propfilename.c_str());
  propfilename += ".prop";
  std::ifstream propfile(propfilename.c_str());
  if (!propfile.is_open()) {
    cerr << "Can't open property file " << propfilename << std::cerr;
    exit(1);
  } else {
    props.resize(ps.size());
    int seq = -1;
    int num_prop = 0;
    std::string line;
    while (!propfile.eof()) {
      /* format of each entry in files is FASTA-like
       * one line: "> id"
       * then: NUM_PROPS lines of properties,
       *       values on each line are space separated
       */
      getline(propfile, line);
      if (line.empty()) {
        continue;
      } else if (line[0] == '>') {
        if (!(seq == -1 || num_prop == NUM_PROPS)) {
          fprintf(stderr, "seq==%d   num_prop==%d\n", seq, num_prop);
          assert(seq==-1 || num_prop == NUM_PROPS);
        }
        ++seq;
        props[seq].resize(ps[seq].size());
        for (size_t i = 0; i < props[seq].size(); ++i) {
          props[seq][i].resize(NUM_PROPS);
        }
        num_prop = 0;
        continue;
      } else {
        /* whether or not has '$' been added to the initial sequences? */
        const int termOffset = Global::termMode == BOTH || Global::termMode
            == POS ? 1 : 0;
        std::stringstream s(line);
        double v;
        int idx = termOffset;
        while (s >> v) {
          props[seq][idx++][num_prop] = v;
        }
        if (!(idx == (int) ps[seq].size() - termOffset)) {
          fprintf(stderr, "Index is %d, should be %d-%d==%d\n", idx,
              static_cast<int> (ps[seq].size()), termOffset,
              static_cast<int> (ps[seq].size() - termOffset));
          assert(idx==(int)ps[seq].size() - termOffset);
        }
        ++num_prop;
      }
    }
    propfile.close();
  }
}

double NNetSupplInfProvider::getPvalue(AbstractKmer const* const kmer,
    const int &seq, const int &startPos) const {

  /* get a reference to the matched sequence */
  const std::vector<Sequence> &ps = MProGlobal::getS()._posSet;

  /* Starting at startPos, add the log of P-values for each match
   * position, skipping the gaps.
   * No values are available for the termini ('$'), so these
   * are ignored.
   */
  int pos = startPos;
  int ignored = 0; // counter for number of $ characters in matched seq
  double sum_log_p = 0;
  for (int i = 0; i < kmer->numMatches() - 1; ++i) {
    if (MProGlobal::getS().alphabet[ps[seq][pos]]
        == MProGlobal::getS().alphabet.startStopChar) {
      ++ignored;
    } else {
      const int pred_index = static_cast<int> (RES * nn->predict(
          props[seq][pos]));
      sum_log_p += log_ptable[pred_index];
    }
    pos += 1 + kmer->gapsAfter(i);
  }
  if (MProGlobal::getS().alphabet[ps[seq][pos]]
      == MProGlobal::getS().alphabet.startStopChar) {
    ++ignored;
  } else {
    const int pred_index =
        static_cast<int> (RES * nn->predict(props[seq][pos]));
    assert(pred_index >= 0 && pred_index <= RES);
    sum_log_p += log_ptable[pred_index];
  }

  /* The final result is the geometric mean of P-values,
   * i.e. the arithmetic mean in logspace.
   */
  const double log_P_combined = sum_log_p / (kmer->numMatches() - ignored);
  return exp(log_P_combined);
}

/* This is completely analogous to the kmer case.
 * Caveat: DNA numbering starts at 1, not 0 as for proteins,
 * for historical reasons.
 */
double NNetSupplInfProvider::getPvalue(motif_columns_type const* const cols,
    const int &seq, const int &startPos) const {
  const std::vector<Sequence> &ps = MProGlobal::getS()._posSet;
  const int seq_pro = seq - 1; // DNA<->Protein numbering
  int ignored = 0; // number of $ characters in matched seq
  double sum_log_p = 0;
  const int offset = -cols->front();
  for (motif_columns_type::const_iterator it = cols->begin(); it != cols->end(); ++it) {
    const int pos_pro = startPos - 1 + *it + offset; // DNA indices start at 1
    //		if  (seq_pro>=(int)props.size() || pos_pro >= (int)props[seq_pro].size()) {
    //			fprintf(
    //					stderr,
    //					"ERROR: index: props[%d][%d]\tsize(props)=%d\tsize(props[%d])=%d\tid=%s\t|%s[%d]|=%d\n",
    //					seq_pro, pos_pro, static_cast<int>(props.size()), seq_pro,
    //					static_cast<int>(props[seq_pro].size()),
    //					Global::negSet->entity[seq]->info[0],
    //					"posSet", seq,
    //					Global::posSet->entity[seq]->n);
    //		}
    if (MProGlobal::getS().alphabet[ps[seq_pro][pos_pro]]
        == MProGlobal::getS().alphabet.startStopChar) {
      /* ignore sequence ends (no values for $) */
      ++ignored;
    } else {
      const int pred_index = static_cast<int> (RES * nn->predict(
          props[seq_pro][pos_pro]));
      assert(pred_index >= 0 && pred_index <= RES);
      sum_log_p += log_ptable[pred_index];
    }
  }

  const double log_P_combined = sum_log_p / (static_cast<int> (cols->size())
      - ignored);
  return exp(log_P_combined);
}
