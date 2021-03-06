// Copyright 2009, Andreas Biegert

#ifndef SRC_PROFILE_INL_H_
#define SRC_PROFILE_INL_H_

#include "profile.h"

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <vector>

#include "exception.h"
#include "log.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
const char* Profile<Alphabet>::kClassID = "Profile";

template<class Alphabet>
Profile<Alphabet>::Profile()
    : logspace_(false) {}

template<class Alphabet>
Profile<Alphabet>::Profile(int num_cols)
    : data_(num_cols, Alphabet::instance().size() + 1, 0.0f),
      logspace_(false) {
  const int any = Alphabet::instance().any();
  for (int i = 0; i < num_cols; ++i)
    data_[i][any] = 1.0f;
}

template<class Alphabet>
Profile<Alphabet>::Profile(FILE* fin)
    : logspace_(false) {
  Read(fin);
}

template<class Alphabet>
Profile<Alphabet>::Profile(const Profile& other,
                                  int index,
                                  int length)
    : data_(length, Alphabet::instance().size() + 1, 0.0f),
      logspace_(other.logspace_) {
  if (index + length > other.num_cols())
    throw Exception("Index=%i and length=%i of sub-profile are out of bounds!",
                    index, length);
  for (int i = 0; i < num_cols(); ++i)
    for (int a = 0; a <= alphabet_size(); ++a)
      data_[i][a] = other[i+index][a];
}

template<class Alphabet>
void Profile<Alphabet>::ReadAll(FILE* fin,
                                       std::vector< shared_ptr<Profile> >* v) {
  while (!feof(fin)) {
    shared_ptr<Profile> p(new Profile(fin));
    v->push_back(p);

    int c = getc(fin);
    if (c == EOF) break;
    ungetc(c, fin);
  }
}

template<class Alphabet>
void Profile<Alphabet>::TransformToLogSpace() {
  if (!logspace_) {
    for (int i = 0; i < num_cols(); ++i)
      for (int a = 0; a <= alphabet_size(); ++a)
        data_[i][a] = log2(data_[i][a]);
    logspace_ = true;
  }
}

template<class Alphabet>
void Profile<Alphabet>::TransformToLinSpace() {
  if (logspace_) {
    for (int i = 0; i < num_cols(); ++i)
      for (int a = 0; a <= alphabet_size(); ++a)
        data_[i][a] = static_cast<float>(pow(2,data_[i][a]));
    logspace_ = false;
  }
}

template<class Alphabet>
void Profile<Alphabet>::Read(FILE* fin) {
//  LOG(DEBUG1) << "Reading profile from stream ...";

  // Check if stream actually contains a serialized profile
  char buffer[kBufferSize];
  while (fgetline(buffer, kBufferSize, fin))
    if (strscn(buffer)) break;
  if (!strstr(buffer, class_id()))
    throw Exception("Bad format: profile does not start with '%s'!",
                    class_id());

  ReadHeader(fin);
  ReadBody(fin);

//  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void Profile<Alphabet>::ReadHeader(FILE* fin) {
  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Read num_cols
  int num_cols = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "LENG")) {
    ptr = buffer;
    num_cols = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'LENG' record!");
  }
  // Read alphabet_size
  int alphabet_size = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "ALPH")) {
    ptr = buffer;
    alphabet_size = strtoi(ptr);
  } else {
    throw Exception("Bad format: profile does not contain 'ALPH' record!");
  }
  if (alphabet_size != Alphabet::instance().size())
    throw Exception("Bad format: profile alphabet_size should be %i but is %i!",
                    Alphabet::instance().size(), alphabet_size);

  Resize(num_cols);
}

template<class Alphabet>
void Profile<Alphabet>::ReadBody(FILE* fin) {
  const int alph_size = alphabet_size();
  char buffer[kBufferSize];
  const char* ptr = buffer;
  int i = 0;

  fgetline(buffer, kBufferSize, fin);  // skip alphabet description line
  while (fgetline(buffer, kBufferSize, fin)
         && buffer[0] != '/' && buffer[1] != '/') {
    ptr = buffer;
    i = strtoi(ptr) - 1;
    for (int a = 0; a < alph_size; ++a) {
      if (logspace_)
        data_[i][a] = static_cast<float>(log2(exp(static_cast<float>(-strtoi_ast(ptr)) / kLogScale)));
      else
        data_[i][a] = static_cast<float>(exp(static_cast<float>(-strtoi_ast(ptr)) / kLogScale));
    }
    data_[i][alph_size] = logspace_ ? 0.0f : 1.0f;
  }
  if (i != num_cols() - 1)
    throw Exception("Bad format: profile has %i columns but should have %i!",
                    i+1, num_cols());
}

template<class Alphabet>
void Profile<Alphabet>::Write(FILE* fout) const {
  fprintf(fout, "%s\n", class_id());
  WriteHeader(fout);
  WriteBody(fout);
}

template<class Alphabet>
void Profile<Alphabet>::WriteHeader(FILE* fout) const {
  fprintf(fout, "NCOLS\t%i\n", num_cols());
  fprintf(fout, "ALPH\t%i\n", alphabet_size());
  fprintf(fout, "LOG\t%i\n", logspace() ? 1 : 0);
}

template<class Alphabet>
void Profile<Alphabet>::WriteBody(FILE* fout) const {
  fputs("PROF\t", fout);
  Alphabet::instance().Write(fout);
  fputc('\n', fout);

  for (int i = 0; i < num_cols(); ++i) {
    fprintf(fout, "%i", i+1);
    for (int a = 0; a < alphabet_size(); ++a) {
      const float log_p = logspace_ ? data_[i][a] : static_cast<float>(log(data_[i][a])/log(2));
      if (log_p == -INFINITY)
        fputs("\t*", fout);
      else
        fprintf(fout, "\t%i", -iround(log_p * kLogScale));
    }
    fputc('\n', fout);
  }
  fputs("//\n", fout);
}

template<class Alphabet>
void Profile<Alphabet>::Print(std::ostream& out) const {
  out << "\t" << Alphabet::instance() << std::endl;

  for (int i = 0; i < num_cols(); ++i) {
    out << i+1;
    for (int a = 0; a < alphabet_size(); ++a)
      out << strprintf("\t%6.4f",
                       logspace_ ? pow(2,data_[i][a]) : data_[i][a]);
    out << std::endl;
  }
}

template<class Alphabet>
void Profile<Alphabet>::Resize(int num_cols) {
  if (num_cols == 0)
    throw Exception("Cannot resize profile to num_cols=%i!", num_cols);
  data_.resize(num_cols, Alphabet::instance().size() + 1, 0.0f);
}

template<class Alphabet>
inline void Reset(Profile<Alphabet>* p) {
  Profile<Alphabet>& profile = *p;
  const int num_cols = profile.num_cols();
  const int alphabet_size = profile.alphabet_size();

  for (int i = 0; i < num_cols; ++i)
    for (int a = 0; a < alphabet_size; ++a)
      profile[i][a] = 0.0f;
}

template<class Alphabet>
bool Normalize(Profile<Alphabet>* p, float value) {
  Profile<Alphabet>& profile = *p;
  const bool logspace = profile.logspace();
  if (logspace) profile.TransformToLinSpace();

  const int num_cols       = profile.num_cols();
  const int alphabet_size  = profile.alphabet_size();
  bool rv = true;

  for (int i = 0; i < num_cols; ++i) {
    float sum = 0.0f;
    for (int a = 0; a < alphabet_size; ++a) sum += profile[i][a];
    if (sum != 0.0f) {
      float fac = value / sum;
      for (int a = 0; a < alphabet_size; ++a) profile[i][a] *= fac;
    } else {
      rv = false;  // couldn't normalize at least one column
    }
  }

  if (logspace) profile.TransformToLogSpace();
  return rv;
}

}  // namespace cs

#endif  // SRC_PROFILE_INL_H_
