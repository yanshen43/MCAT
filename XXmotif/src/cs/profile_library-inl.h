// Copyright 2009, Andreas Biegert

#ifndef SRC_PROFILE_LIBRARY_INL_H_
#define SRC_PROFILE_LIBRARY_INL_H_

#include "profile_library.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>

#include "context_profile-inl.h"
#include "count_profile-inl.h"
#include "exception.h"
#include "profile-inl.h"
#include "log.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
ProfileLibrary<Alphabet>::ProfileLibrary(int num_profiles, int num_cols)
    : num_profiles_(num_profiles),
      num_cols_(num_cols),
      iterations_(0),
      profiles_(),
      logspace_(false) {
  profiles_.reserve(num_profiles);
}

template<class Alphabet>
ProfileLibrary<Alphabet>::ProfileLibrary(FILE* fin)
    : num_profiles_(0),
      num_cols_(0),
      iterations_(0),
      profiles_(),
      logspace_(false) {
  Read(fin);
}

template<class Alphabet>
ProfileLibrary<Alphabet>::ProfileLibrary(
    int num_profiles,
    int num_cols,
    const ProfileInitializer<Alphabet>& profile_init)
    : num_profiles_(num_profiles),
      num_cols_(num_cols),
      iterations_(0),
      profiles_(),
      logspace_(false) {
  profiles_.reserve(num_profiles);
  profile_init.Init(*this);
}

template<class Alphabet>
inline void ProfileLibrary<Alphabet>::clear() {
  profiles_.clear();
  profiles_.reserve(num_profiles());
}

template<class Alphabet>
inline int ProfileLibrary<Alphabet>::AddProfile(
    const Profile<Alphabet>& profile) {
  if (full())
    throw Exception("Profile library contains already %i profiles!",
                    num_profiles());
  if (profile.num_cols() != num_cols())
    throw Exception("Profile to add as state has %i columns but should have %i!",
                    profile.num_cols(), num_cols());

  shared_ptr< ContextProfile<Alphabet> >
    profile_ptr(new ContextProfile<Alphabet>(profiles_.size(), profile));
  profile_ptr->set_prior(1.0f / num_profiles());

  profiles_.push_back(profile_ptr);
  return profiles_.size() - 1;
}

template<class Alphabet>
inline void ProfileLibrary<Alphabet>::TransformToLogSpace() {
  if (!logspace()) {
    for (profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
      (*pi)->TransformToLogSpace();
    logspace_ = true;
  }
}

template<class Alphabet>
inline void ProfileLibrary<Alphabet>::TransformToLinSpace() {
  if (logspace()) {
    for (profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
      (*pi)->TransformToLinSpace();
    logspace_ = false;
  }
}

template<class Alphabet>
void ProfileLibrary<Alphabet>::Read(FILE* fin) {
//  LOG(DEBUG1) << "Reading profile library from stream ...";

  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Check if stream actually contains a serialized HMM
  while (fgetline(buffer, kBufferSize, fin))
    if (strscn(buffer)) break;
  if (!strstr(buffer, "ContextLibrary"))
    throw Exception("Context library does not start with 'ContextLibrary' "
                    "keyword!");

  // Read number of profiles
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "SIZE")) {
    ptr = buffer;
    num_profiles_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'SIZE' record!");
  }
  // Read number of columns
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "LENG")) {
    ptr = buffer;
    num_cols_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'LENG' record!");
  }
  logspace_ = true;

  // Read context profiles
  profiles_.reserve(num_profiles());
  while (!full() && !feof(fin)) {
    shared_ptr< ContextProfile<Alphabet> > p(new ContextProfile<Alphabet>(fin));
    profiles_.push_back(p);
  }
  if (!full())
    throw Exception("Profile library has %i profiles but should have %i!",
                    profiles_.size(), num_profiles());
//  LOG(ERROR) << *this;
}

template<class Alphabet>
void ProfileLibrary<Alphabet>::Write(FILE* fout) const {
  // Write header
  fputs("ProfileLibrary\n", fout);
  fprintf(fout, "NPROF\t%i\n", num_profiles());
  fprintf(fout, "NCOLS\t%i\n", num_cols());
  fprintf(fout, "ITERS\t%i\n", iterations());
  fprintf(fout, "LOG\t%i\n", logspace() ? 1 : 0);

  // Serialize profiles
  for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
    (*pi)->Write(fout);
}

template<class Alphabet>
void ProfileLibrary<Alphabet>::Print(std::ostream& out) const {
  out << "ProfileLibrary" << std::endl;
  out << "Total number of profiles: " << num_profiles() << std::endl;
  out << "Context profile columns:  " << num_cols() << std::endl;
  out << "Clustering iterations:    " << iterations() << std::endl;

  for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
    out << **pi;
}


}  // namespace cs

#endif  // SRC_PROFILE_LIBRARY_INL_H_
