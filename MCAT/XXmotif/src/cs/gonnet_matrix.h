// Copyright 2009, Andreas Biegert

#ifndef CS_GONNET_MATRIX_H_
#define CS_GONNET_MATRIX_H_

#include "amino_acid.h"
#include "substitution_matrix-inl.h"

namespace cs {

// Gonnet substitution matrix
class GonnetMatrix : public SubstitutionMatrix<AminoAcid> {
 public:
  GonnetMatrix();
  virtual ~GonnetMatrix() {}
};

}  // namespace cs

#endif  // CS_GONNET_MATRIX_H_
