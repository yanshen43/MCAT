// Copyright 2009, Andreas Biegert

#include "blosum_matrix.h"

#include "exception.h"

namespace {

const float g_blosum45[] = {
  //A    R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
  0.0181f,
  0.0029f,0.0130f,
  0.0026f,0.0020f,0.0079f,
  0.0027f,0.0020f,0.0031f,0.0132f,
  0.0015f,0.0006f,0.0007f,0.0006f,0.0094f,
  0.0024f,0.0025f,0.0016f,0.0017f,0.0004f,0.0057f,
  0.0038f,0.0031f,0.0022f,0.0047f,0.0007f,0.0032f,0.0131f,
  0.0063f,0.0022f,0.0033f,0.0028f,0.0010f,0.0019f,0.0025f,0.0285f,
  0.0013f,0.0013f,0.0013f,0.0012f,0.0003f,0.0010f,0.0014f,0.0012f,0.0058f,
  0.0036f,0.0016f,0.0015f,0.0014f,0.0008f,0.0013f,0.0017f,0.0019f,0.0007f,0.0124f,
  0.0052f,0.0029f,0.0019f,0.0022f,0.0015f,0.0022f,0.0031f,0.0032f,0.0016f,0.0093f,0.0263f,
  0.0037f,0.0061f,0.0028f,0.0028f,0.0008f,0.0029f,0.0045f,0.0030f,0.0013f,0.0020f,0.0031f,0.0120f,
  0.0016f,0.0010f,0.0007f,0.0006f,0.0004f,0.0008f,0.0009f,0.0011f,0.0006f,0.0023f,0.0041f,0.0011f,0.0026f,
  0.0021f,0.0014f,0.0011f,0.0010f,0.0007f,0.0007f,0.0013f,0.0017f,0.0008f,0.0031f,0.0057f,0.0015f,0.0012f,0.0124f,
  0.0024f,0.0013f,0.0011f,0.0015f,0.0004f,0.0011f,0.0023f,0.0022f,0.0007f,0.0016f,0.0019f,0.0020f,0.0007f,0.0009f,0.0160f,
  0.0062f,0.0026f,0.0031f,0.0028f,0.0011f,0.0024f,0.0032f,0.0049f,0.0013f,0.0023f,0.0032f,0.0033f,0.0010f,0.0017f,0.0020f,0.0104f,
  0.0041f,0.0020f,0.0025f,0.0023f,0.0010f,0.0015f,0.0025f,0.0027f,0.0009f,0.0028f,0.0038f,0.0028f,0.0011f,0.0017f,0.0019f,0.0047f,0.0086f,
  0.0006f,0.0004f,0.0002f,0.0002f,0.0001f,0.0003f,0.0004f,0.0006f,0.0001f,0.0005f,0.0008f,0.0005f,0.0002f,0.0008f,0.0003f,0.0004f,0.0003f,0.0053f,
  0.0017f,0.0014f,0.0009f,0.0011f,0.0004f,0.0010f,0.0012f,0.0014f,0.0012f,0.0019f,0.0031f,0.0015f,0.0009f,0.0034f,0.0007f,0.0015f,0.0013f,0.0008f,0.0066f,
  0.0055f,0.0021f,0.0016f,0.0017f,0.0012f,0.0014f,0.0023f,0.0025f,0.0008f,0.0094f,0.0088f,0.0025f,0.0022f,0.0031f,0.0016f,0.0031f,0.0038f,0.0004f,0.0019f,0.0141f };

const float g_blosum62[] = {
  //A    R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
  0.0215f,
  0.0023f,0.0178f,
  0.0019f,0.0020f,0.0141f,
  0.0022f,0.0016f,0.0037f,0.0213f,
  0.0016f,0.0004f,0.0004f,0.0004f,0.0119f,
  0.0019f,0.0025f,0.0015f,0.0016f,0.0003f,0.0073f,
  0.0030f,0.0027f,0.0022f,0.0049f,0.0004f,0.0035f,0.0161f,
  0.0058f,0.0017f,0.0029f,0.0025f,0.0008f,0.0014f,0.0019f,0.0378f,
  0.0011f,0.0012f,0.0014f,0.0010f,0.0002f,0.0010f,0.0014f,0.0010f,0.0093f,
  0.0032f,0.0012f,0.0010f,0.0012f,0.0011f,0.0009f,0.0012f,0.0014f,0.0006f,0.0184f,
  0.0044f,0.0024f,0.0014f,0.0015f,0.0016f,0.0016f,0.0020f,0.0021f,0.0010f,0.0114f,0.0371f,
  0.0033f,0.0062f,0.0024f,0.0024f,0.0005f,0.0031f,0.0041f,0.0025f,0.0012f,0.0016f,0.0025f,0.0161f,
  0.0013f,0.0008f,0.0005f,0.0005f,0.0004f,0.0007f,0.0007f,0.0007f,0.0004f,0.0025f,0.0049f,0.0009f,0.0040f,
  0.0016f,0.0009f,0.0008f,0.0008f,0.0005f,0.0005f,0.0009f,0.0012f,0.0008f,0.0030f,0.0054f,0.0009f,0.0012f,0.0183f,
  0.0022f,0.0010f,0.0009f,0.0012f,0.0004f,0.0008f,0.0014f,0.0014f,0.0005f,0.0010f,0.0014f,0.0016f,0.0004f,0.0005f,0.0191f,
  0.0063f,0.0023f,0.0031f,0.0028f,0.0010f,0.0019f,0.0030f,0.0038f,0.0011f,0.0017f,0.0024f,0.0031f,0.0009f,0.0012f,0.0017f,0.0126f,
  0.0037f,0.0018f,0.0022f,0.0019f,0.0009f,0.0014f,0.0020f,0.0022f,0.0007f,0.0027f,0.0033f,0.0023f,0.0010f,0.0012f,0.0014f,0.0047f,0.0125f,
  0.0004f,0.0003f,0.0002f,0.0002f,0.0001f,0.0002f,0.0003f,0.0004f,0.0002f,0.0004f,0.0007f,0.0003f,0.0002f,0.0008f,0.0001f,0.0003f,0.0003f,0.0065f,
  0.0013f,0.0009f,0.0007f,0.0006f,0.0003f,0.0007f,0.0009f,0.0008f,0.0015f,0.0014f,0.0022f,0.0010f,0.0006f,0.0042f,0.0005f,0.0010f,0.0009f,0.0009f,0.0102f,
  0.0051f,0.0016f,0.0012f,0.0013f,0.0014f,0.0012f,0.0017f,0.0018f,0.0006f,0.0120f,0.0095f,0.0019f,0.0023f,0.0026f,0.0012f,0.0024f,0.0036f,0.0004f,0.0015f,0.0196f };

const float g_blosum80[] = {
  //A    R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
  0.0252f,
  0.0020f,0.0210f,
  0.0016f,0.0017f,0.0166f,
  0.0018f,0.0013f,0.0037f,0.0262f,
  0.0015f,0.0003f,0.0004f,0.0003f,0.0172f,
  0.0017f,0.0024f,0.0014f,0.0014f,0.0003f,0.0094f,
  0.0028f,0.0023f,0.0019f,0.0048f,0.0003f,0.0035f,0.0208f,
  0.0053f,0.0015f,0.0025f,0.0022f,0.0006f,0.0011f,0.0017f,0.0463f,
  0.0009f,0.0012f,0.0012f,0.0008f,0.0002f,0.0011f,0.0012f,0.0008f,0.0104f,
  0.0027f,0.0010f,0.0007f,0.0008f,0.0011f,0.0007f,0.0010f,0.0009f,0.0004f,0.0220f,
  0.0036f,0.0018f,0.0011f,0.0011f,0.0014f,0.0014f,0.0015f,0.0016f,0.0008f,0.0111f,0.0442f,
  0.0029f,0.0061f,0.0022f,0.0020f,0.0004f,0.0028f,0.0036f,0.0020f,0.0010f,0.0012f,0.0019f,0.0190f,
  0.0011f,0.0006f,0.0004f,0.0003f,0.0004f,0.0007f,0.0006f,0.0005f,0.0003f,0.0025f,0.0052f,0.0007f,0.0053f,
  0.0014f,0.0007f,0.0006f,0.0006f,0.0005f,0.0005f,0.0006f,0.0009f,0.0007f,0.0027f,0.0052f,0.0007f,0.0010f,0.0211f,
  0.0021f,0.0009f,0.0007f,0.0009f,0.0003f,0.0007f,0.0012f,0.0010f,0.0004f,0.0007f,0.0012f,0.0012f,0.0003f,0.0004f,0.0221f,
  0.0064f,0.0020f,0.0029f,0.0024f,0.0010f,0.0017f,0.0026f,0.0034f,0.0010f,0.0015f,0.0021f,0.0026f,0.0007f,0.0010f,0.0014f,0.0167f,
  0.0036f,0.0015f,0.0020f,0.0016f,0.0009f,0.0012f,0.0019f,0.0019f,0.0007f,0.0024f,0.0028f,0.0020f,0.0009f,0.0011f,0.0011f,0.0048f,0.0156f,
  0.0003f,0.0002f,0.0001f,0.0001f,0.0001f,0.0002f,0.0002f,0.0003f,0.0001f,0.0003f,0.0006f,0.0002f,0.0002f,0.0007f,0.0001f,0.0002f,0.0002f,0.0087f,
  0.0011f,0.0007f,0.0006f,0.0005f,0.0003f,0.0005f,0.0006f,0.0006f,0.0016f,0.0013f,0.0020f,0.0008f,0.0005f,0.0046f,0.0003f,0.0009f,0.0008f,0.0010f,0.0148f,
  0.0046f,0.0013f,0.0009f,0.0010f,0.0013f,0.0010f,0.0015f,0.0014f,0.0005f,0.0123f,0.0089f,0.0015f,0.0022f,0.0022f,0.0010f,0.0021f,0.0033f,0.0004f,0.0012f,0.0246f };

}  // namnespace


namespace cs {

BlosumMatrix::BlosumMatrix(Type matrix) {
  switch (matrix) {
    case BLOSUM45:
      Init(g_blosum45);
      break;
    case BLOSUM62:
      Init(g_blosum62);
      break;
    case BLOSUM80:
      Init(g_blosum80);
      break;
    default:
      throw Exception("Unsupported BLOSUM matrix!");
  }
}

void BlosumMatrix::Init(const float* blosum_xx) {
  // Read raw BLOSUM data vector
  int n = 0;
  for (int a = 0; a < size_; ++a)
    for (int b = 0; b <= a; ++b, ++n)
      p_[a][b] = blosum_xx[n];

  // Add uppper right matrix part
  for (int a = 0; a < size_-1; ++a)
    for (int b = a+1; b < size_; ++b)
      p_[a][b] = p_[b][a];

  // Let base class init method do the rest.
  init_from_target_frequencies();
}

BlosumMatrix::Type blosum_matrix_type_from_string(const std::string& s) {
  if (s == "BLOSUM45")
    return BlosumMatrix::BLOSUM45;
  else if (s == "BLOSUM62")
    return BlosumMatrix::BLOSUM62;
  else if (s == "BLOSUM80")
    return BlosumMatrix::BLOSUM80;
  else
    throw Exception("Unable to infer BLOSUM matrix from string '%s'!", s.c_str());
}

}  // namespace cs