#include <cassert>
#include <cmath>

#include "NNode.h"

double NNode::predict(const std::vector<double> &input) const {
  assert(input.size() == weights.size());
  double in = 0;
  for (size_t i = 0; i < weights.size(); ++i) {
    in += weights[i] * input[i];
  }
  return 1.0 / (1 + exp(-(in + bias)));
}
