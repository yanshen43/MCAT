#ifndef _NNODE_H
#define _NNODE_H

/**
 * Node of a neural network.
 */

#include <vector>
#include <cstddef>

class NNode {

  private:

    /* the bias of the node */
    double bias;

    /* weights of each input connection */
    std::vector<double> weights;

  public:
    NNode(const double bias, const std::vector<double> &weights) :
      bias(bias), weights(weights) {
    }
    ;

    /* Does the actual prediction. The vector has to have
     * the same size as the node's weight vector.
     */
    double predict(const std::vector<double> &input) const;

    size_t getNumberOfInputs() const {
      return weights.size();
    }

    double getBias() const {
      return bias;
    }

    const std::vector<double>& getWeights() const {
      return weights;
    }
};

#endif
