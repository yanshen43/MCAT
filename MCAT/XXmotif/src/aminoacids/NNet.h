/**
 * This class provides a simple FF neural network with a single hidden
 * layer and a single output node. It reads the output generated
 * by R after training with the nnet package.
 *
 * TODO cleanup if we decide to use it in the final version
 */

#ifndef _NNET_H
#define _NNET_H

#include <ostream>
#include <string>
#include <vector>

#include "NNode.h"

class NNet {

    friend std::ostream& operator<<(std::ostream&, const NNet&);

  public:
    /* construct NN from R output file */
    NNet(const std::string &filename);

    /* Do the actual prediction for a vector of input
     * values. It has to be of the same size as
     * the number of hidden nodes.
     */
    double predict(const std::vector<double> &input) const;

  private:
    std::vector<NNode> hiddenNodes;

    /* for future extensions, outputNodes is
     * a vector, however only one output node is supported
     * at the moment
     */
    std::vector<NNode> outputNodes;

    /* adds given node as a hidden node,
     * only used during initial construction
     */
    void addHiddenNode(const NNode &node);

    /* adds given node as the output node,
     * only used during initial construction
     */
    void addOutputNode(const NNode &node);

};

std::ostream& operator<<(std::ostream &os, const NNode &node);
std::ostream& operator<<(std::ostream &os, const NNet& net);

#endif
