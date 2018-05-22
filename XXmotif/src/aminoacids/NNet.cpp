#include "NNet.h"
#include "NNode.h"
#include "../aminoacids/utils.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

NNet::NNet(const std::string &filename) {

  // read network file in memory
  std::ifstream instr(filename.c_str());
  if (!instr.is_open()) {
    std::cerr << "Can't open " << filename << std::endl;
    exit(1);
  }
  std::vector<std::string> lines;
  while (!instr.eof()) {
    std::string line;
    getline(instr, line);
    lines.push_back(line);
  }
  instr.close();

  /* first, determine number of input and output nodes in network */
  int numI = std::numeric_limits<int>::min();
  int numH = std::numeric_limits<int>::min();
  for (std::vector<std::string>::const_iterator s = lines.begin(); s
      != lines.end(); ++s) {
    /* find lines matching "i(\d+)->h(\d+)"
     * and capture $1 in numI and $2 in numH
     */
    if (s->find("\"i") != std::string::npos) {
      const size_t startI = 2;
      const size_t endI = s->find("->h") - 1;
      const size_t startH = endI + 4;
      const size_t endH = s->rfind("\"") - 1;
      const std::string Istr(s->substr(startI, endI - startI + 1));
      const std::string Hstr(s->substr(startH, endH - startH + 1));
      numI = std::max(numI, atoi(Istr.c_str()));
      numH = std::max(numH, atoi(Hstr.c_str()));
    }
  }

  /* now, generate numH hidden nodes with numI inputs each */
  int lineNo = 1; // skip first line, contains only unneeded info
  for (int n = 0; n < numH; ++n) {
    std::string dummy;
    double bias;
    std::stringstream l(lines[lineNo++]);
    l >> dummy >> bias;
    std::vector<double> weights(numI);
    for (size_t w = 0; w < weights.size(); ++w) {
      std::stringstream l(lines[lineNo++]);
      l >> dummy;
      l >> weights[w];
    }
    addHiddenNode(NNode(bias, weights));
  }

  /* output node follows next in the file */
  std::stringstream l(lines[lineNo++]);
  std::string dummy;
  double bias;
  l >> dummy >> bias;
  std::vector<double> weights(numH);
  for (size_t w = 0; w < weights.size(); ++w) {
    std::stringstream l(lines[lineNo++]);
    l >> dummy;
    l >> weights[w];
  }
  addOutputNode(NNode(bias, weights));
}

void NNet::addHiddenNode(const NNode &node) {
  if (outputNodes.size() > 0) {
    std::cerr
        << "ERROR: you can't add further hidden nodes after adding output nodes"
        << std::endl;
    exit(1);
  }
  if (hiddenNodes.size() > 0) {
    if (node.getNumberOfInputs() != hiddenNodes[0].getNumberOfInputs()) {
      std::cerr << "ERROR: can't add node with " << node.getNumberOfInputs()
          << " inputs (other nodes have " << hiddenNodes[0].getNumberOfInputs()
          << " inputs)" << std::endl;
      exit(1);
    }
  }
  hiddenNodes.push_back(node);
}

void NNet::addOutputNode(const NNode &node) {
  if (hiddenNodes.size() == 0) {
    std::cerr << "ERROR: add hidden nodes before output nodes" << std::endl;
    exit(1);
  } else if (node.getNumberOfInputs() != hiddenNodes.size()) {
    std::cerr << "ERROR: output node has " << node.getNumberOfInputs()
        << " inputs, must be " << hiddenNodes.size() << std::endl;
    exit(1);
  } else if (outputNodes.size() != 0) {
    std::cerr << "ERROR: currently, only a single output node is allowed"
        << std::endl;
    exit(1);
  }
  outputNodes.push_back(node);
}

double NNet::predict(const std::vector<double> &input) const {
  if (input.size() != hiddenNodes[0].getNumberOfInputs()) {
    std::cerr << "ERROR: wrong number of inputs" << std::endl;
    exit(1);
  }
  /* collect predictions of all hidden nodes */
  std::vector<double> outH(hiddenNodes.size());
  for (size_t i = 0; i < hiddenNodes.size(); ++i) {
    outH[i] = hiddenNodes[i].predict(input);
  }
  /* and feed them into the output node */
  return outputNodes[0].predict(outH);
}

std::ostream& operator<<(std::ostream &os, const NNode &node) {
  os << "(Bias=" << node.getBias() << ", Weights=" << node.getWeights() << ")";
  return os;
}

std::ostream& operator<<(std::ostream &os, const NNet& net) {
  os << "FFNN with " << net.hiddenNodes.size() << " hidden and "
      << net.outputNodes.size() << " output nodes";
  return os;
}
