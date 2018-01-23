/*
 * Copyright (C) 2006-2016  Music Technology Group - Universitat Pompeu Fabra
 *
 * This file is part of Essentia
 *
 * Essentia is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the Free
 * Software Foundation (FSF), either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the Affero GNU General Public License
 * version 3 along with this program.  If not, see http://www.gnu.org/licenses/
 */

#ifndef ESSENTIA_VITERBI_H
#define ESSENTIA_VITERBI_H

#include "algorithm.h"

namespace essentia {
namespace standard {

class Viterbi : public Algorithm {

 private:
  Input<std::vector<std::vector<Real> > > _transitions;
  Input<std::vector<std::vector<Real> > > _emissions;
  Input<std::vector<Real> > _pi;

  Output<std::vector<std::vector<Real> > > _lattProbs;
  Output<std::vector<std::vector<Real> > > _lattPtrs;

 protected:
  bool _useLog;
  int _forcedAlignment;

 public:
  Viterbi() {
    declareInput(_transitions, "transitionMatrix", "the transition matrix from the HMM");
    declareInput(_emissions, "emisionMatrix", "matrix describing de emissions for the states");
    declareInput(_pi, "pi", "vector containing the probabilities for the initial state");

    declareOutput(_lattProbs, "latticeProbabilities", "matrix with the resulting lattice probabilites");
    declareOutput(_lattPtrs, "latticePointers", "matrix with the pointers to compute the best path");
    
  }

  void declareParameters() {
    declareParameter("forcedAlignment", "if 1, restricts the next state search range to the current or next state", "[0,1]", 0);
    declareParameter("log", "if true use log-probabilities", "{true,false}", true);
    //todo: optimization implementation flag
  }

  void configure();
  void compute();

  static const char* name;
  static const char* category;
  static const char* description;

};

} // namespace standard
} // namespace essentia

#include "streamingalgorithmwrapper.h"

namespace essentia {
namespace streaming {

class Viterbi : public StreamingAlgorithmWrapper {

 protected:
  Sink<std::vector<std::vector<Real> > > _transitions;
  Sink<std::vector<std::vector<Real> > > _emissions;
  Sink<std::vector<Real> > _pi;

  Source<std::vector<std::vector<Real> > > _lattProbs;
  Source<std::vector<std::vector<int> > > _lattPtrs;

 public:
  Viterbi() {
    declareAlgorithm("Viterbi");
    declareInput(_transitions, TOKEN, "transitionMatrix");
    declareInput(_emissions, TOKEN, "emisionMatrix");
    declareInput(_pi, TOKEN, "pi");
    declareOutput(_lattProbs, TOKEN, "latticeProbabilities");
    declareOutput(_lattPtrs, TOKEN, "latticePointers");
  }
};

} // namespace streaming
} // namespace essentia

#endif // ESSENTIA_VITERBI_H
