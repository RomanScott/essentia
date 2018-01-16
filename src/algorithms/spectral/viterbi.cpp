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

#include "viterbi.h"
#include "essentiamath.h"

using namespace essentia;
using namespace standard;

const char* Viterbi::name = "Viterbi";
const char* Viterbi::category = "Statistics";
const char* Viterbi::description = DOC("This algorithm computes the viterbi algorithm as in [1].\n" 
"[1] https://en.wikipedia.org/wiki/Viterbi_algorithm\n");


void Viterbi::configure() {
  _useLog = parameter("log").toBool();
}

void Viterbi::compute() {
  const std::vector<std::vector<Real> >& transitions = _transitions.get();
  const std::vector<std::vector<Real> >& emissions = _emissions.get();
  const std::vector<Real>& pi = _pi.get();

  std::vector<std::vector<Real> >& lattProbs = _lattProbs.get();
  std::vector<std::vector<Real> >& lattPtrs = _lattPtrs.get();



  // Sanity checks
  if ( transitions.size() != transitions[0].size() )
    throw EssentiaException("transition matrix should be squared");

  if ( transitions.size() != pi.size() )
    throw EssentiaException("initial probability vector pi should match the size of the transition matrix");

  int nStates = transitions.size();
  int nObs = emissions[0].size();

  // Init lattice and pointer latt
  Real currentPi;

  lattProbs.resize(nStates);
  lattPtrs.resize(nStates);
  
  for(int i=0; i<nStates; i++) {
    _useLog ? currentPi = log(pi[i]) : currentPi = pi[i];
    
    lattProbs[i].resize(nObs);
    lattPtrs[i].resize(nObs);

    _useLog ? lattProbs[i][0] = currentPi + emissions[i][0] : currentPi * emissions[i][0]; // continue in this line
  }

  std::vector<Real> candidates(nStates);
  for(int i=1; i<nObs; i++) {
    for(int j=0; j<nStates; j++) {  // review  for the forced case
      for(int k=0; k<nStates; k++) {
        candidates[k] = _useLog ? lattProbs[k][i-1] + transitions[k][j] + emissions[j][i] : 
          lattProbs[k][i-1] * transitions[k][j] * emissions[j][i];
      }
      lattPtrs[j][i] = (float)argmax(candidates);
      lattProbs[j][i] = candidates[(int)lattPtrs[j][i]];
    }
  }
}
