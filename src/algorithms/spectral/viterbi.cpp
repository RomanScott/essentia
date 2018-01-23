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
  _forcedAlignment = parameter("forcedAlignment").toInt();
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

  int i, j, k, idx, nStates = transitions.size(), nObs = emissions[0].size();
  Real currentPi;

  lattProbs.resize(nStates);
  lattPtrs.resize(nStates);
  
  // Init probabilities and pointer lattices
  for (i=0; i<nStates; i++) {
    _useLog ? currentPi = log(pi[i]) : currentPi = pi[i];
    
    lattProbs[i].resize(nObs);
    lattPtrs[i].resize(nObs);

    _useLog ? lattProbs[i][0] = currentPi + emissions[i][0] : currentPi * emissions[i][0];
  }

  std::vector<Real> candidates(nStates);
  int firstCandidate = 0;
  int lastCandidate = nStates;

  // Viterbi loop
  for (i=1; i<nObs; i++) {
    for (j=0; j<nStates; j++) {
      
      // Restricted search for the forced alignmet case 
      if (_forcedAlignment == 1) {
        firstCandidate = std::max(j-2, 0);
        lastCandidate = j+1;
        candidates.resize(lastCandidate - firstCandidate);
      }

      for (k=firstCandidate, idx=0; k<lastCandidate; k++, idx++) {
        _useLog ? candidates[idx] = lattProbs[k][i-1] + transitions[k][j] : 
          candidates[idx] = lattProbs[k][i-1] * transitions[k][j];
      }

      // Updating the lattices
      lattPtrs[j][i] = (float)(argmax(candidates) + firstCandidate);
      lattProbs[j][i] = candidates[(int)lattPtrs[j][i] - firstCandidate];
      _useLog ? lattProbs[j][i] += emissions[j][i] : 
        lattProbs[j][i] *= emissions[j][i];
    }
  }
}
