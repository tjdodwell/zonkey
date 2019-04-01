#ifndef ZONKEY_MCMC_GET_ESS_HH
#define ZONKEY_MCMC_GET_ESS_HH


#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <fstream>

#include <unsupported/Eigen/FFT>

namespace Zonkey {
  namespace MCMC {

double getESS(std::vector<double>& myChain){

  int numSamples = myChain.size(); // Number of current samples

  assert(numSamples > 0); // Must have at least one sample

  if (numSamples < 2) // Just return 0.0;
    return 0.0;

  double chainMean = std::accumulate(myChain.begin(), myChain.end(),0.0) / numSamples;

  Eigen::FFT<double> fft;

  int tmax    = floor(numSamples / 2);
  double Stau = 1.5;


  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> freqVec;

  Eigen::Matrix<std::complex<double>, Eigen::Dynamic,1> timeVec = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>::Zero(numSamples + tmax);

  for (int i = 0; i < numSamples; i++) {
    timeVec(i) = std::complex<double>(myChain[i] - chainMean, 0.0);
  }


  fft.fwd(freqVec, timeVec);

  // compute out1*conj(out2) and store in out1 (x+yi)*(x-yi)
  double real;

  for (int i = 0; i < numSamples + tmax; i++) {
    real       = freqVec(i).real() * freqVec(i).real() + freqVec(i).imag() * freqVec(i).imag();
    freqVec(i) = std::complex<double>(real, 0.0);
  }


  // now compute the inverse fft to get the autocorrelation (stored in timeVec)
  fft.inv(timeVec, freqVec);

  for (int i = 0; i < tmax + 1; ++i) {
    timeVec(i) = std::complex<double>(timeVec(i).real() / double(numSamples - i), 0.0);
  }

// the following loop uses ideas from "Monte Carlo errors with less errors." by Ulli Wolff to figure out how far we
// need to integrate
//int MaxLag = 0;
double Gint = 0;
int    Wopt = 0;
for (int i = 1; i < tmax + 1; ++i) {
  Gint += timeVec(i).real() / timeVec(0).real(); //

  double tauW;
  if (Gint <= 0) {
    tauW = 1.0e-15;
  } else {
    tauW = Stau / log((Gint + 1) / Gint);
  }
  double gW = exp(-double(i) / tauW) - tauW / sqrt(double(i) * numSamples);

  if (gW < 0) {
    Wopt = i;
    tmax = min(tmax, 2 * i);
    break;
  }
}

// correct for bias
double CFbbopt = timeVec(0).real();
for (int i = 0; i < Wopt + 1; ++i) {
  CFbbopt += 2 * timeVec(i + 1).real();
}

CFbbopt = CFbbopt / numSamples;

for (int i = 0; i < tmax + 1; ++i) {
  timeVec(i) += std::complex<double>(CFbbopt, 0.0);
}

// compute the normalized autocorrelation
double scale = timeVec(0).real();
for (int i = 0; i < Wopt; ++i) {
  timeVec(i) = std::complex<double>(timeVec(i).real() / scale, 0.0);
}

double tauint = 0;
for (int i = 0; i < Wopt; ++i) {
  tauint += timeVec(i).real();
}

tauint -= 0.5;

// return the effective sample size
double frac = 1.0 / (2.0 * tauint);
frac = fmin(1.0, frac);
return numSamples * frac;

}

}

}
