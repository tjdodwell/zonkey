#ifndef ZONKEY_MCMC_SINGLE_CHAIN_HH
#define ZONKEY_MCMC_SINGLE_CHAIN_HH


#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <fstream>

#include <unsupported/Eigen/FFT>

using namespace Eigen;
using namespace std;

namespace Zonkey {
  namespace MCMC {

  template<class Link>
  class SingleChain{

  public:

      SingleChain(){  }


      void addLink(Link& newLink){

        theChain.push_back(newLink);

      }

      Link& operator[] (const int index){
        return theChain[index];
      }


      void inline resize(int N){  theChain.resize(N); }

      int inline size(){ return theChain.size(); }


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

}


  private:


      std::vector<Link> theChain;

  };

}
}
#endif /* chain_h */
