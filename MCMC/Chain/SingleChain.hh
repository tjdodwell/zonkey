#ifndef ZONKEY_MCMC_SINGLE_CHAIN_HH
#define ZONKEY_MCMC_SINGLE_CHAIN_HH

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include "getEss.hh"

using namespace Eigen;
using namespace std;

#include <vector>

namespace Zonkey {
  namespace MCMC {

  template<class Link>
  class SingleChain{

  public:

      SingleChain(){ }


      void addLink(Link& newLink){  theChain.push_back(newLink);  }

      Link& back(){return theChain.back(); }

      Link& operator[] (const int index){
        return theChain[index];
      }

      void inline resize(int N){  theChain.resize(N); }

      int inline size(){ return theChain.size(); }

      Eigen::VectorXd EffectiveSampleSizes(){
        int numParam = theChain[0].size(); // Number of Parameters
        Eigen::VectorXd ESS(numParam);
        for (int j = 0; j < numParam; j++){
          std::vector<double> vals(theChain.size());
          for (int i = 0; i < theChain.size(); i++){
            vals[i] = theChain[i].getTheta(j);
          }
          ESS(j) = getESS(vals);
        }
        return ESS;
      }

      double getMaxESS(){
        Eigen::VectorXd ESS = this->EffectiveSampleSizes();
        double max = ESS.maxCoeff();
        return max;
      }




  private:

    std::vector<Link> theChain;


  };

}
}
#endif /* chain_h */
