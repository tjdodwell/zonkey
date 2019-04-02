#ifndef ZONKEY_MCMC_SINGLE_CHAIN_HH
#define ZONKEY_MCMC_SINGLE_CHAIN_HH

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using namespace Eigen;
using namespace std;

namespace Zonkey {
  namespace MCMC {

  template<class Link>
  class SingleChain{

  public:

      SingleChain(){ }


      // void addLink(Link& newLink){
      //
      //   theChain.push_back(newLink);
      //
      // }
      //
      // Link& operator[] (const int index){
      //   return theChain[index];
      // }
      //
      //
      // void inline resize(int N){  theChain.resize(N); }
      //
      // int inline size(){ return theChain.size(); }
      //
      // Eigen::VectorXd EffectiveSampleSizes(){
      //   int numParam = theChain[0].size(); // Number of Parameters
      //   Eigen::VectorXd ESS(numParam);
      //   for (int j = 0; j < numParam; j++){
      //     std::vector<double> vals(theChain.size());
      //     for (int i = 0; i < theChain.size(); i++){
      //       vals[i] = theChain[i].getTheta(j);
      //     }
      //     ESS(j) = getESS(vals);
      //   }
      //
      // double getMaxESS(){
      //   Eigen::VectorXd ESS = EffectiveSampleSizes();
      //   return EffectiveSampleSizes.maxCoeff();
      // }




  private:


      std::vector<Link> theChain;

  };

}
}
#endif /* chain_h */
