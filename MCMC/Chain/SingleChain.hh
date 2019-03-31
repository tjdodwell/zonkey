#ifndef ZONKEY_MCMC_SINGLE_CHAIN_HH
#define ZONKEY_MCMC_SINGLE_CHAIN_HH

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


  private:


      std::vector<Link> theChain;

  };

}
}
#endif /* chain_h */
