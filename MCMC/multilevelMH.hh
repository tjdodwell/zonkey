


#include <iostream>
#include <cstdio>
#include <ctime>


namespace Zonkey {
  namespace MCMC {


template<class MultiChain, class ForwardModel, class CoarseProposal>
class MultiLevelMetropolisHastings{

  public:

    MultiLevelMetropolisHastings(int L_, ForwardModel & F_, CoarseProposal & coarseProposal_):
      L(L_),
      F(F_),
      coarseProposal(coarseProposal_){

        hierarchy.resize(L + 1);
        // Setup a hierarchy of Multilevel chains
        for (int i = 0; i < L + 1; i++){
          hierarchy[i].setLevel(i);
        }



      }

    //MultiLevelMetropolisHastings(int L_, ForwardModel& F_, CoarseProposal& coarseProposal_):
      //F(F_),
      //coarseProposal(coarseProposal_){

        /*L = L_;

        hierarchy.resize(L + 1);
        // Setup a hierarchy of Multilevel chains
        for (int i = 0; i < L + 1; i++){
          hierarchy[i].setLevel(i);
        }
        isBurnt_in.resize(L + 1);
        // Setup initial Samples fro each level
        NStar.resize(L + 1);
        NStar[0] = 1000;
        for (int i = 1; i < L + 1; i++){
          NStar[i] = std::max(NStar[i-1]/2,10);
        }

        //

        payOff.resize(L + 1);*/

    //}

    private:

      int L;
      ForwardModel F;
      CoarseProposal coarseProposal;


      std::vector<MultiChain> hierarchy;

};

}
}
