#ifndef ZONKEY_MCMC_MH_HH
#define ZONKEY_MCMC_MH_HH

using namespace Eigen;
using namespace std;

namespace Zonkey {
  namespace MCMC {

  template<typename Link, typename Chain, typename PROPOSAL, typename FowardModel>
  class MetropolisHastings{

  public:

      MetropolisHastings(ForwardModel & F, Chain & markovChain_):
        F(F_),
        markovChain(markovChain_)
        {  }


      void inline burnin(int est_ACT, int factor = 10){
        run(est_ACT);
        int currentEss = markovChain.getEss();
        double moreSamples = currentEss - factor;
        if (moreSamples > 0){ run(moreSamples); }
      }


      void inline run(int numSamples){
          if(markovChain.size() < 1){ // If this is the first sample
            Link markovChains.push_back(F.prior());
            F.apply(markovChain.back());
            numSamples -= 1;
          }
          for (int i = 0; i < numSamples; i++){
            Link theta_p = proposal.apply(markovChain.back()); // Make a proposal
            F.apply(theta_p); // Apply forward Model
            markovChain.push_back(proposal.acceptReject(markovChain.back(),  theta_p)); // Accept / Reject Step
          }
      }

      int inline size(){
        return markovChain.size();
      }



  private:

      ForwardModel F;
      Chain markovChain;

  };

}
}
#endif /* chain_h */
