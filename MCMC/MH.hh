#ifndef ZONKEY_MCMC_MH_HH
#define ZONKEY_MCMC_MH_HH

using namespace Eigen;
using namespace std;

#include <vector>

namespace Zonkey {
  namespace MCMC {

  template<typename Link, typename Chain, typename PROPOSAL, typename ForwardModel>
  class MetropolisHastings{

  public:

      MetropolisHastings(ForwardModel & F_, PROPOSAL& proposal_, Chain & markovChain_):
        F(F_),
        proposal(proposal_),
        markovChain(markovChain_)
        {  }


      void inline burnin(int est_ACT, int factor = 10){
        run(est_ACT);
        int currentEss = markovChain.getMaxESS();
        double moreSamples = currentEss - factor;
        if (moreSamples > 0){ run(moreSamples); }
      }


      void inline run(int numSamples){
          if(markovChain.size() < 1){ // If this is the first sample
            markovChain.push_back (F.samplePrior());
            F.apply(markovChain.back());
            numSamples -= 1;
          }
          for (int i = 0; i < numSamples; i++){
            Link theta_p = proposal.apply(markovChain.back()); // Make a proposal
            F.apply(theta_p); // Apply forward Model


            // Accept / Reject Step
            bool accept = proposal.acceptReject(markovChain.back(),  theta_p)

            if(accept){
              theta_p.setAccepted(1);
              markovChain.push_back(theta_p);
            }
            else {
              markovChain.push_back(markovChain.back());
            }
          }
      }

      int inline size(){
        return markovChain.size();
      }



  private:

      ForwardModel F;
      Chain markovChain;
      PROPOSAL proposal;

  };

}
}
#endif /* chain_h */
