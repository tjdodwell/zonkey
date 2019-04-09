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
        run(est_ACT,"Initial Burnin . . . ");
        int currentEss = markovChain.getMaxESS() + 1; // Rounding up
        double moreSamples = est_ACT - factor * currentEss;
        if (moreSamples > 0){ run(moreSamples,"More burnin samples . . . "); }
      }


      void inline run(int numSamples, string printout = "Computing Samples ..."){

          if(markovChain.size() < 1){ // If this is the first sample
            Eigen::VectorXd theta_fP = F.samplePrior();
            Link firstPoint(theta_fP);
            markovChain.addLink(firstPoint);
            F.apply(markovChain.back());
            numSamples -= 1;
          }

            int x = 0;
            cout << printout << numSamples << "Samples" << endl;


          for (int i = 0; i < numSamples; i++){

            Link theta_p = proposal.apply(markovChain.back()); // Make a proposal

            F.apply(theta_p);  // Apply forward Model

            // Accept / Reject Step
            bool accept = proposal.acceptReject(markovChain.back(),  theta_p);

            if(accept){ theta_p.setAccepted(1); markovChain.addLink(theta_p);

            std::cout << "Do I accept?" << std::endl; }
            else {  markovChain.addLink(markovChain.back());  }
            x++;

            int num = 100 * x / numSamples;
            cout << "\r" << setw(-20) << printProg(num) << " " << num << "% completed." << flush;

          }

          cout << " " << endl;
      }

      int inline size(){
        return markovChain.size();
      }

      Chain getChain(){
        return markovChain;
      }



  private:


      string printProg(int x){
          string s;
          s="[";
          for (int i=1;i<=(100/2);i++){
              if (i<=(x/2) || x==100)
                  s+="=";
              else if (i==(x/2))
                  s+=">";
              else
                  s+=" ";
          }

          s+="]";
          return s;
      }

      ForwardModel F;
      Chain markovChain;
      PROPOSAL proposal;

  };

}
}
#endif // ZONKEY_MCMC_MH_HH
