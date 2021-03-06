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
        {  

          bestObserved.setlogPhi(-10.0e6); // This will be on the fine chain only since that is all we

        }


      void inline burnin(int est_ACT, int level = 0, int factor = 10){
        run(est_ACT,level,"Initial Burnin . . . ");
        int currentEss = markovChain.getMaxESS() + 1; // Rounding up
        double moreSamples = factor * currentEss - est_ACT;
        if (moreSamples > 0){ run(moreSamples,level, "More burnin samples . . . "); }
        burninSamples = markovChain.size();
        markovChain.setBurninLength(burninSamples);
      }

      void inline setStart(Eigen::VectorXd & xi, int level = 0){
        Link firstPoint(xi);
        F.apply(firstPoint, level);
        markovChain.addLink(firstPoint,1);
      }

      void inline setStart(Link & firstPoint, int level = 0){
        markovChain.addLink(firstPoint);
      }


      void inline run(int numSamples, int level = 0, string printout = "Computing Samples ...", bool verb = true){

          if(markovChain.size() < 1){ // If this is the first sample
            Eigen::VectorXd theta_fP = F.samplePrior(level);
            Link firstPoint(theta_fP);
            F.apply(firstPoint,level);
            markovChain.addLink(firstPoint,1);
            numSamples -= 1;
          }

            int x = 0;
            if (verb){
              cout << printout << numSamples << "Samples" << endl;
            }

          bool best_seen = false;

          for (int i = 0; i < numSamples; i++){

            Link lastLink = markovChain.back();

            Link theta_p = proposal.apply(lastLink); // Make a proposal

            F.apply(theta_p,level);  // Apply forward Model

            // 

            if (theta_p.getlogPhi() > bestObserved.getlogPhi()){
              bestObserved = theta_p;
              best_seen = true;
            }


            // Accept / Reject Step
            bool accept = proposal.acceptReject(lastLink,  theta_p);

            if(accept){ 

              markovChain.addLink(theta_p,1);

            }
            else {

              markovChain.addLink(lastLink,0);  

            }
            

            x++;

            if(verb){
              int num = 100 * x / numSamples;
              cout << "\r" << setw(-20) << printProg(num) << " " << num << "% completed." << flush;
            }

          }
          if(verb){
            cout << " " << endl;
          }

          if(best_seen){
            markovChain.setbestObserved(bestObserved);
            best_seen = false;
          }
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
      int burninSamples;

      Link bestObserved;

  };

}
}
#endif // ZONKEY_MCMC_MH_HH
