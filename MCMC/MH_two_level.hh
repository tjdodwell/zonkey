#ifndef ZONKEY_MCMC_MH_TWO_LEVEL_HH
#define ZONKEY_MCMC_MH_TWO_LEVEL_HH

using namespace Eigen;
using namespace std;


#include <random>
#include <cassert>
#include <cmath>

#include <vector>

namespace Zonkey {
  namespace MCMC {

  template<typename Link, typename Chain, typename PROPOSAL, typename ForwardModel>
  class MetropolisHastings_TwoLevel{

  public:

    typedef typename Zonkey::MCMC::MetropolisHastings<Link,Chain,PROPOSAL,ForwardModel> SUBCHAIN_MH;

      MetropolisHastings_TwoLevel(ForwardModel & F_, PROPOSAL& proposal_, Chain & fineChain_, Chain & coarseChain_, int subSamplingRate_ = 10, int level_ = 1):
        F(F_),
        proposal(proposal_),
        coarseChain(coarseChain_),
        fineChain(fineChain_),
        maxLevel(level_),
        subSamplingRate(subSamplingRate_)
        {  

          std::cout << " *** Initialising Two Level MetropolisHastings . . . ";
         // initialise fineChain & courseChain
         Eigen::VectorXd thetaC = F.samplePrior(maxLevel - 1);
         Link coarseStart(thetaC);
         F.apply(coarseStart, maxLevel - 1);
         coarseChain.addLink(coarseStart,1);
         Eigen::VectorXd thetaF = F.samplePrior(maxLevel);
         Link fineStart(thetaF);
         F.apply(fineStart, maxLevel);
         fineChain.addLink(fineStart,1);
         std::cout << "All done!" << std::endl;
        }


      void inline burnin(int est_ACT, int level = 0, int factor = 10){
        run(est_ACT,level,"Initial Burnin . . . ");
        int currentEss = fineChain.getMaxESS() + 1; // Rounding up
        double moreSamples = factor * currentEss - est_ACT;
        if (moreSamples > 0){ run(moreSamples,level, "More burnin samples . . . "); }
        burninSamples = fineChain.size();
        fineChain.setBurninLength(burninSamples);
      }

      void inline setStart(Eigen::VectorXd & xi, int level = 0){
        Link firstPoint(xi);
        F.apply(firstPoint, level);
        if (level == 0){ coarseChain.addLink(firstPoint,1); }
        else{ fineChain.addLink(firstPoint,1);}
      }

      void inline setStart(Link & firstPoint, int level = 0){
        if (level == 0){ coarseChain.addLink(firstPoint,1); }
        else{ fineChain.addLink(firstPoint,1);}
      }

      void inline setBurnin(int burninLength_){ burninLength = burninLength_; }

      void inline run(int numSamples, int level = 0, string printout = "Computing Samples ...", int verb = 1){

            int x = 0;
            if (verb > 0){
              cout << printout << numSamples << " Samples" << endl;
            }

          bool best_seen = false;

          for (int i = 0; i < numSamples; i++){

            if(verb > 1){

              std::cout << "# Sample " << i << std::endl;

            }

            
            // **** Subchain run for coarse proposal

            Chain subChain;
            Link currentF = fineChain.back();
            Link currentC = coarseChain.back(); // Last Element in the Chain

            if(verb > 1){

              std::cout << "Fine Chain - First Mode = " << currentF.getTheta(0) << std::endl;
              std::cout << "Fine Chain - Log Phi = " << currentF.getlogPhi() << std::endl;
               std::cout << "Coarse Chain - First Mode = " << currentC.getTheta(0) << std::endl;
               std::cout << "Coarse Chain - Log Phi = " << currentC.getlogPhi() << std::endl;
              std::cout << " " << std::endl;
            }

            
            SUBCHAIN_MH subMCMC(F,proposal,subChain);

            subMCMC.setStart(currentC);

            subMCMC.run(subSamplingRate,maxLevel - 1,"",0); // We run Markov Chain on level l-1
            auto thisSubChain = subMCMC.getChain(); // Extract the chain
            Link PropC = thisSubChain.back(); // this is the proposal on the coarse modes.
            
            coarseChain.addLink(PropC,1); // Note that coarse proposal is always added to chain

            // More modes? (Not for now)

            Link PropF = PropC;

            F.apply(PropF,maxLevel);  // Apply forward Model on the fine


            if(verb > 1){
              std::cout << "Fine Chain Proposal - First Mode = " << PropF.getTheta(0) << std::endl;
              std::cout << "Fine Chain - Log Phi = " << PropF.getlogPhi() << std::endl;
               std::cout << "Coarse Chain - First Mode = " << PropC.getTheta(0) << std::endl;
               std::cout << "Coarse Chain - Log Phi = " << PropC.getlogPhi() << std::endl;
              std::cout << " " << std::endl;
            }

            // 

            /*if (theta_p.getlogPhi() > bestObserved.getlogPhi()){
              bestObserved = theta_p;
              best_seen = true;
            }*/


            // Accept / Reject Step on the fine
            bool accept = (*this).acceptReject(currentF,currentC,PropF,PropC);

            if(verb > 1){ std::cout << "Proposal Accepted = " << accept << std::endl; }


            // If Accept

            if(accept){     // Accept
              fineChain.addLink(PropF,1);
            }
            else{           // Reject
              fineChain.addLink(currentF,0);
            }

            x++;

            if(verb > 0){
              int num = 100 * x / numSamples;
              cout << "\r" << setw(-20) << printProg(num) << " " << num << "% completed." << flush;
            }

          }
          if(verb > 0){
            cout << " " << endl;
          }

         /* if(best_seen){
            fineChain.setbestObserved(bestObserved);
            best_seen = false;
          }*/
      }

      int inline size(bool isCoarse = false){
        if(isCoarse){ return coarseChain.size();  }
        else{ return fineChain.size();  }
      }

      Chain getChain(bool isCoarse = false){
        if(isCoarse){ return coarseChain; }
        else{ return fineChain; }
      }

      Eigen::VectorXd getY(int i){
        Eigen::VectorXd Qc = coarseChain[i].getQ();
        Eigen::VectorXd Qf = fineChain[i].getQ();
        return Qf - Qc;
      }

      std::vector<double> getY_All(int j){

        std::vector<double> All(fineChain.size() - burninLength);

        for (int i = 0; i < fineChain.size() - burninLength; i++){
          auto Y = (*this).getY(i + burninLength);
          All[i] = Y(j);
        }
        
        return All;
      }

      Eigen::VectorXd getY_mean_variance(int j){

        double meanY = 0.0;
        double meanY2 = 0.0;

        int N = fineChain.size() - burninLength;

        for (int i = 0; i < fineChain.size() - burninLength; i++){
          auto Y = (*this).getY(i + burninLength);
          meanY += Y(j);
          meanY2 += Y(j) * Y(j);
        }
        meanY /= N;
        meanY2 /= N;

        Eigen::VectorXd res(2);

        res(0) = meanY;

        res(1) = meanY2 - meanY * meanY;
        
        return res;
      }

      bool acceptReject(Link& u, Link& uC, Link& v, Link& vC, bool verb = false){
          bool accept = false;
          std::random_device rd;
          std::mt19937 gen(rd());
          std::uniform_real_distribution<> dis(0.0, 1.0);
          double logUniform = std::log(dis(gen));
          double logalpha = v.getlogPhi() + uC.getlogPhi() - u.getlogPhi() - vC.getlogPhi();
          if (logUniform < logalpha){ accept = true; }
          return accept;
        }


      int inline getBurnin(){ return burninLength; }


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
      Chain fineChain, coarseChain;
      PROPOSAL proposal;
      int burninSamples, maxLevel, subSamplingRate;
      int burninLength;

      Link bestObserved;

  };

}
}
#endif // ZONKEY_MCMC_MH_HH
