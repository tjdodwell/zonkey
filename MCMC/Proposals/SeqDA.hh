#ifndef ZONKEY_MCMC_DELAYED_ACCEPTANCE_HH
#define ZONKEY_MCMC_DELAYED_ACCEPTANCE_HH

#include "../MH.hh"

namespace Zonkey {
  namespace MCMC {

    template<typename LINK, typename CHAIN, typename PROPOSAL, typename ForwardModel>
    class SeqDA{

      public:

        SeqDA(PROPOSAL& myProposal_, ForwardModel& F_, int level_ = 1, int subChain_length_ = 10, bool adaptive_ = false) :
          subChain_length(subChain_length_),
          myProposal(myProposal_),
          F(F_),
          level(level_),
          adaptive(adaptive_){

            count = 0;



        }

        LINK apply(LINK& currentState){

          Eigen::VectorXd xi = currentState.getTheta(); // Obtain current point

          CHAIN markovChain; // Setup a sub chain

          MetropolisHastings<LINK,CHAIN,PROPOSAL,ForwardModel> myMCMC(F,myProposal,markovChain); // Setup up subChain

          myMCMC.setStart(xi,level-1); // Initate start point

          myMCMC.run(subChain_length,level-1,"",false);

          auto theChain = myMCMC.getChain(); // extract the chain

          logCoarse = theChain[0].getlogPhi(); // Obtain Log logCoarse

          auto prop = theChain.back();

          logCoarseProp = prop.getlogPhi();

          return prop; // Return Proposal from SeqDA

        }

        bool acceptReject(LINK& u,  LINK& v, bool verb = false){

          bool accept = false;


          std::random_device rd;
          std::mt19937 gen(rd());
          std::uniform_real_distribution<> dis(0.0, 1.0);

          double logalpha = std::log(dis(gen));

          double logtestProbability = v.getlogPhi(false) + logCoarse - u.getlogPhi(false) - logCoarseProp;


          if (logalpha < logtestProbability){
            accept = true;
          }

          if(verb){

            std::cout << "logalpha = " << logalpha << std::endl;
            std::cout << "logtestProbability = " << logtestProbability << std::endl;
            std::cout << "Do we accept ? = " << accept << std::endl;
            std::cout << "proposal ll = " << v.getlogPhi(false) << std::endl;
            std::cout << "existing ll = " << u.getlogPhi(false) << std::endl;

          }

          return accept;
        }

        void updateParameters(Eigen::VectorXd & newParam){ param = newParam; }

        Eigen::VectorXd getParameters(){  return param; }

      private:

        Eigen::VectorXd param; // Step size for random step

        int subChain_length;

        PROPOSAL& myProposal;

        ForwardModel& F;

        double logCoarse, logCoarseProp;

        int count, level;

        int adaptive;


    };


}
}


#endif
