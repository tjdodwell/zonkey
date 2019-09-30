#ifndef ZONKEY_MCMC_HIERARCHICAL_DELAYED_ACCEPTANCE_HH
#define ZONKEY_MCMC_HIERARCHICAL_DELAYED_ACCEPTANCE_HH

#include "../MH.hh"

#inlcuyde "SeqDA.hh"

namespace Zonkey {
  namespace MCMC {

    template<typename LINK, typename CHAIN, typename PROPOSAL, typename ForwardModel>
    class SeqHIERARCHICAL_DA{

      public:

        SeqHIERARCHICAL_DA(PROPOSAL& myProposal_, std::vector<ForwardModel> & F_, int L_, std::vector<double>& subChain_lengths) :
          subChain_lengths(subChain_lengths_),
          myProposal(myProposal_),
          F(F_),
          L(L_){

            count = 0;
            sig_pcn = myProposal.getScaling();


        }

        LINK apply(LINK& currentState){

          Eigen::VectorXd xi = currentState.getTheta();


          



          MetropolisHastings<LINK,CHAIN,,ForwardModel> myMCMC(F,myProposal,markovChain);


          if

          MetropolisHastings<LINK,CHAIN,PROPOSAL,ForwardModel> myMCMC(F,myProposal,markovChain);

          myMCMC.setStart(xi,level-1);

          myMCMC.run(subChain_length,level-1,"",false);

          auto theChain = myMCMC.getChain();

          logCoarse = theChain[0].getlogPhi();
          myMCMC.run(subChain_length-1,level-1,"",false);

          auto mc = myMCMC.getChain();

          auto prop = mc.back();
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
            std::cout << "proposal ll = " << v.getlogPhi() << std::endl;
            std::cout << "existing ll = " << u.getlogPhi() << std::endl;

          }

          return accept;
        }


      private:


        PROPOSAL& myProposal;

        std::vector<ForwardModel>& F;

        double logCoarse, logCoarseProp;

        int count;

        int L; // Maximum number of Levels


    };


}
}


#endif
