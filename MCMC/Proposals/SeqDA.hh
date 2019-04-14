#ifndef ZONKEY_MCMC_DELAYED_ACCEPTANCE_HH
#define ZONKEY_MCMC_DELAYED_ACCEPTANCE_HH

namespace Zonkey {
  namespace MCMC {

    template<typename Link, typename Chain, typename PROPOSAL, typename ForwardModel>
    class SeqDA{

      public:

        SeqDA(PROPOSAL& myProposal_, ForwardModel& F_, int subChain_length_ = 10) :
          subChain_length(subChain_length_),
          myProposal(myProposal_),
          F(F_){


        }

        LINK apply(LINK& currentState){
          Eigen::VectorXd xi = currentState.getTheta();
          CHAIN markovChain; // Setup a sub chain
          Zonkey::MCMC::MetropolisHastings<LINK,CHAIN,PROPOSAL,ForwardModel> myMCMC(F,myProposal,markovChain);
          myMCMC.setStart(xi,0);
          myMCMC.run(subChain_length,0);
          LINK prop = myMCMC.back()
          return prop; // Return Proposal from SeqDA
        }

        bool acceptReject(LINK& u,  LINK& v, bool verb = false){

          bool accept = false;

          std::random_device rd;
          std::mt19937 gen(rd());
          std::uniform_real_distribution<> dis(0.0, 1.0);

          double logalpha = std::log(dis(gen));

          double logtestProbability = v.getlogPhi(false) + u.getlogPhi(true) - u.getlogPhi(false) - v.getlogPhi(true);

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

        void updateParameters(Eigen::VectorXd & newParam){ param = newParam; }

        Eigen::VectorXd getParameters(){  return param; }

      private:

        Eigen::VectorXd param; // Step size for random step

        int subChain_length;

        PROPOSAL& myProposal;
        ForwardModel& F;


    };





#endif
