#ifndef ZONKEY_MCMC_ADAPTIVE_DELAYED_ACCEPTANCE_HH
#define ZONKEY_MCMC_ADAPTIVE_DELAYED_ACCEPTANCE_HH

#include "../MH.hh"

namespace Zonkey {
  namespace MCMC {

    template<typename LINK, typename CHAIN, typename PROPOSAL, typename ForwardModel>
    class SeqAdaptiveDA{

      public:

        SeqAdaptiveDA(PROPOSAL& myProposal_, ForwardModel& F_, int level_ = 1, int subChain_length_ = 10, bool adaptive_ = false) :
          subChain_length(subChain_length_),
          myProposal(myProposal_),
          F(F_),
          level(level_),
          adaptive(adaptive_){

            count = 1;

            int dimensionData = F.getDIMData();

            muB.resize(dimensionData);

            SigmaE = F.getSigma()

            SigmaB.resize(dimensionData,  dimensionData);

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

          count += 1;


          std::random_device rd;
          std::mt19937 gen(rd());
          std::uniform_real_distribution<> dis(0.0, 1.0);

          double logalpha = std::log(dis(gen));

          double logtestProbability = v.getlogPhi(false) + logCoarse - u.getlogPhi(false) - logCoarseProp;


          if (logalpha < logtestProbability){
            accept = true;
          }

          if(accept){

            // Update Mean and Covariance - According to Equation Equation (7) - Cui et al. Adaptive delayed acceptance MH algorithm

            Eigen::VectorXd B = v.getB();

            mu = (1 / (count+1)) * (count * muB + B);

            SigmaB = (1/count) * ( (count - 1) * SigmaB + (B * B.transpose()) - count * (muB * muB.transpose()))

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

        Eigen::VectorXd muB;

        Eigen::MatrixXd SigmaE, SigmaB;


    };


}
}


#endif
