#ifndef ZONKEY_MCMC_SEQ_TWO_LEVEL_HH
#define ZONKEY_MCMC_SEQ_TWO_LEVEL_HH

#include "../Chain/Link.hh"


#include <random>
#include <cassert>
#include <cmath>


namespace Zonkey {
  namespace MCMC {

    template<typename LINK, typename CHAIN>
    class SeqTwoLevel{

      public:

        SeqTwoLevel(CHAIN & coarseChain_, Eigen::VectorXd & param_,int subSampleRate_ = 1):
          param(param_),
          subSampleRate(subSampleRate_),
          coarseChain(coarseChain_){

            // Get Number of Parameters in Coarse Chain Links

            dim_coarse = coarseChain[0].getSDIM();

            LINK hackforDIM;

            dim_fine = hackforDIM.getSDIM();

            subSampledCount = subSampleRate;
        
        }

        LINK apply(const LINK& currentState){
          
          Eigen::VectorXd xi = currentState.getTheta(); // This is the fine state

          Eigen::VectorXd xip(xi.size()); // This is a container for the proposal.

          auto coarseModes = coarseChain[subSampledCount].getTheta();

          //std::cout << subSampledCount << std::endl;

          std::random_device rd;
          std::normal_distribution<double> dis(0.0,1.0);
          std::mt19937 gen(rd());

          for (int i = 0; i < xi.size(); i++){
          //  if(i < dim_coarse){
              xip(i) =  coarseModes(i);
          /*  }
            else{
              std::cout << "Should not get here" << std::endl;
              xip(i) = std::sqrt(1 - param(0) * param(0)) * xi(i) + param(0) * sig * dis(gen);
            }*/
          }

          LINK prop;
          prop.setTheta(xip);

          Eigen::VectorXd Qc = coarseChain[subSampledCount].getQ();

          prop.setQoI(Qc,true);

          return prop; // Return Proposal
        }

        bool acceptReject(LINK& u,  LINK& v, bool verb = false){

          bool accept = false;

          std::random_device rd;
          std::mt19937 gen(rd());
          std::uniform_real_distribution<> dis(0.0, 1.0);

          double logUniform = std::log(dis(gen));

          auto vC = coarseChain[subSampledCount];

          auto uC = coarseChain[subSampledCount - subSampleRate]; 

          // I need to update this for MLMC accept reject step.

          double logalpha = v.getlogPhi() + uC.getlogPhi() - u.getlogPhi() - vC.getlogPhi();

          if (logUniform < logalpha){
            accept = true;
          }

          if(verb){

            std::cout << "logalpha = " << logalpha << std::endl;
            //std::cout << "logtestProbability = " << logtestProbability << std::endl;
            std::cout << "Do we accept ? = " << accept << std::endl;
            std::cout << "proposal ll = " << v.getlogPhi() << std::endl;
            std::cout << "existing ll = " << u.getlogPhi() << std::endl;

          }

          subSampledCount += subSampleRate;

          return accept;
        }

        void setCounter(int newCount){  counter = newCount;}

        double getScaling(){ return sig;}

        void updateParameters(Eigen::VectorXd & newParam){ param = newParam; }

        Eigen::VectorXd getParameters(){  return param; }

      private:

        Eigen::VectorXd param; // Step size for random step

        double optimal_alpha, sig, log_sig_new;

        int counter;

        double gamma_scaling;

        bool optimize_scaling;

        CHAIN & coarseChain;

        int subSampleRate;

        int dim_coarse, dim_fine, subSampledCount;


    };


  }
}


#endif
