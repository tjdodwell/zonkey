#ifndef ZONKEY_MCMC_SEQ_DRAM_HH
#define ZONKEY_MCMC_SEQ_DRAM_HH

#include "../Chain/Link.hh"


#include <random>
#include <cassert>
#include <cmath>


namespace Zonkey {
  namespace MCMC {

    template<typename LINK>
    class SeqDRAM{

      public:

        SeqDRAM(Eigen::VectorXd param_, Eigen::VectorXd & qcov_):
          param(param_),
          qcov(qcov_){
          assert(param.size() == 2 && "param input wrong length for SeqPCN, should be length 2"); //

          LLT<Eigen::MatrixXd> lltOfA(qcov); // compute the Cholesky decomposition of proposal covariance
          R = lltOfA.matrixL();


        }

        LINK apply(const LINK& currentState){
          Eigen::VectorXd xi = currentState.getTheta();
          Eigen::VectorXd zeta(xi.size());
          std::random_device rd;
          std::normal_distribution<double> dis(0.0,1.0);
          std::mt19937 gen(rd());
          for (int i = 0; i < xi.size(); i++){
            zeta(i) = dis(gen);
          }
          auto xip = xi + zeta * R;     // a new proposal
          LINK prop;
          prop.setTheta(xip);
          return prop; // Return Proposal from SeqRandomWalk
        }

        bool acceptReject(LINK& u,  LINK& v){

          bool accept = false;

          std::random_device rd;
          std::mt19937 gen(rd());
          std::uniform_real_distribution<> dis(0.0, 1.0);

          double logalpha = std::log(dis(gen));

          double logtestProbability = (v.getlogPhi() + v.getlogPi0()) - (u.getlogPhi() + u.getlogPi0());

          if (logalpha < logtestProbability){
            accept = true;
          }

          return accept;
        }

        void updateParameters(Eigen::VectorXd & newParam){ param = newParam; }

        Eigen::VectorXd getParameters(){  return param; }

      private:

        Eigen::VectorXd param; // Step size for random step
        Eigen::VectorXd R, cov;


    };


  }
}


#endif
