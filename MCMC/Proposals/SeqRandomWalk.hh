#ifndef ZONKEY_MCMC_SEQ_RANDOM_WALK_HH
#define ZONKEY_MCMC_SEQ_RANDOM_WALK_HH

#include "../Chain/Link.hh"


#include <random>
#include <cassert>
#include <cmath>


namespace Zonkey {
  namespace MCMC {

    template<typename LINK>
    class SeqRandomWalk{

      public:

        SeqRandomWalk(Eigen::VectorXd param_):param(param_){
          assert(param.size() == 1 && "param input wrong length for SeqRandomWalk, should be length 1"); //
        }

        LINK apply(const LINK& currentState){
          Eigen::VectorXd xi = currentState.getTheta();
          Eigen::VectorXd xip(xi.size());
          std::random_device rd;
          std::normal_distribution<double> dis(0.0,1.0);
          std::mt19937 gen(rd());
          for (int i = 0; i < xi.size(); i++){
            xip(i) = xi(i) + param(0) * dis(gen);
          }
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


    };


  }
}


#endif
