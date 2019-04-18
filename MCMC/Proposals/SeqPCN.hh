#ifndef ZONKEY_MCMC_SEQ_PCN_HH
#define ZONKEY_MCMC_SEQ_PCN_HH

#include "../Chain/Link.hh"


#include <random>
#include <cassert>
#include <cmath>


namespace Zonkey {
  namespace MCMC {

    template<typename LINK>
    class SeqPCN{

      public:

        SeqPCN(Eigen::VectorXd & param_):param(param_){
          assert(param.size() == 2 && "param input wrong length for SeqPCN, should be length 2"); //

          gamma_scaling = -0.5;

          optimal_alpha = 0.235;

          sig = param(1);

          counter = 0;



        }

        LINK apply(const LINK& currentState){
          Eigen::VectorXd xi = currentState.getTheta();
          Eigen::VectorXd xip(xi.size());
          std::random_device rd;
          std::normal_distribution<double> dis(0.0,1.0);
          std::mt19937 gen(rd());
          for (int i = 0; i < xi.size(); i++){
            xip(i) = std::sqrt(1 - param(0) * param(0)) * xi(i) + param(0) * sig * dis(gen);
          }
          LINK prop;
          prop.setTheta(xip);
          return prop; // Return Proposal from SeqRandomWalk
        }

        bool acceptReject(LINK& u,  LINK& v, bool verb = false){

          bool accept = false;

          std::random_device rd;
          std::mt19937 gen(rd());
          std::uniform_real_distribution<> dis(0.0, 1.0);

          double logUniform = std::log(dis(gen));

          double logalpha = v.getlogPhi() - u.getlogPhi();

          if (logUniform < logalpha){
            accept = true;

            std::cout << "We ACCEPT! " << std::endl;
          }

          // Optimal Scaling

          counter += 1;

          double alpha = std::min(std::exp(logalpha),1.0);

          double gamma = std::pow((double) counter, gamma_scaling);

          double log_sig_new = std::log(sig) + gamma * (alpha - optimal_alpha);

          sig = std::exp(log_sig_new);

          std::cout << alpha << std::endl;


          if(verb){

            std::cout << "logalpha = " << logalpha << std::endl;
            //std::cout << "logtestProbability = " << logtestProbability << std::endl;
            std::cout << "Do we accept ? = " << accept << std::endl;
            std::cout << "proposal ll = " << v.getlogPhi() << std::endl;
            std::cout << "existing ll = " << u.getlogPhi() << std::endl;

          }

          return accept;
        }

        double getScaling(){return sig;}

        void updateParameters(Eigen::VectorXd & newParam){ param = newParam; }

        Eigen::VectorXd getParameters(){  return param; }

      private:

        Eigen::VectorXd param; // Step size for random step

        double optimal_alpha, sig, log_sig_new;

        int counter;

        double gamma_scaling;


    };


  }
}


#endif
