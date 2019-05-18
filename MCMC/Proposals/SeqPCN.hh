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

        SeqPCN(Eigen::VectorXd & param_, bool optimize_scaling_ = true):
          param(param_),
          optimize_scaling(optimize_scaling_){
          assert(param.size() == 2 && "param input wrong length for SeqPCN, should be length 2"); //a

          std::cout << "Setting up P" << std::endl;

          gamma_scaling = -0.5;

          optimal_alpha = 0.235;

          sig = param(1);

          counter = 0;

          std::cout << counter << std::endl;



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
          }

          // Optimal Scaling


          if (optimize_scaling){

            counter += 1;

            double alpha = std::min(std::exp(logalpha),1.0);

            double gamma = std::pow((double) counter, gamma_scaling);

            log_sig_new = std::log(sig) + gamma * (alpha - optimal_alpha);

            std::cout << log_sig_new << std::endl;

            sig = std::exp(log_sig_new);

          }

          if(verb){

            std::cout << "logalpha = " << logalpha << std::endl;
            //std::cout << "logtestProbability = " << logtestProbability << std::endl;
            std::cout << "Do we accept ? = " << accept << std::endl;
            std::cout << "proposal ll = " << v.getlogPhi() << std::endl;
            std::cout << "existing ll = " << u.getlogPhi() << std::endl;

          }

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


    };


  }
}


#endif
