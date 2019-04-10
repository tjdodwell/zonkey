#ifndef ZONKEY_MCMC_SEQ_AM_HH
#define ZONKEY_MCMC_SEQ_AM_HH

#include "../Chain/Link.hh"


#include <random>
#include <cassert>
#include <cmath>

/*

Adaptive MCMC

In AM the covariance of the proposal distribution depends is buitl from the history of the chain.

After an initial non-adaptation period, the proposal to be centered at the current position of the Markov chain, X_t, and we set the covariance to be

C_t = s_d Cov(X_0, ... , X_d) + s_d epsilon I_d

where sd is a parameter that depends only on the dimension d of the state space where π is defined

ε > 0 is a constant that we may choose very small compared to the size of S.

Here I_d denotes the d-dimensional identity matrix.

In order to start the adaptation procedure an arbitrary strictly positive definite initial covariance, C0, is chosen according to a priori knowledge (which may be quite poor).

*/


namespace Zonkey {
  namespace MCMC {

    template<int STOCHASTIC_DIM, typename LINK>
    class SeqAM{

      public:

        SeqAM(Eigen::VectorXd param_):param(param_){

          assert(param.size() == 2 && "param input wrong length for SeqPCN, should be length 2"); //

          sd = (double) (2.4 * 2.4) / STOCHASTIC_DIM;

          epsilon = 0.0001;

        }

        void computeCov(LINK& newLink){

          Eigen::Vector Xk = newLink.getTheta();

          xBar.push_back((1./(k+1)) * (k * xBar.back()  + Xk));




        }

        LINK apply(const LINK& currentState){
          Eigen::VectorXd xi = currentState.getTheta();
          Eigen::VectorXd xip(xi.size());
          std::random_device rd;
          std::normal_distribution<double> dis(0.0,1.0);
          std::mt19937 gen(rd());
          for (int i = 0; i < xi.size(); i++){
            xip(i) = std::sqrt(1 - param(0) * param(0)) * xi(i) + param(0) * param(1) * dis(gen);
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

          double logalpha = std::log(dis(gen));

          double logtestProbability = v.getlogPhi() - u.getlogPhi();

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


        double sd, epsilon;


    };


  }
}


#endif
