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



Notes Greedy Start on update Covariance on those samples which are accepted.
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

          alpha_star = 0.25;

          eps = -0.5;

          C = sd * MatrixXd::Identity(STOCHASTIC_DIM,STOCHASTIC_DIM);

        }

        void updateC(LINK& newLink){

          /*

          TH(:,i) = X;
        %Update scale
        idx = (i+(m-1)*mrun+burnin);    % Global index
        gamma=idx^eps;
        s =exp(log(s)+gamma*(acc_rate-optimal));
        avg_acc = avg_acc*(i-1)/i + acc_rate/i;
        %Update Sigma
        delta_X=X-X_bar;    X_bar = X_bar+gamma*delta_X;
        Sdelta=gamma*(delta_X*delta_X')-((idx-1)^eps)*V;
        V = V+Sdelta;   L = tril(V,-1);V=diag(diag(V))+L+L';
        [V2,errV]=chol(V,'lower');

          */

          int k = xBar.size()

          Cold = C;

          // Calculate bar{X}_k

          Eigen::VectorXd Xk = newLink.getTheta();

          if (k < 2 * STOCHASTIC_DIM){
            if (k < 2){
              xBar.push_back(Xk);
            }
            else{
              xBar.push_back((1./(k+1)) * (k * xBar.back()  + Xk));
            }
          }
          else{

            xBar.push_back((1./(k+1)) * (k * xBar.back()  + Xk));

            Eigen::MatrixXd mat1 = xBar[k-2].transpose() * xBar[k-2];

            Eigen::MatrixXd mat2 = xBar[k-1].transpose() * xBar[k-1];

            Eigen::MatrixXd mat3 = Xk.transpose() * Xk;

            Eigen::MatrixXd mat4 = epsilon * MatrixXd::Identity(STOCHASTIC_DIM,STOCHASTIC_DIM);

            C = ((k-1) / k) * Cold + (sd / k) * ( k * mat1 - (k + 1) * mat2 + mat3 + mat4);

          }

        }

        LINK apply(const LINK& currentState){

          Eigen::VectorXd xi = currentState.getTheta();

          updateC(currentState);

          // Make Proposal
          Eigen::VectorXd xip(xi.size());
          std::random_device rd;
          std::normal_distribution<double> dis(0.0,1.0);
          std::mt19937 gen(rd());
          Eigen::VectorXd gauss(xi.size());
          for (int i = 0; i < xi.size(); i++){
            gauss(i) =  dis(gen);
          }
          // Compute Chol of  C
          Eigen::MatrixXd L =  C.llt().matrixL();
          xip = xi  + std::sqrt(sd) * (L * gauss);
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

          // In this step we update scaling.

          // Update Scaling for proposal

          double gamma = (double) std::pow(xBar.size(),eps);

          double log_sd = std::log(sd) + gamma * (std::exp(logalpha) - alpha_star);

          sd = std::exp(sd);

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

        Eigen::MatrixXd C, Cnew;


        std::vector<Eigen::VectorXd> Xbar;

        double sd, epsilon, eps, alpha_star;


    };


  }
}


#endif
