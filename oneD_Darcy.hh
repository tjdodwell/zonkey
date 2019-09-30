#ifndef ZONKEY_MCMC_ONED_DARCY_HH
#define ZONKEY_MCMC_ONED_DARCY_HH

// TO DO: ADD BASE CLASS FOR MODEL IN WHICH EACH MODEL INHERITS CORE FUNCTIONS

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include <random>

namespace Zonkey {
  namespace Models {

  template<typename Link, int STOCHASTIC_DIM>
  class oneD_Darcy{

  public:

      oneD_Darcy(int nelem = 10, double mean_ = 0.0, double sig_ = 1.0, bool isParallel = true) :
        mean(mean_),
        sig(sig_){
          mu.resize(nelem);
          for (int i = 0; i < nelem; i++){
            mu(i) = mean;
          }

          Matrix2d Sigma = Matrix2d::Identity(nelem,nelem);

          Sigma *= sig;

          invSigma = Sigma.inverse();
          C = Sigma.llt().matrixL();
      }

      double logPrior(Link & u){
        Eigen::VectorXd xi = u.getTheta();
        Eigen::VectorXd x = xi - mu;
        return -0.5 * (x.transpose() * invSigma) * x;
      } // sample logPrior

      Eigen::VectorXd samplePrior(int level = 0){
        std::random_device rd;
        std::normal_distribution<double> dis(0.0,1.0);
        std::mt19937 gen(rd());
        Eigen::VectorXd z(STOCHASTIC_DIM);
        for (int i = 0; i < STOCHASTIC_DIM; i++){
          z(i) = dis(gen);
        }
        return C * z + mu; // Sample from prior - note Sigma = C' * C

      } // samplePrior

      void apply(Link & u, int level = 0){
        Eigen::VectorXd xi = u.getTheta();

        // This is where the model will go



        u.setlogPhi(-ld);
        u.setlogPi0(this->logPrior(u));
      }

      void buildSnapShots(){



      }






  private:
    double a, b; // Parameters
    Matrix2d invSigma, C;
    Eigen::VectorXd mu;
  };

}
}
#endif /* end Zonkey::Models::Rosenbrock */
