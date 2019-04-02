#ifndef ZONKEY_MCMC_ROSENBROCK_HH
#define ZONKEY_MCMC_ROSENBROCK_HH


#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

namespace Zonkey {
  namespace Models {

  template<typename LINK>
  class Rosenbrock{

  public:

      STOCHASTIC_DIM = 2;

      Rosenbrock(double a_ = 1.0, double b_ = 100.0, double mean = 0.0, double sig = 1.0) :
        a(a_),
        b(b_){
          Eigen::VectorXd mu(STOCHASTIC_DIM);
          mu(0) = mean; mu(1) = mean;
          Matrix2d Sigma;
          Sigam << sig, 0,
                   0, sig;
          invSigma = Sigma.inverse()
          C = llt.compute(Sigma);
      }

      double logPrior(Link & u){
        Eigen::VectorXd xi = u.getTheta();
        Eigen::VectorXd x = xi - mu;
        return -0.5 * (x.tranpose() * invSigma) * x
      } // sample logPrior

      Eigen::VectorXd samplePrior(){
        std::random_device rd;
        std::normal_distribution<double> dis(0.0,1.0);
        std::mt19937 gen(rd());
        Eigen::VectorXd z(STOCHASTIC_DIM);
        for (int i = 0; i < STOCHASTIC_DIM; i++){
          z(i) = dis(gen);
        }
        return C * z + mu; // Sample from prior - note Sigma = C' * C
      } // samplePrior

      double logDensity(LINK & u){
        Eigen::VectorXd xi = u.getTheta();
        return pow((a - xi[0]),2) + b * pow(xi[1] - xi[0] * x[0],2);
      }

      Eigen::VectorXd grad(LINK & u){
        Eigen::VectorXd xi = u.getTheta();
        Eigen::VectorXd grad(2);
        grad[0] = 2.0 * (a - xi[0]) + 4.0 * b * (xi[1] - xi[0] * xi[0]) * xi[0];
        grad[1] = 2.0 * b * (xi[1] - xi[0] * xi[0]);
      }

  private:
    double a, b; // Parameters
    Matrix2d invSigma, C;
  };

}
}
#endif /* end Zonkey::Models::Rosenbrock */
