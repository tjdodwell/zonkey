#ifndef ZONKEY_MCMC_ROSENBROCK_HH
#define ZONKEY_MCMC_ROSENBROCK_HH


using namespace Eigen;
using namespace std;

namespace Zonkey {
  namespace Models {

  class Rosenbrock{

  public:

      Rosenbrock(double a_ = 1.0, double b_ = 100.0) :
        a(a_),
        b(b_){    }

      double apply(const Link & u){
        Eigen::VectorXd xi = u.getTheta();
        return pow((a - xi[0]),2) + b * pow(xi[1] - xi[0] * x[0],2);
      }

      Eigen::VectorXd grad(const Link & u){
        Eigen::VectorXd xi = u.getTheta();
        Eigen::VectorXd grad(2);
        grad[0] = 2.0 * (a - xi[0]) + 4.0 * b * (xi[1] - xi[0] * xi[0]) * xi[0];
        grad[1] = 2.0 * b * (xi[1] - xi[0] * xi[0]);
      }

  private:

    double a, b; // Parameters
    

  };

}
}
#endif /* chain_h */
