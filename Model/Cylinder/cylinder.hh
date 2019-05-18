


#include <iostream>
#include <fstream>
#include <random>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace Zonkey {
  namespace Models {


template <typename Link, int STOCHASTIC_DIM, typename GRID>
class Cylinder{

  public:

    Cylinder(GRID& grid_): grid(grid_){

      

    }


    Eigen::VectorXd samplePrior(int level = 0){
      std::random_device rd;
      std::normal_distribution<double> dis(0.0,1.0);
      std::mt19937 gen(rd());
      Eigen::VectorXd z(STOCHASTIC_DIM);
      for (int i = 0; i < STOCHASTIC_DIM; i++){
        z(i) = dis(gen);
      }
      return z; // Sample from prior - note Sigma = C' * C
    } // samplePrior


    void apply(Link& u, int level = 0, bool plotSolution = false, bool setasData = false){

 

    }


  private:

    Eigen::VectorXd thetaObs, Fobs;
    Eigen::MatrixXd obsCoord;

    int Nobs;

    double sigf;

    GRID& grid;

};

}
}
