#ifndef ZONKEY_MCMC_SHM_HH
#define ZONKEY_MCMC_SHM_HH


#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include <random>

namespace Zonkey {
  namespace Models {

  template<typename Link, int STOCHASTIC_DIM>
  class ImpactOscillator{

  public:

      ImpactOscillator(double a_ = 1.0, double b_ = 100.0, double mean = 0.0, double sig = 1.0) :
        a(a_),
        b(b_){
          mu.resize(STOCHASTIC_DIM);
          mu(0) = mean; mu(1) = mean;
          Matrix2d Sigma;
          Sigma << sig, 0,
                   0, sig;
          invSigma = Sigma.inverse();
          C = Sigma.llt().matrixL();

          setLevels();

          int numDataPoints = 10;

          totalTime = 1.0;

          double sigf = 1.0;

          data_time.resize(numDataPoints);
          data_val.resize(numDataPoints);

          // This is the data generated using an adaptive high-fidelity time stepper





          // This is generated from data values, plus some Gaussian Noise






      }

      void setLevels(int numSteps0 = 10, double factor = 2.0;){
        numSteps.resize(L);
        numSteps[0] = numSteps0;
        for (int i = 1; i < numSteps.size(); i++){
          numSteps[i] = numSteps[i-1] * factor;
        }

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

        // Multilevel Model for Impact Oscillator
        Eigen::VectorXd xi = u.getTheta();
        std::vector<double> x(numSteps + 2), time(numSteps + 2);

        double initialDisplacement = xi(0);
        double initialVelocity = xi(1);

        dt = totalTime / numSteps[level];
        x[1] = initialDisplacement;
        x[0] = x[1] - dt * initialVelocity;

        Eigen::VectorXd F(data_points); // Container for forward model.
        time[0] = -dt;
        time[1] = 0.0;
        int k = 0;
        for (int t = 2; t < numSteps + 2; t++ ){
          time[t] = time[t-1] + dt;
          x[t] = 2 * x[t-1] - x[t-2] + dt * dt * forcing(time[t-1],x[t-1],xi);
          if ( (time[t-1] < data_time[k]) && (time[t] > data_time[k]) ){
            // There is a data recording within this time step
            double data_dt = data_time[k] - x[t-1]; // this is the time step to the data point
            F(k) = 2 * x[t-1] - x[t-2] + data_dt * data_dt * forcing(time[t-1],x[t-1],xi);
            k += 1;
          } // Record forward Model
        } // loop through time step

        // Compute logLikelihood
        Eigen::VectorXd misMatch(numDataPoints);
        for (int k = 0; k < numDataPoints; k++){  misMatch(k) = (F(k) - data_val[k]) * (F(k) - data_val[k]); }
        double logLikelihood = -0.5 * misMatch.transpose() * ( invSigma * misMatch );

        // Save information about proposal
        u.setlogPhi(logLikelihood);
        u.setlogPi0(this->logPrior(u));

      }

      double inline forcing(double t, double x, Eigen::VectorXd & xi){
        // This model assumes linear spring and one-sided Hertz Contact
        return xi(2) * x - xi(3) * std::pow(x,3./2.);
      }

      Eigen::VectorXd grad(Link & u){
        Eigen::VectorXd xi = u.getTheta();
        Eigen::VectorXd grad(2);
        return grad;
      }

  private:

      double a, b; // Parameters
      Matrix2d invSigma, C;
      Eigen::VectorXd mu;

      double totalTime;

      std::vector<int> numSteps;
  };

}
}
#endif /* end Zonkey::Models::ImpactOscillator */
