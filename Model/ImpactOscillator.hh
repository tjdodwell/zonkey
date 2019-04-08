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

      ImpactOscillator(){

          L = 3; // Max number of Levels

          mu.resize(STOCHASTIC_DIM);

          mu(0) = 0.0; // initial displacement
          mu(1) = 0.1; // initial velocity
          mu(2) = 1.0; // stiffness
          mu(3) = 100.0; // contact stiffness
          mu(4) = 0.1; // Amplitude of osciallations of contact plate
          mu(5) = 1.0; // Frequency of contact plate

          // Diagonal MatrixXd

          double sig = 1.0;

          MatrixXd Sigma = MatrixXd::Identity(6,6);

          Sigma *= sig;

          invSigma = Sigma.inverse();
          C = Sigma.llt().matrixL();

          setLevels();

          numDataPoints = 10;

          totalTime = 20.0;

          double sigf = 0.0229;

          data_time.resize(numDataPoints);
          data_val.resize(numDataPoints);

          // This is the data generated using an adaptive high-fidelity time stepper in Matlab

          double data_step = totalTime / numDataPoints;
          data_time[0] = data_step;
          for (int i = 1; i < numDataPoints; i++){
            data_time[i] = data_time[i-1] + data_step;
          }

          // This is generated from data values, plus some Gaussian Noise

          Eigen::VectorXd gauss_noise(numDataPoints);

          gauss_noise(0) = 0.5377;
          gauss_noise(1) = 1.8339;
          gauss_noise(2) = -2.2588;
          gauss_noise(3) = 0.8622;
          gauss_noise(4) = 0.3188;
          gauss_noise(5) = -1.3077;
          gauss_noise(6) = -0.4336;
          gauss_noise(7) = 0.3426;
          gauss_noise(8) = 3.5784;
          gauss_noise(9) = 2.7694;

          //

          data_val(0) = 0.0908;
          data_val(1) = -0.1594;
          data_val(2) = -0.4049;
          data_val(3) = 0.4645;
          data_val(4) = -0.4233;
          data_val(5) = -0.3248;
          data_val(6) = 0.3536;
          data_val(7) = -0.3671;
          data_val(8) = -0.1754;
          data_val(9) = 0.0643;

          MatrixXd Sigmad = MatrixXd::Identity(numDataPoints, numDataPoints);

          Sigmad *= sigf;

          invSigmad = Sigmad.inverse();

          data_val = data_val + invSigmad * gauss_noise;


      }

      void setLevels(int numSteps0 = 10, double factor = 2.0){
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

        double dt = totalTime / numSteps[level];
        x[1] = initialDisplacement;
        x[0] = x[1] - dt * initialVelocity;

        Eigen::VectorXd F(numDataPoints); // Container for forward model.
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
        double logLikelihood = -0.5 * misMatch.transpose() * ( invSigmad * misMatch );

        // Save information about proposal
        u.setlogPhi(logLikelihood);
        u.setlogPi0(this->logPrior(u));

      }

      double inline forcing(double t, double x, Eigen::VectorXd & xi){
        // This model assumes linear spring and one-sided Hertz Contact
        double offsetOfContactSurface = x - xi(4) * std::sin(xi(5) * t);
        return - xi(2) * x - xi(3) * std::real(std::pow(offsetOfContactSurface,3./2.));
      }

      Eigen::VectorXd grad(Link & u){
        Eigen::VectorXd xi = u.getTheta();
        Eigen::VectorXd grad(2);
        return grad;
      }

  private:

      MatrixXd invSigma, C, invSigmad;
      Eigen::VectorXd mu, data_val;
      double totalTime;
      std::vector<int> numSteps;

      std::vector<double> data_time;

      int L, numDataPoints;

  };

}
}
#endif /* end Zonkey::Models::ImpactOscillator */
