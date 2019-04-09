#ifndef ZONKEY_MCMC_SHM_HH
#define ZONKEY_MCMC_SHM_HH


#include <Eigen/Dense>

#include <boost/lambda/lambda.hpp>

#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>

#include "Impact_ode.hh"

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

          MatrixXd Sigma = MatrixXd::Zero(6,6);

          Sigma(0,0) = 0.02;
          Sigma(1,1) = 0.02;
          Sigma(2,2) = 0.1;
          Sigma(3,3) = 1.0;
          Sigma(4,4) = 0.02;
          Sigma(5,5) = 0.1;


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

          // data value

          data_val(0) = 0.0419;
          data_val(1) = -0.1127;
          data_val(2) = -0.3168;
          data_val(3) = 0.0383;
          data_val(4) = -0.0244;
          data_val(5) = -0.3644;
          data_val(6) = 0.0182;
          data_val(7) = 0.0338;
          data_val(8) = -0.4053;
          data_val(9) = -0.0163;

          MatrixXd Sigmad = MatrixXd::Identity(numDataPoints, numDataPoints);

          Sigmad *= sigf;

          invSigmad = Sigmad.inverse();

          data_val = data_val + invSigmad * gauss_noise;


      }

      void setLevels(int numSteps0 = 100, double factor = 2.0){
        numSteps.resize(L);
        numSteps[0] = numSteps0;
        for (int i = 1; i < numSteps.size(); i++){
          numSteps[i] = numSteps[i-1] * factor;
        }

      }

      void testRun(Link & u){

        // Multilevel Model for Impact Oscillator
        Eigen::VectorXd xi = u.getTheta();

        using namespace std;
    using namespace boost::numeric::odeint;


      //[ state_initialization
      state_type x(2);
      x[0] = xi(0); // start at x=1.0, p=0.0
      x[1] = xi(1);
      //]

    //[ integrate_observ
    vector<state_type> x_vec;
    vector<double> times;


    //[ define_const_stepper
    runge_kutta4< state_type > stepper;
    size_t steps = integrate_const( stepper , harmonic_oscillator , x , 0.0 , 10.0 , 0.1, push_back_state_and_time( x_vec , times ) );

    /* output */
    for( size_t i=0; i<=steps; i++ )
    {
        cout << times[i] << '\t' << x_vec[i][0] << '\t' << x_vec[i][1] << '\n';
    }
    //]
    //]


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

        std::cout << "I am sampling the prior = " << C * z + mu << std::endl;
      } // samplePrior

      void apply(Link & u, int level = 0, bool plotSolution = false){

        // Multilevel Model for Impact Oscillator
        Eigen::VectorXd xi = u.getTheta();

        double logLikelihood = 0.0;

        if (xi(2) < 0 || xi(3) < 0){ // First we need to check the bounds on proposal
          logLikelihood = 10e9;
        }
        else{

        //

            std::vector<double> x(numSteps[level] + 2), time(numSteps[level] + 2);

            double initialDisplacement = xi(0);
            double initialVelocity = xi(1);

            double dt = totalTime / numSteps[level];
            x[1] = initialDisplacement;
            x[0] = x[1] - dt * initialVelocity;

            Eigen::VectorXd F(numDataPoints); // Container for forward model.
            time[0] = -dt;
            time[1] = 0.0;
            int k = 0;
            for (int t = 2; t < numSteps[level] + 2; t++ ){
              time[t] = time[t-1] + dt;
              double v = (x[t-1] - x[t-2]) / dt;
              x[t] = 2 * x[t-1] - x[t-2] + dt * dt * forcing(time[t-1],v,x[t-1],xi);
              if ( (time[t-1] < data_time[k]) && (time[t] > data_time[k]) ){
                // There is a data recording within this time step
                double data_dt = data_time[k] - x[t-1]; // this is the time step to the data point
                F(k) = 2 * x[t-1] - x[t-2] + data_dt * data_dt * forcing(time[t-1],v,x[t-1],xi);
                k += 1;
              } // Record forward Model
            } // loop through time step

            // Compute logLikelihood
            Eigen::VectorXd misMatch(numDataPoints);
            for (int k = 0; k < numDataPoints; k++){  misMatch(k) = (F(k) - data_val(k)) * (F(k) - data_val(k)); }
            logLikelihood = -0.5 * misMatch.transpose() * ( invSigmad * misMatch );

            if(plotSolution){

              std::cout << "=== Time" << std::endl;
              for (int i = 1; i < numSteps[level] + 2; i++){
                std::cout << time[i] << std::endl;
              }
              std::cout << "=== F" << std::endl;
              for (int i = 1; i < numSteps[level] + 2; i++){  std::cout << x[i] << std::endl; }
            }

        }

        // Save information about proposal
        u.setlogPhi(logLikelihood);
        u.setlogPi0(this->logPrior(u));

      }

      double inline forcing(double t, double x, double v, Eigen::VectorXd & xi){
        // This model assumes linear spring and one-sided Hertz Contact
        double offsetOfContactSurface = x - xi(4) * std::sin(xi(5) * t);
        return - xi(2) * x;// - 100.0 * std::real(std::pow(offsetOfContactSurface,3./2.)) - v;
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
