
#include <random>

#include "KLFunctions.hh"

template<int DIM>
class RandomField{

public:

    RandomField(){

      std::cout << "Initialising Random Field!...." << std::endl;

    }

    void setup(double L_ = 1.0, double sigKl_ = 1.0, double correlation_length_ = 0.7, int S_Dim = 100){

      L = L_;
      sigKl = sigKl_;
      correlation_length = correlation_length_;
      Stochastic_Dim = S_Dim;

      // Computed Quantities
      xi.resize(Stochastic_Dim);
      lambda.resize(Stochastic_Dim);

      oneD_Stochastic_dim = std::ceil(std::pow(Stochastic_Dim,(double) 1./DIM));

      std::cout << oneD_Stochastic_dim << std::endl;

      freq.resize(oneD_Stochastic_dim); lam1D.resize(oneD_Stochastic_dim);

      rootFinder(oneD_Stochastic_dim, correlation_length / L, freq);

      evaluate1DeigenValues(correlation_length / L, lam1D, freq);


      index.resize(DIM);

      if (DIM == 2){
        construct2DeigenValues(lam1D,lambda,index);
      }
      else{
        construct3DeigenValues(lam1D,lambda,index);
      }




    }

    int getStochasticDim(){return Stochastic_Dim; }

    void writeToVtk(){


    }

    void setXi(std::vector<double>& xi_){
      for (int i = 0; i < xi_.size(); i++){ xi[i] = xi_[i]; }
    }

    void inline sample_prior(int rank, int seed = 1){

      using Dune::PDELab::Backend::native;

      // user_define_prior return a single (independent) sample from prior distribution
        double * xi_local = new double[Stochastic_Dim];
        if(rank == 0){
          std::random_device rd;
          std::normal_distribution<double> dis(0.0,1.0);
          std::mt19937 gen(seed);
          for (int k = 0; k < Stochastic_Dim; k++){  xi_local[k] = dis(gen);  }
        } // end if rank 0
         MPI_Bcast(xi_local,Stochastic_Dim,MPI_DOUBLE,0,MPI_COMM_WORLD);
        //Eigen::VectorXd xi = Eigen::VectorXd::Zero(Stochastic_Dim);
        for (int k = 0; k < Stochastic_Dim; k++){  xi[k] = xi_local[k];  }
        delete xi_local;
    }

    double inline evalPhi(double x, int i){
      double phi = 0.0;
      double omega = freq[i];
      x /= L;
      double l = correlation_length / L;

      double norm = std::sin(2.0*omega)*(0.25*(l*l*omega - 1.0/omega)) - 0.5 * l * std::cos(2.0 * omega) + 0.5 * ( 1 + l + l*l*omega*omega);
      norm = 1.0 / std::sqrt(norm);
      phi = norm * (std::sin(x * omega) + l * omega * std::cos(x * omega) );
      return phi;
    }

    double inline getPerm(const Dune::FieldVector<double,DIM> & x){
      // Compute Random Field at a point
        double perm = 0.0;
        for (int j = 0; j < Stochastic_Dim; j++){
            double contribution2perm = sigKl * std::sqrt(lambda[j])*xi[j];

            for (int i = 0; i < DIM; i++){
              contribution2perm *= evalPhi(x[i],index[i][j]);
            }
            perm += contribution2perm;
        }
        return std::exp(perm); // Return random_field
    }



private:

    double L, sigKl, correlation_length;
    int Stochastic_Dim, oneD_Stochastic_dim;
    std::vector<std::vector<int>> index;
    std::vector<double> lambda, xi ,freq, lam1D;


};
