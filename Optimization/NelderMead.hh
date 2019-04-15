#ifndef NELDERMEADSOLVER_H_
#define NELDERMEADSOLVER_H_

#include "assert.h"

using namespace Eigen;
using namespace std;

namespace Zonkey{

	namespace Optim{


vector<int> sortIndx(const vector<double> &v) {

  // initialize original index locations
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](int i1, int i2) {return v[i1] < v[i2];});

  return idx;
}

template<class MODEL, int N>
class NelderMead{

	/* NelderMead - Implementation of a Derivative-Free Optimiser

		returns argmin( f.eval(x) ) - x vector of dimension dim

		Requires:

		MODEL - public functions:
			eval(x) - evaluates model f at x
			samplePrior() - returns x from prior distribution

		Optional:

		a, gamma, rho & sig - method parameters

		Example Usage:

		NelderMead<QUADRATIC> myOptimiser(f);

		f.minimize(true);

		--

		For details see - https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method

		--

		Written by Dr T. Dodwell, University of Exeter, t.dodwell@exeter.ac.uk
		2/3/18

	*/

public:

	// Constructors

	NelderMead(MODEL &f_, double a_ = 1.0, double g_ = 2.0, double r_ = 0.5, double s_ =0.5): 
		f(f_), a(a_), gamma(g_), rho(r_), sig(s_){
			xm.resize(N);
			X.resize(N + 1);  F.resize(N + 1); id.resize(N+1);
		}

	// Minimize

	VectorXd maximize(bool verb = false, double tol = 1e-8, int jmax = 100){
		VectorXd x = minimize(verb,tol,jmax,-1.0);
		return x;
	}

	VectorXd minimize(bool verb = false, double tol = 1e-8, int jmax = 100, double factor = 1.0){

		// Initial Evaluation of N samples


		for (int i = 0; i < N + 1; i++){
			X[i] = f.samplePrior();
			F[i] = f.eval(X[i]);
		}	

		reOrder(); // re-order initial samples
		
		j = 0; // initialise number of steps

		if (verb){std::cout << "Starting Nelder-Mead Optimiser ..." << std::endl;}
		

		if(verb){std::cout << "Step j = " << j << ", fval = " << F[id[0]] << " error = " << error() << std::endl;}

		
		maxjExceed = false;

		bool converged = check(tol);

		while(converged == false){ // if converged is false

			computeCentroid();

			// 3) Reflection

			VectorXd xr = xm + a * (xm - X[id[N]]);
			double fxr = factor * f.eval(xr);


			if (fxr < F[id[N-1]] && fxr > F[id[0]]){ 

				X[id[N]] = xr;	F[id[N]] = fxr;

			}
			else{ // F[0] < fxr < F[id[N-1]]
				
				// 4) Expansion
				
				if(fxr < F[id[0]]){ // if the Best so far
					VectorXd xe = xm + gamma * (xr - xm);
					double fxe = factor * f.eval(xe);
					
					if ( fxe < fxr){	X[id[N]] = xe;	F[id[N]] = fxe;	}
					else{	X[id[N]] = xr;	F[id[N]] = fxr;	}
				}
				else{ // 5)  Contraction
					assert(fxr > F[id[N-1]]);
					VectorXd xc = xm + rho * (X[id[N]] - xm);
					double fxc = factor * f.eval(xc);

					if (fxc < F[id[N]]){
						X[id[N]] = xc;	
						F[id[N]] = fxc;}
					else{
						for (int i = 1; i < N + 1; i++){
							X[id[i]] = X[id[0]] + sig * (X[id[i]] - X[id[0]]);
							F[id[i]] = factor * f.eval(X[id[i]]);
						}
					}
			} // end step 4
			
		} // end step 3
		

		reOrder(); // re-order samples

		converged = check(tol); // check if tolerance is met

		j++; // Add step increment

		if(verb){std::cout << "Step j = " << j << ", fval = " << F[id[0]] << " error = " << error() << std::endl;}

		

		if (j > jmax){	
			converged = true;	
			maxjExceed = true;
		} // check if number steps to many

		

	} // end while

	if(maxjExceed){
		std::cout << "*** Max Number of Iterations Exceed, solution not minimised. Error = " << error() << std::endl;
	}
	
	return X[id[0]];

	} // end minimize()

	void printX(){
		for (int i = 0; i < N+1; i++){
			std::cout<< " Trial point " << i << " F = " << F[i] << std::endl;
			std::cout<< " Order of Best " << id[i] << std::endl;
			std::cout  << X[i] << std::endl;
			std::cout << " " << std::endl;
		}
	}

	bool check(double tol){
		bool con = false;
		double err = error();
		if (err < tol){ con = true; }
		return con;
	}

	double error(){
		// Convergence defined according to "Compact Numerical Methods for Computers" 2nd Ed., by John C. Nash
		bool con = false;
		double mean = std::accumulate(F.begin(), F.end(), 0.0) / (N + 1);
		double var = 0.0;
		for (int i = 0; i < N + 1; i++){
			var += (F[i] - mean ) * (F[i] - mean);
		}
		var /= N;
		return std::sqrt(var);
	}

	int numStep(){ return j; }

	bool converged(){ return maxjExceed; }

private:

	void reOrder(){
		id = sortIndx(F);

	}

	void computeCentroid(){
		for (int i =0; i < N; i++){ xm(i) = 0.0; }

		for (int i = 0; i < N; i++){xm += X[id[i]];}
		xm /= N;
	}

	MODEL &f;
	VectorXd xm;
	double a, gamma, rho, sig;
	std::vector<VectorXd> X;
	std::vector<double> F;
	std::vector<int> id;
	int j;
	bool maxjExceed;
	double Fold;
};

} // end Optim

} // end of MarkovPP Namespace

#endif