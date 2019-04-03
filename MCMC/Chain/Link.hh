#ifndef ZONKEY_MCMC_LINK_HH
#define ZONKEY_MCMC_LINK_HH


namespace Zonkey {
  namespace MCMC {


    #include <assert.h>

using namespace Eigen;
using namespace std;

template<int STOCHASTIC_DIM, int numQoI>
class Link{

	// Defines a Class which contains single link in Markov Chain
	// Link<STOCHASTIC_DIM,numQoI> mylink()
	// or in SEQ_MarkovChain
	// std::vector<Link<STOCHASTIC_DIM,numQoI>> chain(10) - generates a chain of 10 links, in which each call the default constructor.

	// May need to contatin more information i.e. gradients??
public:

	int accept;


	Link(): accept(0), theta(STOCHASTIC_DIM), Q(numQoI) {}; // Default Constructor

	Link(Eigen::VectorXd& theta_): theta(theta_), accept(0), Q(numQoI) {};

	int getSDIM(){return STOCHASTIC_DIM;}
	VectorXd getTheta() const{return theta;}

  double getTheta(int i) const{return theta(i); }
	int getAccept(){	return accept;}
	VectorXd getQ(){ 	return Q;}

	double getlogPhi(bool isCoarse = false){ if(isCoarse){return logPhi_Coarse;} else {return logPhi;} }

	double getlogPi0(){ return logPi0; }

	void setlogPi0(double val){ logPi0 = val; }

	void setlogPhi_Coarse(double val) { logPhi_Coarse = val;}

	double getlogPhi_Coarse() { return logPhi_Coarse; }


	void operator=(Link &u){
		// Operator Overrides class with u
		VectorXd old = u.getTheta();
		(*this).setTheta(old);
		(*this).setlogPhi(u.getlogPhi());
		//VectorXd tmpQ = u.getQ();
		//(*this).setQoI(tmpQ);
	}

	void setTheta(VectorXd & vals){
		assert(vals.size() == STOCHASTIC_DIM);
		theta = vals;
	}


	void setQoI(VectorXd & vals){
		assert(vals.size() == numQoI);
		Q = vals;
	}

	void setAccepted(int val){accept = val;} // If accepted

	void setlogPhi(double val,bool coarseModel = false){
		if(coarseModel){
			logPhi_Coarse = val;
		}
		else{
			logPhi = val;
		}
		}

    int size(){return theta.size(); }


private:

	double logPhi, logPi0;

	double logPhi_Coarse;

	VectorXd theta;
	VectorXd Q;

};


  }
}


#endif
