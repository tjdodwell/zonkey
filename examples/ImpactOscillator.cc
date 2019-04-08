



#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
//#include <windows.h>
#include <string>

#include <Eigen/Dense>

#include "../MCMC/Chain/Link.hh"
#include "../MCMC/Chain/SingleChain.hh"
#include "../MCMC/Chain/MultiChain.hh"

#include "../MCMC/Proposals/SeqRandomWalk.hh"
#include "../MCMC/Proposals/SeqPCN.hh"

#include "../Model/ImpactOscillator.hh"

#include "../MCMC/multilevelMH.hh"

using namespace std;


int main()
{

  // getParameters
  double a = 1.0;
  double b = 10.0;
  double mu = 0.0;
  double sig = 1.0;

  int burninSamples = 100;
  int burningFactor = 10;

  int Nsamples = 10000;

  Eigen::VectorXd randomWalk_parameters(1);
  randomWalk_parameters(0) = 0.1;

  Eigen::VectorXd PCN_parameters(2);
  PCN_parameters(0) = 0.3;
  PCN_parameters(1) = 1.0;

  const int STOCHASTIC_DIM = 6;

  typedef Zonkey::MCMC::Link<STOCHASTIC_DIM,0> LINK;

  typedef Zonkey::MCMC::SingleChain<LINK> CHAIN;

  typedef Zonkey::MCMC::MultiChain<LINK,CHAIN> MULTICHAIN;

  typedef Zonkey::MCMC::SeqPCN<LINK> PROPOSAL;
    PROPOSAL myProposal(PCN_parameters);

  typedef Zonkey::Models::ImpactOscillator<LINK,STOCHASTIC_DIM> MODEL;
    MODEL F;




  return 0;
}
