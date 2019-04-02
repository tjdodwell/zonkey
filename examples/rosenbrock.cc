



#include <iostream>
#include <Eigen/Dense>

#include "../MCMC/Chain/Link.hh"
#include "../MCMC/Chain/SingleChain.hh"

#include "../MCMC/Proposals/SeqRandomWalk.hh"

#include "../Model/Rosenbrock.hh"

#include "../MCMC/MH.hh"


int main()
{

  // getParameters
  double a = 1.0;
  double b = 100.0;
  double mu = 0.0;
  double sig = 1.0;

  Eigen::VectorXd randomWalk_parameters(1);
  randomWalk_parameters(0) = 0.3;

  std::cout << "== Starting Rosenbrock Example ==" << std::endl;

  const int STOCHASTIC_DIM = 2;

  typedef Zonkey::MCMC::Link<STOCHASTIC_DIM,0> LINK;

  typedef Zonkey::MCMC::SingleChain<LINK> CHAIN;
    CHAIN markovChain;

  typedef Zonkey::MCMC::SeqRandomWalk<LINK> PROPOSAL;
    PROPOSAL myProposal(randomWalk_parameters);

  typedef Zonkey::Models::Rosenbrock<LINK,STOCHASTIC_DIM> MODEL;
    MODEL F(a,b,mu,sig);

  Zonkey::MCMC::MetropolisHastings<LINK,CHAIN,PROPOSAL,MODEL> myMCMC(F,myProposal,markovChain);

  myMCMC.burnin(1000,10);


  return 0;
}
