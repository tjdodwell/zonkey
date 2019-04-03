



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
  double b = 10.0;
  double mu = 0.0;
  double sig = 1.0;

  int burninSamples = 100;
  int burningFactor = 10;

  int Nsamples = 10000;

  Eigen::VectorXd randomWalk_parameters(1);
  randomWalk_parameters(0) = 0.1;

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

  myMCMC.burnin(burninSamples,burningFactor);

  myMCMC.run(Nsamples);

  auto theChain = myMCMC.getChain();

  std::cout << "Effective Sample size / Samples = " << theChain.getMaxESS() << " / " << theChain.size() << std::endl;


  return 0;
}
