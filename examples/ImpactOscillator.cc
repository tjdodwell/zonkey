



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

#include "../MCMC/MH.hh"


using namespace std;


int main()
{

  int burninSamples = 1000;
  int burningFactor = 10;

  int Nsamples = 20000;

  Eigen::VectorXd randomWalk_parameters(1);
  randomWalk_parameters(0) = 0.1;

  Eigen::VectorXd PCN_parameters(2);
  PCN_parameters(0) = 0.3;
  PCN_parameters(1) = 0.2;

  const int STOCHASTIC_DIM = 5;

  typedef Zonkey::MCMC::Link<STOCHASTIC_DIM,0> LINK;

  typedef Zonkey::MCMC::SingleChain<LINK> CHAIN;
    CHAIN markovChain;

  typedef Zonkey::MCMC::MultiChain<LINK,CHAIN> MULTICHAIN;

  typedef Zonkey::MCMC::SeqPCN<LINK> PROPOSAL;
    PROPOSAL myProposal(PCN_parameters);

  typedef Zonkey::Models::ImpactOscillator<LINK,STOCHASTIC_DIM> MODEL;
    MODEL F;

  Zonkey::MCMC::MetropolisHastings<LINK,CHAIN,PROPOSAL,MODEL> myMCMC(F,myProposal,markovChain);

  myMCMC.burnin(burninSamples,burningFactor);

  myMCMC.run(Nsamples);

  auto theChain = myMCMC.getChain();

/*  Eigen::VectorXd mu(STOCHASTIC_DIM);

  mu(0) = 0.0; // initial displacement
  mu(1) = 0.1; // initial velocity
  mu(2) = 1.0; // stiffness
  mu(3) = 0.1; // Amplitude of osciallations of contact plate
  mu(4) = 1.0; // Frequency of contact plateAmplitude of osciallations of contact plate

  Eigen::VectorXd start = F.samplePrior();

  LINK testLink(start);

  F.apply(testLink,0,true);

  LINK proposal = myProposal.apply(testLink);

  F.apply(proposal,0,true);

  std::cout << myProposal.acceptReject(testLink,proposal,true) << std::endl;*/

  //F.apply(testLink,0,true);

  std::cout << "Acceptance Ratio " << theChain.acceptRatio() << std::endl;

  std::cout << "Effective Sample size / Samples = " << theChain.getMinESS() << " / " << theChain.size() << std::endl;


  return 0;
}
