



#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
//#include <windows.h>
#include <string>

#include <Eigen/Dense>

#include "../MCMC/Chain/Link.hh"
#include "../MCMC/Chain/SingleChain.hh"

#include "../MCMC/Proposals/SeqRandomWalk.hh"
#include "../MCMC/Proposals/SeqPCN.hh"

#include "../Model/Rosenbrock.hh"

#include "../MCMC/MH.hh"

using namespace std;

string printProg(int x){
    string s;
    s="[";
    for (int i=1;i<=(100/2);i++){
        if (i<=(x/2) || x==100)
            s+="=";
        else if (i==(x/2))
            s+=">";
        else
            s+=" ";
    }

    s+="]";
    return s;
}



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

  const int STOCHASTIC_DIM = 2;

  typedef Zonkey::MCMC::Link<STOCHASTIC_DIM,0> LINK;

  typedef Zonkey::MCMC::SingleChain<LINK> CHAIN;
    CHAIN markovChain;

//  typedef Zonkey::MCMC::SeqRandomWalk<LINK> PROPOSAL;
//    PROPOSAL myProposal(randomWalk_parameters);
  typedef Zonkey::MCMC::SeqPCN<LINK> PROPOSAL;

    PROPOSAL myProposal(PCN_parameters);

  typedef Zonkey::Models::Rosenbrock<LINK,STOCHASTIC_DIM> MODEL;
    MODEL F(a,b,mu,sig);

  Zonkey::MCMC::MetropolisHastings<LINK,CHAIN,PROPOSAL,MODEL> myMCMC(F,myProposal,markovChain);

  myMCMC.burnin(burninSamples,burningFactor);

  myMCMC.run(Nsamples);

  auto theChain = myMCMC.getChain();

  std::cout << "Effective Sample size / Samples = " << theChain.getMaxESS() << " / " << theChain.size() << std::endl;

  srand(time(0));
    int x = 0;
    cout << "Working ..." << endl;
    for(int i = 0; i<100; i++){
        int r = rand() % 1000;
        x++;
        cout << "\r" << setw(-20) << printProg(x) << " " << x << "% completed." << flush;
    }

    cout << endl << endl << "Operation completed successfully.\n" << flush;


  return 0;
}
