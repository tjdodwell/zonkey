#ifndef ZONKEY_MCMC_PROPOSAL_HH
#define ZONKEY_MCMC_PROPOSAL_HH

#include "../Chain/Link.hh"


namespace Zonkey {
  namespace PDELab {


    class Proposal{

      typename typedef Zonkey::MCMC::Link LINK;

      public:

        Proposal(){}

        Zonkey::MCMC::link make(const LINK& currentState){

          LINK prop(currentState);

        }


      private:


    }


  }
}


#endif
