// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include <dune/common/exceptions.hh> // We use exceptions
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>

// dune-istl included by pdelab
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>

#include <dune/pdelab/common/functionutilities.hh>

#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/newton/newton.hh>


#include <dune/typetree/treepath.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include "/Users/td336/zonkey/Model/Darcy.hh"
#include "/Users/td336/zonkey/MCMC/Chain/Link.hh"
#include "/Users/td336/zonkey/MCMC/Chain/SingleChain.hh"
#include "/Users/td336/zonkey/MCMC/Proposals/SeqPCN.hh"

#include "/Users/td336/zonkey/MCMC/MH.hh"

#include "/Users/td336/zonkey/MCMC/Proposals/SeqDA.hh"

int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

        int Nsamples = 1000;

        const int M0 = 8;

        const int maxL = 4;

        const int dim = 2;

        const int maxR = 169;

        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;

        Dune::FieldVector<DF,dim> L;

        L[0] = 1.0;
        L[1] = 1.0;

        std::array<int,dim> N;
        N[0] = M0 + 1;
        N[1] = M0 + 1;

        Grid grid(L,N);

        grid.globalRefine(maxL);

        Eigen::VectorXd PCN_parameters(2);
          PCN_parameters(0) = 0.333;
          PCN_parameters(1) = 0.100;

        const int STOCHASTIC_DIM = 169;

        typedef Zonkey::MCMC::Link<STOCHASTIC_DIM,1> LINK;

        typedef Zonkey::MCMC::SingleChain<LINK> CHAIN;
          CHAIN markovChain;

        typedef Zonkey::MCMC::SeqPCN<LINK> PROPOSAL;
          PROPOSAL myProposal(PCN_parameters);

        typedef Zonkey::Models::Darcy<LINK,STOCHASTIC_DIM,Grid> MODEL;
          MODEL F(grid);

        typedef Zonkey::MCMC::SeqDA<LINK,CHAIN,PROPOSAL,MODEL> DA;
          DA seqDA(myProposal,F);

       // Zonkey::MCMC::MetropolisHastings<LINK,CHAIN,PROPOSAL,MODEL> myMCMC(F,myProposal,markovChain);

        Zonkey::MCMC::MetropolisHastings<LINK,CHAIN,DA,MODEL> myMCMC(F,seqDA,markovChain);

        // Test the model

        /*Eigen::VectorXd xi = F.samplePrior();

        std::cout << xi << std::endl;

        LINK testSample(xi);

        int level = 0;

        F.apply(testSample,level,true);*/

        myMCMC.burnin(1000,1);

        myMCMC.run(Nsamples,1);

        auto theChain = myMCMC.getChain();

        std::cout << "Effective Sample size / Samples = " << theChain.getMaxESS() << " / " << theChain.size() << std::endl;

        std::cout << "Acceptance Ratio = " << theChain.acceptRatio() << std::endl;


    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
