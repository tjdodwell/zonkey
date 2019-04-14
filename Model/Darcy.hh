
/*


const int dim=2;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<DF,dim> L;
        L[0] = ptree.get("grid.structured.LX",(double)1.0);
        L[1] = ptree.get("grid.structured.LY",(double)1.0);
        std::array<int,dim> N;
        N[0] = ptree.get("grid.structured.NX",(int)10);
        N[1] = ptree.get("grid.structured.NY",(int)10);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N));
        gridp->globalRefine(refinement);
        typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();
        if (degree==1) {
          typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }


*/

#include <iostream>
#include <fstream>

#include "dune_includes/problem.hh"
#include "dune_includes/randomField/RandomField.hh"
using namespace std;


template <typename GRID>
class Darcy{

  public:

    Darcy(){


      int nobs = 5;
      int Nobs = nobs * nobs;

      int maxR = 169;

      MatrixXd

      obsCoord.resize(Nobs,2);

      VectorXd

      Fobs.resize(Nobs);

      // Construct Coordinates locations of the data

      VectorXd coords(nobs);

      double dx = 1.0 / (nobs + 1);

      coords(0) = dxcoords;
      for (int i = 1; i < nobs; i++){
      	coords(i) = coords(i-1) + dxcoords;
      }

      int kk = 0;
      for (int i = 0; i < nobs; i++){
      	for (int j = 0; j < nobs; j++){
      		obsCoord(kk,0) = coords(i);
      		obsCoord(kk,1) = coords(j);
      		kk += 1;
      	}
      }

      // == Read data in from file
      ifstream input("data.txt");
      thetaObs.resize(maxR);
      for (int i = 0; i < maxR; i++) {
        input >> thetaObs(i);
      }

      // == Read in obs

      Fobs.resize(Nobs);

      Fobs(0) = 0.103861;
      Fobs(1) = 0.133822;
      Fobs(2) = 0.162300;
      Fobs(3) = 0.163523;
      Fobs(4) = 0.127226;

      Fobs(5) = 0.233255;
      Fobs(6) = 0.315055;
      Fobs(7) = 0.391402;
      Fobs(8) = 0.380278;
      Fobs(9) = 0.320238;

      Fobs(10) = 0.466618;
      Fobs(11) = 0.515041;
      Fobs(12) = 0.544738;
      Fobs(13) = 0.569354;
      Fobs(14) = 0.562334;

      Fobs(15) = 0.602408;
      Fobs(16) = 0.630779;
      Fobs(17) = 0.650833;
      Fobs(18) = 0.678896;
      Fobs(19) = 0.701951;

      Fobs(20) = 0.772371;
      Fobs(21) = 0.783087;
      Fobs(22) = 0.802645;
      Fobs(23) = 0.812591;
      Fobs(24) = 0.823080;

    }


    void apply(Link & u, int level = 0, bool plotSolution = false){

      // Setting up Model

      typedef GRID::LevelGridView GV;
      GV gv = grid.levelGridView(level);

      // dimension and important types
      const int dim = GV::dimension;
      typedef double RF;                   // type for computations


      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
        FEM fem(gv);





      typedef GenericEllipticProblem<GV,RF,RandomField> PROBLEM;
        PROBLEM problem(gv,field);

      // Make grid function space
      typedef Dune::PDELab::ConformingDirichletConstraints CON;
      typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
      GFS gfs(gv,fem);
      gfs.name("p");

      // Assemble constraints
      typedef typename GFS::template
        ConstraintsContainer<RF>::Type CC;
      CC cc;
      Dune::PDELab::constraints(b,gfs,cc); // assemble constraints
      std::cout << "constrained dofs=" << cc.size() << " of "
                << gfs.globalSize() << std::endl;

      // A coefficient vector
      using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
      Z z(gfs); // initial value

      // Make a grid function out of it
      typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
      ZDGF zdgf(gfs,z);

      // Fill the coefficient vector
      Dune::PDELab::interpolate(g,gfs,z);

      // Make a local operator
      typedef NonlinearPoissonFEM<Problem<RF>,FEM> LOP;
      LOP lop(problem);

      // Make a global operator
      typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
      int degree = ptree.get("fem.degree",(int)1);
      MBE mbe((int)pow(1+2*degree,dim));
      typedef Dune::PDELab::GridOperator<
        GFS,GFS,  /* ansatz and test space */
        LOP,      /* local operator */
        MBE,      /* matrix backend */
        RF,RF,RF, /* domain, range, jacobian field type*/
        CC,CC     /* constraints for ansatz and test space */
        > GO;
      GO go(gfs,cc,gfs,cc,lop,mbe);

      // Select a linear solver backend
      typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
      LS ls(100,2);






      u.setlogPhi(logLikelihood);




    }


  private:

    Eigen::VectorXd thetaObs, Fobs;


    RandomField field;








}
