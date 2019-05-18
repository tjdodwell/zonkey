


#include <iostream>
#include <fstream>

#include "dune_includes/problem.hh"
#include "dune_includes/problemQoI.hh"
#include "dune_includes/randomField/RandomField.hh"
#include "dune_includes/QoI.hh"




#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace Zonkey {
  namespace Models {


template <typename Link, int STOCHASTIC_DIM, typename GRID>
class Darcy{

  public:

    Darcy(GRID& grid_): grid(grid_){

      // == Parameters from the Section 4. Dodwell et al. 2015

      sigf = 0.031622776601684;

      // == Setup Random Field

      double L = 1.0;
      double sigKl = 1.0;
      double correlation_length = 0.5;
      int maxR = 169;

      field.setup(L,sigKl,correlation_length,maxR);


      int nobs = 5;
      Nobs = nobs * nobs;



      obsCoord.resize(Nobs,2);

      Fobs.resize(Nobs);

      // Construct Coordinates locations of the data

      VectorXd coords(nobs);

      double dx = 1.0 / (nobs + 1);
      coords(0) = dx;
      for (int i = 1; i < nobs; i++){
      	coords(i) = coords(i-1) + dx;
      }

      int kk = 0;
      for (int i = 0; i < nobs; i++){
      	for (int j = 0; j < nobs; j++){
      		obsCoord(kk,0) = coords(i);
      		obsCoord(kk,1) = coords(j);
      		kk += 1;
      	}
      }

    /*  // == Read data in from file
      ifstream input("data.txt");
      thetaObs.resize(maxR);
      for (int i = 0; i < maxR; i++) {
        input >> thetaObs(i);
      }*/

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

    void inline plot_field(const Eigen::VectorXd& theta,  std::string title = "randomfield", int level = -1){

      std::vector<double> xi_vec(theta.size());

      for (int i = 0; i < theta.size(); i++){
        xi_vec[i] = theta(i);
      }

      field.setXi(xi_vec);

      if(level < 0){
        level = grid.maxLevel();
      }

    typedef typename GRID::LevelGridView GV;
        GV gv = grid.levelGridView(level); // Get finest grid

    typedef typename GV::Grid::ctype Coord;

    typedef typename GV::Grid::ctype e_ctype;
    typedef Dune::PDELab::QkDGLocalFiniteElementMap<e_ctype,double,0,2> Q0;

    std::cout << "My Grid Size " << gv.size(0) << std::endl;

    //typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,double,0> Q0;
        Q0 fem;
    typedef Dune::PDELab::ConformingDirichletConstraints CON;
    typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
    typedef Dune::PDELab::GridFunctionSpace<GV, Q0, CON, VBE> GFS;
        GFS gfs(gv,fem); gfs.name("random_field");

    typedef typename Dune::PDELab::Backend::impl::BackendVectorSelector<GFS,double>::Type U;
        U k(gfs,0.0);

    using Dune::PDELab::Backend::native;

    for (const auto& eit : elements(gv)){
      int id = gv.indexSet().index(eit);
      native(k)[id] = field.getPerm(eit.geometry().center());
    }

    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
        DGF xdgf(gfs,k);

    // Write solution to VTK
          Dune::VTKWriter<GV> vtkwriter(gv);
          typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
          auto adapt = std::make_shared<ADAPT>(xdgf,title);
          vtkwriter.addVertexData(adapt);
          vtkwriter.write(title);
    
}


    Eigen::VectorXd samplePrior(int level = 0){
      std::random_device rd;
      std::normal_distribution<double> dis(0.0,1.0);
      std::mt19937 gen(rd());
      Eigen::VectorXd z(STOCHASTIC_DIM);
      for (int i = 0; i < STOCHASTIC_DIM; i++){
        z(i) = dis(gen);
      }
      return z; // Sample from prior - note Sigma = C' * C
    } // samplePrior


    void apply(Link& u, int level = 0, bool plotSolution = false, bool setasData = false){

      // Unwrap Xi

      Eigen::VectorXd xi = u.getTheta();

      std::vector<double> xi_vec(xi.size());

      for (int i = 0; i < xi.size(); i++){
        xi_vec[i] = xi(i);
      }

      field.setXi(xi_vec);


      // Setting up Model
      typedef typename GRID::LevelGridView GV;
      GV gv = grid.levelGridView(level);

      // dimension and important types
      const int dim = GV::dimension;
      typedef double RF;                   // type for computations
      typedef typename GRID::ctype DF;


      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
        FEM fem(gv);


      typedef GenericEllipticProblem<GV,RF,RandomField<2> > PROBLEM;
        PROBLEM problem(gv,field);

      typedef GenericEllipticProblemQoI<GV,RF,RandomField<2>> PROB_QOI;
        PROB_QOI qp(gv,field);

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PROBLEM> BCType;
        BCType bctype(gv,problem);

        typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PROB_QOI> BCType_qoi;
          BCType_qoi bctype_qoi(gv,qp);

      // Make grid function space
      typedef Dune::PDELab::ConformingDirichletConstraints CON;
      typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
      GFS gfs(gv,fem);
      gfs.name("p");

      typedef Dune::PDELab::Backend::Vector<GFS,RF> V;
        V x(gfs,0.0);

      // Extract domain boundary constraints from problem definition, apply trace to solution vector
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PROBLEM> G;
        G g(gv,problem);
      Dune::PDELab::interpolate(g,gfs,x);

      // Assemble constraints
      typedef typename GFS::template
        ConstraintsContainer<RF>::Type CC;
      CC cc;
      Dune::PDELab::constraints(bctype,gfs,cc); // assemble constraints

      CC cc_qoi;
      Dune::PDELab::constraints(bctype_qoi,gfs,cc_qoi); // assemble constraints

      // LocalOperator for given problem
      typedef Dune::PDELab::ConvectionDiffusionFEM<PROBLEM,FEM> LOP;
      LOP lop(problem);

      // Make a global operator
      typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
      MBE mbe(9);

      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC> GO;
      GO go(gfs,cc,gfs,cc,lop,mbe);

      // Select a linear solver backend
      typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
      LS ls(100,0);

      // Assemble and solve linear problem
      typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
      SLP slp(go,ls,x,1e-10,1e-99,0);
      slp.apply(); // here all the work is done!

      typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
        DGF xdgf(gfs,x);

      Eigen::VectorXd F(Nobs);

      for (int i = 0; i < Nobs; i++){ // For each observation

        Dune::FieldVector<double,2> point(0.0);
        point[0] = obsCoord(i,0);
        point[1] = obsCoord(i,1);

        Dune::PDELab::GridFunctionProbe<DGF> probe(xdgf,point);

        Dune::FieldVector<double,1> val(0);

        probe.eval(val);

        F(i) = val[0];

        if (setasData){
          Fobs(i) = F(i);
        }

      }

      if(setasData){std::cout << Fobs << std::endl;}


      /*typedef Dune::PDELab::QoI<PROBLEM,FEM> QOI;
        QOI qoi_lop(problem);

      typedef Dune::PDELab::GridOperator<GFS,GFS,QOI,MBE,RF,RF,RF,CC,CC> QGO;
        QGO qgo(gfs,cc,gfs,cc,qoi_lop,mbe);


      using Dune::PDELab::Backend::native;

      typedef Dune::PDELab::Backend::Vector<GFS,RF> V;
        V q(gfs,0.0);

      qgo.residual(x,q);

      Eigen::VectorXd Q(1);
      Q(0) = 0.0;
      for (int i = 0; i < native(q).size(); i++){
        Q(0) += native(q)[i];
      }

      std::cout << " " << std::endl;
      std::cout << Q << std::endl;*/
      Eigen::VectorXd Q(1);
      Q(0) = xi(0);
      u.setQoI(Q);

      // Compute logLikelihood
      Eigen::VectorXd misMatch(Nobs);

      double logLikelihood = 0.0;
      for (int k = 0; k < Nobs; k++){
        logLikelihood -= (F(k) - Fobs(k)) * (F(k) - Fobs(k))  / (sigf * sigf);
      }

      u.setlogPhi(logLikelihood);


      if(plotSolution){
        std::cout << "logDensity = " << logLikelihood << std::endl;
        // Write solution to VTK
          Dune::VTKWriter<GV> vtkwriter(gfs.gridView());
          typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
          auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
          vtkwriter.addVertexData(adapt);
          vtkwriter.write("solution");
      }



    }


  private:

    Eigen::VectorXd thetaObs, Fobs;
    Eigen::MatrixXd obsCoord;

    int Nobs;

    double sigf;

    GRID& grid;

    RandomField<2> field;

};

}
}
