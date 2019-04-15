


#include <iostream>
#include <fstream>

#include "dune_includes/problem.hh"
#include "dune_includes/randomField/RandomField.hh"

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


template <typename Link, int STOCHASTIC_DIM, typename GRID>
class Darcy{

  public:

    Darcy(GRID& grid_): grid(grid_){


      int nobs = 5;
      int Nobs = nobs * nobs;

      int maxR = 169;

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


    void apply(Link & u, int level = 0, bool plotSolution = false){

      // Setting up Model
      typedef typename GRID::LevelGridView GV;
      GV gv = grid.levelGridView(level);

      // dimension and important types
      const int dim = GV::dimension;
      typedef double RF;                   // type for computations
      typedef typename GRID::ctype DF;


      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
        FEM fem(gv);


      typedef GenericEllipticProblem<GV,RF,RandomField> PROBLEM;
        PROBLEM problem(gv,field);

      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PROBLEM> BCType;
        BCType bctype(gv,problem);

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
      LS ls(100,2);

      // Assemble and solve linear problem
      typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
      SLP slp(go,ls,x,1e-10);
      slp.apply(); // here all the work is done!


      // Need to loop over elements

      std::vector<int> elem(Nobs);

      std::vector<Dune::FieldVector<double,2>> localCoords(Nobs);




      for (int i = 0; i < Nobs; i++){
          int id_min = 0;
          double dmin = 10.0;
          for (const auto& eit : elements(gv)){
                int id = gv.indexSet().index(eit);
                auto point = eit.geometry().center();
                double d = std::sqrt(std::pow(obsCoord(i,0)-point[0],2)+std::pow(obsCoord(i,1)-point[1],2));
                if (d < dmin){
                  d = dmin; id_min = id;
                  Dune::FieldVector<double,2> obs(0.0);
                  obs[0] = obsCoord(i,0); obs[1] = obsCoord(i,1);
                  localCoords[i] = eit.geometry.local(obs);
                }
          }
          elem[i] = id_min;
        }

        // make local function space
        typedef Dune::PDELab::LocalFunctionSpace<GFS> CLFS;
          CLFS lfsu(gfs);
        typedef Dune::PDELab::LFSIndexCache<CLFS> CLFSCache;
          CLFSCache clfsCache(lfsu);
        std::vector<double> u(lfsu.maxSize());

        typedef typename X::template ConstLocalView<CLFSCache> XView;
        XView xView(x);

        typedef typename CLFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianType;



              for (const auto& eg : elements(gv)){

                // bind solution x to local element
                lfsu.bind(eg);
                clfsCache.update();
                xView.bind(clfsCache);
                xView.read(u);



                // select quadrature rule
                auto geo = eg.geometry();
                const Dune::QuadratureRule<double,dim>& rule = Dune::QuadratureRules<double,dim>::rule(geo.type(),1);

                // loop over quadrature points
                for (const auto& ip : rule)
                {
                  // Evaluate Jacobian
                    std::vector<JacobianType> js(lfsu.size());
                    lfsu.finiteElement().localBasis().evaluateJacobian(ip.position(),js);

                    // Transform gradient to real element
                    auto jac = geo.jacobianInverseTransposed(ip.position());
                    std::vector<Dune::FieldVector<double,dim> > gradphi(lfsu.size());

                    for (int i=0; i < lfsu.size(); i++){
                        gradphi[i] = 0.0;
                        jac.umv(js[i][0],gradphi[i]);
                    }

                    Dune::FieldMatrix<double,dim,3> G(0.0);

                    for (int i = 0; i < lfsu.size(); i++){
                      for (int j = 0; j < dim; j++){
                        G[j][i] = gradphi[i][j];
                      }
                    }
                    std::vector<double> gradp(dim), flux(dim);

                    G.mv(u,gradp); // compute pressure gradient  grad(p) = gradphi * p
                    Kij.mv(gradp,flux); // Compute flux = - Perm * G

                    for (int d=0; d < dim; d++){
                        flux_all[is.index(eg)][d] = flux[d];
                    }

                } // end For each integration point

              } // end For each element







      //u.setlogPhi(logLikelihood);



    }


  private:

    Eigen::VectorXd thetaObs, Fobs;
    Eigen::MatrixXd obsCoord;

    GRID& grid;

    RandomField field;

};
