// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONFEM_QOI_HH
#define DUNE_PDELAB_LOCALOPERATOR_CONVECTIONDIFFUSIONFEM_QOI_HH

#include<vector>
#include<type_traits>

#include<dune/common/deprecated.hh>

#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

namespace Dune {
  namespace PDELab {

    /** a local operator for solving the linear convection-diffusion equation with standard FEM
     *
     * \f{align*}{
     *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
     *                                              u &=& g \mbox{ on } \partial\Omega_D \\
     *                (b(x) u - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
     *                        -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O
     * \f}
     * Note:
     *  - This formulation is valid for velocity fields which are non-divergence free.
     *  - Outflow boundary conditions should only be set on the outflow boundary
     *
     * \tparam T model of ConvectionDiffusionParameterInterface
     */
    template<typename T, typename FiniteElementMap>
    class QoI :
      public Dune::PDELab::NumericalJacobianApplyVolume<QoI<T,FiniteElementMap> >,
      public Dune::PDELab::NumericalJacobianApplyBoundary<QoI<T,FiniteElementMap> >,
      public Dune::PDELab::NumericalJacobianVolume<QoI<T,FiniteElementMap> >,
      public Dune::PDELab::NumericalJacobianBoundary<QoI<T,FiniteElementMap> >,
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename T::Traits::RangeFieldType>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = false};
      enum { doAlphaBoundary = true };

      QoI (T& param_, int intorderadd_=0)
        : param(param_), intorderadd(intorderadd_)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        /*// Define types
        using RF = typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSU::Traits::SizeType;

        // dimensions
        const int dim = EG::Entity::dimension;

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometry
        auto geo = eg.geometry();

        // evaluate diffusion tensor at cell center, assume it is constant over elements
        auto ref_el = referenceElement(geo);
        auto localcenter = ref_el.position(0,0);
        auto tensor = param.A(cell,localcenter);


        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        Dune::FieldVector<RF,dim> gradu(0.0);
        Dune::FieldVector<RF,dim> Agradu(0.0);

        // Transformation matrix
        typename EG::Geometry::JacobianInverseTransposed jac;

        // loop over quadrature points
        auto intorder = intorderadd+2*lfsu.finiteElement().localBasis().order();
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // update all variables dependent on A if A is not cell-wise constant
            if (!Impl::permeabilityIsConstantPerCell<T>(param))
            {
              tensor = param.A(cell, ip.position());
            }

            // evaluate basis functions
            auto& phi = cache.evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());

            // evaluate u
            RF u=0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              u += x(lfsu,i)*phi[i];

            // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
            auto& js = cache.evaluateJacobian(ip.position(),lfsu.finiteElement().localBasis());

            // transform gradients of shape functions to real element
            jac = geo.jacobianInverseTransposed(ip.position());
            for (size_type i=0; i<lfsu.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            // compute gradient of u
            gradu = 0.0;
            for (size_type i=0; i<lfsu.size(); i++)
              gradu.axpy(x(lfsu,i),gradphi[i]);


            // compute A * gradient of u
            tensor.mv(gradu,Agradu);

            // evaluate velocity field, sink term and source term
            auto b = param.b(cell,ip.position());
            auto c = param.c(cell,ip.position());
            auto f = param.f(cell,ip.position());

            // integrate (A grad u)*grad phi_i - u b*grad phi_i + c*u*phi_i
            RF factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_type i=0; i<lfsu.size(); i++)
              r.accumulate(lfsu,i,( Agradu[0]*gradphi[i][0])*factor);
          }*/

      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {

      }

      // boundary integral
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {

        // dimensions
        const int dim = IG::Entity::dimension;

        // Define types
        using RF = typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType;
        using size_type = typename LFSV::Traits::SizeType;

        // Reference to the inside cell
        const auto& cell_inside = ig.inside();

        // Get geometry
        auto geo = ig.geometry();

        // Get geometry of intersection in local coordinates of cell_inside
        auto geo_in_inside = ig.geometryInInside();

        // evaluate boundary condition type
        auto ref_el = referenceElement(geo_in_inside);
        auto local_face_center = ref_el.position(0,0);
        auto intersection = ig.intersection();
        auto bctype = param.bctype(intersection,local_face_center);



        // loop over quadrature points and integrate normal flux
        auto intorder = intorderadd+2*lfsu_s.finiteElement().localBasis().order();

        //

        // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu_s.size());
        Dune::FieldVector<RF,dim> gradp(0.0);
        Dune::FieldVector<RF,dim> Agradp(0.0);


        // Transformation matrix
        typename IG::Geometry::JacobianInverseTransposed jac;



        // Unwrap solution on edge

        auto coords = geo.center();


        if(coords[0] <  1e-6){

              Dune::FieldVector<double,2> local(0.5);

              auto k = param.A(cell_inside,local);

              double h = std::abs(geo.corner(0)[1] - geo.corner(1)[1]);

              auto flux = 0.5 * k[0][0] * (x_s(lfsu_s,1) + x_s(lfsu_s,3));

              /*std::cout << "The flux is = " << flux << std::endl;

              std::cout << x_s(lfsu_s,1) << std::endl;
              std::cout << x_s(lfsu_s,3) << std::endl; */

              
              // Share flux equally between two nodes

              r_s.accumulate(lfsu_s, 0, 0.5 * flux); 
              
              r_s.accumulate(lfsu_s, 2, 0.5 * flux);

        }  // end if on right hand boundary
      } // end alpha_boundary

      // jacobian contribution from boundary
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_s) const
      {

      }


      //! set time in parameter class
      void setTime (double t)
      {
        param.setTime(t);
      }

    private:
      T& param;
      int intorderadd;
      using LocalBasisType = typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType;
      Dune::PDELab::LocalBasisCache<LocalBasisType> cache;
    };


  } // end namespace PDELab
} // end namespace Dune
#endif
