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
      enum { doAlphaVolume = false };
      enum { doAlphaBoundary = true };

      QoI (T& param_, int intorderadd_=0)
        : param(param_), intorderadd(intorderadd_)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {

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

        // skip rest if we are on Dirichlet boundary
        if (bctype==ConvectionDiffusionBoundaryConditions::Dirichlet) return;

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


        if(coords[0] > 1.0 - 1e-6){



                  for (const auto& ip : quadratureRule(geo,intorder))
                    {
                      // position of quadrature point in local coordinates of element
                      auto local = geo_in_inside.global(ip.position());

                      auto tensor = param.A(cell_inside,local);

                      // evaluate shape functions (assume Galerkin method)
                      auto& phi = cache.evaluateFunction(local,lfsu_s.finiteElement().localBasis());

                      // evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
                      auto& js = cache.evaluateJacobian(local,lfsu_s.finiteElement().localBasis());

                      // transform gradients of shape functions to real element
                      jac = geo.jacobianInverseTransposed(ip.position());
                      for (size_type i=0; i<lfsu_s.size(); i++)
                        jac.mv(js[i][0],gradphi[i]);

                      // compute gradient of u
                      gradp = 0.0;
                      for (size_type i=0; i<lfsu_s.size(); i++)
                        gradp.axpy(x_s(lfsu_s,i),gradphi[i]);

                      // compute A * gradient of p
                      tensor.mv(gradp,Agradp);


                      auto factor = ip.weight()*geo.integrationElement(ip.position());

                      for (size_type i=0; i<lfsu_s.size(); i++)
                        r_s.accumulate(lfsu_s,i,-Agradp[0]*phi[i]*factor);

                    }

        } //
      }

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
