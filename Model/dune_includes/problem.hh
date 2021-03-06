/*
 * Defining a Darcy problem with alternating layers of permeability and a high contrast
 */

template<typename GV, typename RF, typename RandomField>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  GenericEllipticProblem(GV& gv_, RandomField& field_)
  : gv(gv_),
    field(field_){


  }



  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    const typename Traits::DomainType xglobal = e.geometry().global(x);

    //field.evaluate(xglobal, coeff);

  //  std::cout << "This is rank = " << gv.comm().rank() << " try it here " << xglobal << std::endl;

    double coeff = field.getPerm(xglobal);

  //  std::cout << "Otherside  rank = " << gv.comm().rank() << std::endl;


    //std::cout << coeff << std::endl;


    //field.evaluate(xglobal, coeff);

    typename Traits::PermTensorType I;
    I[0][0] = coeff;
    I[0][1] = 0;
    I[1][0] = 0;
    I[1][1] = coeff;
    /*I[0][2] = 0;
    I[2][0] = 0;
    I[2][2] = coeff;
    I[1][2] = 0;
    I[2][1] = 0;*/

    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    /*typename Traits::DomainType xglobal = e.geometry().global(x);
    return - std::exp(-std::sqrt(std::pow(xglobal[0] - 0.25, 2) + std::pow(xglobal[1] - 0.25, 2)))
           + std::exp(-std::sqrt(std::pow(xglobal[0] - 0.75, 2) + std::pow(xglobal[1] - 0.75, 2)));
    */

    typename Traits::DomainType xglobal = e.geometry().global(x);

      return 1.0;

  }

  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);

    if(xglobal[0] < 1.0e-6 || xglobal[0] > 1.0 - 1e-6){
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    }
    else{
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    }


  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    /*if (xglobal[1] < 1E-10 && xglobal[0] > 0.75)
      return std::cos(xglobal[0] * 8* 3.141592);
    else if (xglobal[1] > 1.0-1E-10 && xglobal[0] < 0.25)
      return std::sin(xglobal[0] * 8 * 3.141592);
    else if (xglobal[0] < 1E-10)
      return std::sin(xglobal[1] * 2 * 3.141592) + 1.0;
    else if (xglobal[0] > 1.0 - 1E-10)
      return std::cos(xglobal[1] * 2 * 3.141592) - 1.0;
    else
      return 0.0;*/

    if(xglobal[0] > 1.0 - 1e-6){
      return 1.0;
    }
    else{
      return 0.0;
    }
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  void setTheta(std::vector<double>& xi){  field.setXi(xi); }

private:
  GV& gv;
  RandomField& field;

};
