template<typename DF, typename RF, unsigned int dimension>
class GridTraits
{
  public:

    enum {dim = dimension};

    using RangeField  = RF;
    using Scalar      = Dune::FieldVector<RF,1>;
    using DomainField = DF;
    using Domain      = Dune::FieldVector<DF,dim>;
};
