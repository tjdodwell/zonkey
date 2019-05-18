#ifndef GRIDTRANSFORMATION__HH
#define GRIDTRANSFORMATION__HH

    //! Grid Transformation
    /*!
     * Gridtransformation is a user-defined function which
     * gives the transformations of points x in the base grid (generated using LayerCake()) to the final
     * structural geometry
     * \f$ y = F(x) : \Omega -> \hat \Omega. \f$
     * This function not only allows us to define more complex
     * geometry, but grade the mesh towards areas of interest and introduce geometric defects.
     */
    template <int dim>
    class GridTransformation
    : public Dune :: AnalyticalCoordFunction< double, dim, dim, GridTransformation <dim> >{
      typedef GridTransformation This;
      typedef Dune :: AnalyticalCoordFunction< double, dim, dim, This > Base;

    public:
      typedef typename Base :: DomainVector DomainVector;
      typedef typename Base :: RangeVector RangeVector;

      GridTransformation(double inner_radius_): inner_radius(inner_radius_){  }

      void evaluate ( const DomainVector &x, RangeVector &y ) const {
          // Grid Transformation for Cylinder

          y[0] = (inner_radius + x[1]) * std::cos(x[0]);
          y[1] = (inner_radius + x[1]) * std::sin(x[0]);
          y[2] = x[2]; 
      }

    private:

      double inner_radius;
      
    };

#endif
