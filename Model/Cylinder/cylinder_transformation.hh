


// GeometryGrid for a cylinder.

// Make a Grid which is 



double radius = 40.0;
double thickness = 4.0; // mm
double Length = 120.0; // mm


typedef Dune::YaspGrid<dim> Grid;
typedef Grid::ctype DF;

Dune::FieldVector<DF,dim> L;

L[0] = 2.0 * M_PI;
L[1] = thickness;
L[2] = Length;

Dune::FieldVector<double,3> inline gridTransformation(const FieldVec & x) const{
	// Grid Transformation for Cylinder
	y = Dune::FieldVector<double,3> y(0.0);
	y[0] = (inner_radius + x[1]) * std::cos(x[0]);
	y[1] = (inner_radius + x[1]) * std::sin(x[0]);
	y[2] = x[2];
	return y; 
}

        

        Grid grid(L,N);
        
        grid.globalRefine(maxL);