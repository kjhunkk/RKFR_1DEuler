#pragma once
#include "DataType.h"
#include "Quadrature.h"
#include "Grid.h"

class FRbasis
{
public:
	// Constructor
	FRbasis(std::shared_ptr<Grid>);

	// Destructor
	~FRbasis();

public:
	// Functions
	// Return solution points(Gauss-Lobatto)
	// p.m. point number, polynomial order
	// r.t. solution points in [-1,1]
	real_t getSolutionPoint(int_t, int_t);

	// N-th order Lagrange polynomial basis in Gauss-Lobatto points
	// p.m. x coordinate, basis number, polynomial order
	real_t LagrangeP(real_t, int_t, int_t);

	// N-th order Lagrange polynomial basis derivative in Gauss-Lobatto points
	// p.m. x coordinate, basis number, polynomial order
	real_t diff_LagrangeP(real_t, int_t, int_t);

	// Convert from physical coordinate to [-1,1]
	// p.m. cell number, physical coordinate
	real_t toCompCoord(int_t, real_t);

protected:
	std::shared_ptr<Grid> _grid;
	real_t _sizeX;
};