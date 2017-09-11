#pragma once
#include "DataType.h"

class Quadrature
{
protected:
	Quadrature();
	~Quadrature();

public:
	// Gauss-Legendre guadrature x and weight / p.m. index
	static real_t gauss2X(int_t);
	static real_t gauss3X(int_t);
	static real_t gauss4X(int_t);
	static real_t gauss5X(int_t);

	static real_t lobatto2X(int_t);
	static real_t lobatto3X(int_t);
	static real_t lobatto4X(int_t);
	static real_t lobatto5X(int_t);

	static real_t gauss2W(int_t);
	static real_t gauss3W(int_t);
	static real_t gauss4W(int_t);
	static real_t gauss5W(int_t);

	// Legendre polynomial / p.m. order(0~), x coordinate
	static real_t legendreP(int_t, real_t);

private:
	// Quadrature variable
	static Quadrature _quad;
};

// Gauss-Legendre quadrature macro
// index : 0~
#define Gauss2_X(index) Quadrature::gauss2X(index)
#define Gauss2_W(index) Quadrature::gauss2W(index)

#define Gauss3_X(index) Quadrature::gauss3X(index)
#define Gauss3_W(index) Quadrature::gauss3W(index)

#define Gauss4_X(index) Quadrature::gauss4X(index)
#define Gauss4_W(index) Quadrature::gauss4W(index)

#define Gauss5_X(index) Quadrature::gauss5X(index)
#define Gauss5_W(index) Quadrature::gauss5W(index)

// Gauss-Lobatto points macro
// index : 0~
#define Lobatto2_X(index) Quadrature::lobatto2X(index)
#define Lobatto3_X(index) Quadrature::lobatto3X(index)
#define Lobatto4_X(index) Quadrature::lobatto4X(index)
#define Lobatto5_X(index) Quadrature::lobatto5X(index)

// Legendre polynomial macro
// index : 0~
#define LegendreP(order, x) Quadrature::legendreP(order, x)