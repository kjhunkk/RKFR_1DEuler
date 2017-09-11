#include "FRbasis.h"

FRbasis::FRbasis(std::shared_ptr<Grid> grid)
{
	_grid = grid;
	_sizeX = grid->getSizeX();
}

FRbasis::~FRbasis()
{

}

real_t FRbasis::getSolutionPoint(int_t n, int_t N)
{
	switch (N)
	{
	case 0: return 0.0;
	case 1: return Lobatto2_X(n);
	case 2: return Lobatto3_X(n);
	case 3: return Lobatto4_X(n);
	case 4: return Lobatto5_X(n);
	default:
		ERROR("exceeds maximum polynomial number");
		return 0.0;
	}
}

real_t FRbasis::LagrangeP(real_t x, int_t n, int_t N)
{
	switch (N)
	{
	case 0: return 1.0;
	case 1:
		switch (n)
		{
		case 0: return 0.5*(1.0 - x);
		case 1: return 0.5*(1.0 + x);
		default:
			ERROR("exceeding basis number");
			return 0.0;
		}
	case 2:
		switch (n)
		{
		case 0: return 0.5*(-x + pow(x, 2.0));
		case 1: return 1.0 - pow(x, 2.0);
		case 2: return 0.5*(x + pow(x, 2.0));
		default:
			ERROR("exceeding basis number");
			return 0.0;
		}
	case 3:
		switch (n)
		{
		case 0: return 0.125*(-1.0 + x + 5.0*pow(x, 2.0) - 5.0*pow(x, 3.0));
		case 1: return 0.625*(1.0 - 2.236067978*x - pow(x, 2.0) + 2.236067978*pow(x, 3.0));
		case 2: return 0.625*(1.0 + 2.236067978*x - pow(x, 2.0) - 2.236067978*pow(x, 3.0));
		case 3: return 0.125*(-1.0 - x + 5.0*pow(x, 2.0) + 5.0*pow(x, 3.0));
		default:
			ERROR("exceeding basis number");
			return 0.0;
		}
	default:
		ERROR("exceeds maximum order");
		return 0.0;
	}
}

real_t FRbasis::diff_LagrangeP(real_t x, int_t n, int_t N)
{
	switch (N)
	{
	case 0: return 0.0;
	case 1:
		switch (n)
		{
		case 0: return -0.5;
		case 1: return 0.5;
		default:
			ERROR("exceeding basis number");
			return 0.0;
		}
	case 2:
		switch (n)
		{
		case 0: return x - 0.5;
		case 1: return -2.0*x;
		case 2: return x + 0.5;
		default:
			ERROR("exceeding basis number");
			return 0.0;
		}
	case 3:
		switch (n)
		{
		case 0: return 0.125*(1.0 + 10.0*x - 15.0*pow(x, 2.0));
		case 1: return 0.625*(-2.236067978 - 2.0*x + 6.708203933*pow(x,2.0));
		case 2: return 0.625*(2.236067978 - 2.0*x - 6.708203933*pow(x, 2.0));
		case 3: return 0.125*(-1.0 + 10.0*x + 15.0*pow(x, 2.0));
		default:
			ERROR("exceeding basis number");
			return 0.0;
		}
	default:
		ERROR("exceeds maximum order");
		return 0.0;
	}
}

real_t FRbasis::toCompCoord(int_t icell, real_t x)
{
	return 2.0*(x - _grid->getCell()[icell]->getPosX()) / _sizeX;
}