#include "Quadrature.h"

Quadrature::Quadrature()
{

}

Quadrature::~Quadrature()
{

}

real_t Quadrature::gauss2X(int_t index)
{
	switch (index)
	{
	case 0: return -CONST13*sqrt(3.0);
	case 1: return CONST13*sqrt(3.0);
	default: return 0.0;
	}
}

real_t Quadrature::gauss3X(int_t index)
{
	switch (index)
	{
	case 0: return -0.2*sqrt(15.0);
	case 1: return 0.0;
	case 2: return 0.2*sqrt(15.0);
	default: return 0.0;
	}
}

real_t Quadrature::gauss4X(int_t index)
{
	switch (index)
	{
	case 0: return -0.861136311594053;
	case 1: return -0.339981043584856;
	case 2: return 0.339981043584856;
	case 3: return 0.861136311594053;
	default: return 0.0;
	}
}

real_t Quadrature::gauss5X(int_t index)
{
	switch (index)
	{
	case 0: return -0.906179845938664;
	case 1: return -0.538469310105683;
	case 2: return 0.0;
	case 3: return 0.538469310105683;
	case 4: return 0.906179845938664;
	default: return 0.0;
	}
}

real_t Quadrature::lobatto2X(int_t index)
{
	switch (index)
	{
	case 0: return -1.0;
	case 1: return 1.0;
	default: return 0.0;
	}
}

real_t Quadrature::lobatto3X(int_t index)
{
	switch (index)
	{
	case 0: return -1.0;
	case 1: return 0.0;
	case 2: return 1.0;
	default: return 0.0;
	}
}

real_t Quadrature::lobatto4X(int_t index)
{
	switch (index)
	{
	case 0: return -1.0;
	case 1: return -0.447213595499958;
	case 2: return 0.447213595499958;
	case 3: return 1.0;
	default: return 0.0;
	}
}

real_t Quadrature::lobatto5X(int_t index)
{
	switch (index)
	{
	case 0: return -1.0;
	case 1: return -0.654653670707977;
	case 2: return 0.0;
	case 3: return 0.654653670707977;
	case 4: return 1.0;
	default: return 0.0;
	}
}

real_t Quadrature::gauss2W(int_t index)
{
	return 1.0;
}

real_t Quadrature::gauss3W(int_t index)
{
	switch (index)
	{
	case 0: return CONST59;
	case 1: return CONST89;
	case 2: return CONST59;
	default: return 0.0;
	}
}

real_t Quadrature::gauss4W(int_t index)
{
	switch (index)
	{
	case 0: return 0.347854845137454;
	case 1: return 0.652145154862546;
	case 2: return 0.652145154862546;
	case 3: return 0.347854845137454;
	default: return 0.0;
	}
}

real_t Quadrature::gauss5W(int_t index)
{
	switch (index)
	{
	case 0: return 0.236926885056189;
	case 1: return 0.478628670499366;
	case 2: return 0.568888888888889;
	case 3: return 0.478628670499367;
	case 4: return 0.236926885056189;
	default: return 0.0;
	}
}

real_t Quadrature::legendreP(int_t Pn, real_t x)
{
	switch (Pn)
	{
	case 0: return 1.0;
	case 1: return x;
	case 2: return -0.5 + 1.5*pow(x, 2.0);
	case 3: return -1.5*x + 2.5*pow(x, 3.0);
	default:
	{
		ERROR("exceed maximum order");
		return 0.0;
	}
	}
}

Quadrature Quadrature::_quad;