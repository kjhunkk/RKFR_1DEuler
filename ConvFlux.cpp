#include "ConvFlux.h"

ConvFlux::ConvFlux(int_t polyOrder, Type correctFunc)
{
	_polyOrder = polyOrder;
	_correctFunc = correctFunc;
}

ConvFlux::~ConvFlux()
{

}

real_t ConvFlux::correctFunc(bool D, real_t z) const
{
	// Right/Left boundary correct function
	if (_correctFunc.compare("DG") == 0)
		return correctDG(D, z);
	else if (_correctFunc.compare("SG") == 0)
		return correctSG(D, z);
	else if ((_correctFunc.compare("2")*_correctFunc.compare("Lump,Lo")) == 0)
		return correct2(D, z);
	else if (_correctFunc.compare("Lo") == 0)
		return correctLo(D, z);
	else if (_correctFunc.compare("Ga") == 0)
		return correctGa(D, z);
	else
	{
		ERROR("cannot find correction function");
		return 0.0;
	}
}

real_t ConvFlux::DcorrectFunc(bool D, real_t z) const
{
	// Right/Left boundary correct function
	if (_correctFunc.compare("DG") == 0)
		return DcorrectDG(D, z);
	else if (_correctFunc.compare("SG") == 0)
		return DcorrectSG(D, z);
	else if ((_correctFunc.compare("2")*_correctFunc.compare("Lump,Lo")) == 0)
		return DcorrectLo(D, z);
	else if (_correctFunc.compare("Lo") == 0)
		return DcorrectLo(D, z);
	else if (_correctFunc.compare("Ga") == 0)
		return DcorrectGa(D, z);
	else
	{
		ERROR("cannot find correction function");
		return 0.0;
	}
}

real_t ConvFlux::correctDG(bool D, real_t z) const
{
	switch (_polyOrder)
	{
	case 1: return 0.5*(-0.5 - z + 1.5*pow(z, 2.0));
	case 2: return -0.5*(0.5 - 1.5*z - 1.5*pow(z, 2.0) + 2.5*pow(z, 3.0));
	case 3: return 0.5*(0.375 + 1.5*z - 3.75*pow(z, 2.0) - 2.5*pow(z, 3.0) + 4.375*pow(z, 4.0));
	case 4: return -0.5*(-0.375 + 1.875*z + 3.75*pow(z, 2.0) - 8.75*pow(z, 3.0) - 4.375*pow(z, 4.0) + 7.875*pow(z, 5.0));
	default:
		ERROR("exceeds maximum polynomial order");
		return 0.0;
	}
}

real_t ConvFlux::DcorrectDG(bool D, real_t z) const
{
	if (D == false)
	{
		switch (_polyOrder)
		{
		case 0: return -0.5;
		case 1: return -0.5 + 1.5*z;
		case 2: return 0.75 + 1.5*z - 3.75*pow(z, 2.0);
		case 3: return 0.75 - 3.75*z - 3.75*pow(z, 2.0) + 8.75*pow(z, 3.0);
		case 4: return -0.9375 - 3.75*z + 13.125*pow(z, 2.0) + 8.75*pow(z, 3.0) - 19.6875*pow(z, 4.0);
		default:
			ERROR("exceeds maximum polynomial order");
			return 0.0;
		}
	}
	else
	{
		switch (_polyOrder)
		{
		case 0: return 0.5;
		case 1: return 0.5 + 1.5*z;
		case 2: return -0.75 + 1.5*z + 3.75*pow(z, 2.0);
		case 3: return -0.75 - 3.75*z + 3.75*pow(z, 2.0) + 8.75*pow(z, 3.0);
		case 4: return 0.9375 - 3.75*z - 13.125*pow(z, 2.0) + 8.75*pow(z, 3.0) + 19.6875*pow(z, 4.0);
		default:
			ERROR("exceeds maximum polynomial order");
			return 0.0;
		}
	}
}

real_t ConvFlux::correctSG(bool D, real_t z) const
{
	return 0.0;
}

real_t ConvFlux::DcorrectSG(bool D, real_t z) const
{
	return 0.0;
}

real_t ConvFlux::correct2(bool D, real_t z) const
{
	return 0.0;
}

real_t ConvFlux::Dcorrect2(bool D, real_t z) const
{
	return 0.0;
}

real_t ConvFlux::correctLo(bool D, real_t z) const
{
	return 0.0;
}

real_t ConvFlux::DcorrectLo(bool D, real_t z) const
{
	return 0.0;
}

real_t ConvFlux::correctGa(bool D, real_t z) const
{
	return 0.0;
}

real_t ConvFlux::DcorrectGa(bool D, real_t z) const
{
	return 0.0;
}