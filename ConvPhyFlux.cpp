#include "ConvPhyFlux.h"

ConvPhyFlux::ConvPhyFlux()
{

}

ConvPhyFlux::~ConvPhyFlux()
{

}

std::vector<real_t> ConvPhyFlux::phyFlux(std::vector<real_t> u)
{
	std::vector<real_t> flux(3, 0.0);
	flux[0] = u[1];
	flux[1] = 0.5*(3.0 - GAMMA)*pow(u[1], 2.0) / u[0] + (GAMMA - 1)*u[2];
	flux[2] = GAMMA*u[1] * u[2] / u[0] - 0.5*(GAMMA - 1)*pow(u[1], 3.0) / pow(u[0], 2.0);
	return flux;
}

std::vector<real_t> ConvPhyFlux::phyCharSpeed(std::vector<real_t> u)
{
	real_t a = phyA(u);
	real_t velocity = u[1] / u[0];
	std::vector<real_t> lambda(3, 0.0);
	lambda[0] = velocity - a;
	lambda[1] = velocity;
	lambda[2] = velocity + a;
	return lambda;
}

real_t ConvPhyFlux::phyA(std::vector<real_t> u)
{
	return sqrt(GAMMA*(GAMMA - 1)*(u[2] / u[0] - 0.5*pow(u[1] / u[0], 2.0)));
}

ConvPhyFlux ConvPhyFlux::_phyFlux;