#include "ConvFluxSW.h"

ConvFluxSW::ConvFluxSW(int_t polyOrder, Type correctFunc)
	: ConvFlux(polyOrder, correctFunc)
{

}

ConvFluxSW::~ConvFluxSW()
{

}

std::vector<real_t> ConvFluxSW::computeFlux(std::vector<real_t> begin, std::vector<real_t> end) const
{
	std::vector<real_t> flux_plus(3, 0.0), flux_minus(3, 0.0);
	std::vector<real_t> lambda_plus(3, 0.0), lambda_minus(3, 0.0);
	real_t a_plus = PHY_SOUND_SPEED(begin);
	real_t a_minus = PHY_SOUND_SPEED(end);
	real_t u_plus = begin[1] / begin[0];
	real_t u_minus = end[1] / end[0];

	// Computing characteristic velocity
	lambda_plus = PHY_CHAR_SPEED(begin);
	lambda_minus = PHY_CHAR_SPEED(end);
	for (int_t i = 0; i < 3; ++i)
	{
		lambda_plus[i] = lambda_plus[i] >= 0 ? lambda_plus[i] : 0.0;
		lambda_minus[i] = lambda_minus[i] <= 0 ? lambda_minus[i] : 0.0;
	}
	
	// Computing flux
	flux_plus[0] = lambda_plus[0] + 2.0*(GAMMA - 1)*lambda_plus[1] + lambda_plus[2];
	flux_minus[0] = lambda_minus[0] + 2.0*(GAMMA - 1)*lambda_minus[1] + lambda_minus[2];

	flux_plus[1] = (u_plus - a_plus)*lambda_plus[0] + 2.0*(GAMMA - 1)*u_plus*lambda_plus[1] + (u_plus + a_plus)*lambda_plus[2];
	flux_minus[1] = (u_minus - a_minus)*lambda_minus[0] + 2.0*(GAMMA - 1)*u_minus*lambda_minus[1] + (u_minus + a_minus)*lambda_minus[2];

	flux_plus[2] = (begin[2] / begin[0] + pow(a_plus, 2.0) / GAMMA - u_plus*a_plus)*lambda_plus[0] + (GAMMA - 1)*pow(u_plus, 2.0)*lambda_plus[1] + (begin[2] / begin[0] + pow(a_plus, 2.0) / GAMMA + u_plus*a_plus)*lambda_plus[2];
	flux_minus[2] = (end[2] / end[0] + pow(a_minus, 2.0) / GAMMA - u_minus*a_minus)*lambda_minus[0] + (GAMMA - 1)*pow(u_minus, 2.0)*lambda_minus[1] + (end[2] / end[0] + pow(a_minus, 2.0) / GAMMA + u_minus*a_minus)*lambda_minus[2];

	std::vector<real_t> flux(3, 0.0);
	for (int_t iflux = 0; iflux < 3; ++iflux)
	{
		flux_plus[iflux] *= 0.5*begin[0] / GAMMA;
		flux_minus[iflux] *= 0.5*end[0] / GAMMA;
		flux[iflux] = flux_plus[iflux] + flux_minus[iflux];
	}


	return flux;
}