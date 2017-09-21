#include "ConvFluxRoe.h"

ConvFluxRoe::ConvFluxRoe(int_t polyOrder, Type correctFunc)
	: ConvFlux(polyOrder, correctFunc)
{

}

ConvFluxRoe::~ConvFluxRoe()
{

}

std::vector<real_t> ConvFluxRoe::computeFlux(std::vector<real_t> begin, std::vector<real_t> end) const
{
	std::vector<real_t> comp_flux(3, 0.0);
	std::vector<real_t> lambda(3, 0.0);
	std::vector<real_t> wave_strength(3, 0.0);
	std::vector<std::vector<real_t>> K; // Averaged eigen vectors
	K.resize(3);
	for (int_t i = 0; i < K.size(); ++i)
		K[i].resize(3);
	real_t u, H, a; // Roe averaged variables

	// Compute Roe averaged variables
	u = (sqrt(end[0])*begin[1] + sqrt(begin[0])*end[1]) / (begin[0] * sqrt(end[0]) + end[0] * sqrt(begin[0]));
	H = (sqrt(begin[0])*computeH(begin[0], begin[1], begin[2]) + sqrt(end[0])*computeH(end[0], end[1], end[2])) / (sqrt(begin[0]) + sqrt(end[0]));
	a = sqrt((GAMMA - 1)*(H - 0.5*u*u));

	// Compute averaged eigenvalues
	lambda[0] = u - a;
	lambda[1] = u;
	lambda[2] = u + a;

	// Compute averaged eigen vectors
	K[0] = { 1, u - a, H - u*a };
	K[1] = { 1, u, 0.5*u*u };
	K[2] = { 1, u + a, H + u*a };

	// Compute wave strength
	wave_strength[1] = (GAMMA - 1)*((end[0] - begin[0])*(H - u*u) + u*(end[1] - begin[1]) - (end[2] - begin[2])) / pow(a, 2.0);
	wave_strength[0] = 0.5*((end[0] - begin[0])*(u + a) - (end[1] - begin[1]) - a*wave_strength[1]) / a;
	wave_strength[2] = end[0] - begin[0] - (wave_strength[0] + wave_strength[1]);

	// Harten-Hyman Entropy fix
	lambda = entropyFix(lambda, wave_strength, u, H, a, begin, end);

	// Compute Roe flux
	comp_flux = PHY_FLUX(begin);
	for (int_t chari = 0; chari < K.size(); ++chari)
		if (lambda[chari] <= 0.0)
			for (int_t sysi = 0; sysi < K[chari].size(); ++sysi)
				comp_flux[sysi] += wave_strength[chari] * lambda[chari] * K[chari][sysi];

	return comp_flux;
}

real_t ConvFluxRoe::computeH(real_t u1, real_t u2, real_t u3) const
{
	return GAMMA*(u3 / u1) - 0.5*(GAMMA - 1)*pow(u2 / u1, 2.0);
}

std::vector<real_t> ConvFluxRoe::entropyFix(std::vector<real_t> lambda, std::vector<real_t> wave_strength, real_t u, real_t H, real_t a, std::vector<real_t> begin, std::vector<real_t> end) const
{
	real_t u_starL, u_starR;
	real_t a_starL, a_starR;
	real_t a_L, a_R;
	a_L = sqrt(GAMMA*(GAMMA - 1)*(begin[2] / begin[0] - 0.5*pow(begin[1] / begin[0], 2.0)));
	a_R = sqrt(GAMMA*(GAMMA - 1)*(end[2] / end[0] - 0.5*pow(end[1] / end[0], 2.0)));
	u_starL = (begin[1] + wave_strength[0] * (u - a)) / (begin[0] + wave_strength[0]);
	a_starL = sqrt((GAMMA*(GAMMA - 1)*(begin[2] + wave_strength[0] * (H - u*a) - 0.5*(begin[0] + wave_strength[0])*u_starL*u_starL)) / (begin[0] + wave_strength[0]));
	u_starR = (end[1] - wave_strength[2] * (u + a)) / (end[0] - wave_strength[2]);
	a_starR = sqrt((GAMMA*(GAMMA - 1)*(end[2] - wave_strength[2] * (H + u*a) - 0.5*(end[0] - wave_strength[2])*u_starR*u_starR)) / (end[0] - wave_strength[2]));
	real_t lambda1L, lambda1R, lambda3L, lambda3R;
	lambda1L = begin[1] / begin[0] - a_L;
	lambda1R = u_starL - a_starL;
	lambda3L = u_starR + a_starR;
	lambda3R = end[1] / end[0] + a_R;

	std::vector<real_t> fix_lambda = lambda;
	// Left transonic rarefaction
	if ((lambda1L < 0.0) && (0.0 < lambda1R))
		fix_lambda[0] = lambda1L*(lambda1R - lambda[0]) / (lambda1R - lambda1L);

	// Right transonic rarefaction
	if ((lambda3L < 0.0) && (0.0 < lambda3R))
		fix_lambda[2] = lambda3R*(lambda[2] - lambda3L) / (lambda3R - lambda3L);

	return fix_lambda;
}