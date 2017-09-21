#pragma once
#include "DataType.h"
#include "ConvFlux.h"
#include "ConvPhyFlux.h"

class ConvFluxRoe : public ConvFlux
{
public:
	// Constructor / p.m. polynomial order, correction function type
	ConvFluxRoe(int_t, Type);

	// Destructor
	virtual ~ConvFluxRoe();

public:
	// Functions
	// Compute flux / p.m. begin, end
	virtual std::vector<real_t> computeFlux(std::vector<real_t>, std::vector<real_t>) const;

	// Compute total enthalpy / p.m. rho, rho*u, E
	real_t computeH(real_t, real_t, real_t) const;

	// Harten-Hyman Entropy fix / p.m. averaged eigenvalues, wave strength, Roe averaged variables(u,H,a), begin, end
	std::vector<real_t> entropyFix(std::vector<real_t>, std::vector<real_t>, real_t, real_t, real_t, std::vector<real_t>, std::vector<real_t>) const;
};