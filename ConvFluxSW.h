#pragma once
#include "DataType.h"
#include "ConvFlux.h"
#include "ConvPhyFlux.h"

class ConvFluxSW : public ConvFlux
{
public:
	// Constructor / p.m. polynomial order, correction function type
	ConvFluxSW(int_t, Type);

	// Destructor
	virtual ~ConvFluxSW();

public:
	// Functions
	// Compute flux / p.m. begin, end
	virtual std::vector<real_t> computeFlux(std::vector<real_t>, std::vector<real_t>) const;
};