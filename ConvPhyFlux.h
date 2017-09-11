#pragma once
#include "DataType.h"

class ConvPhyFlux
{
protected:
	ConvPhyFlux();
	~ConvPhyFlux();

public:
	// Calculate Physical flux / p.m. solution variable
	static std::vector<real_t> phyFlux(std::vector<real_t>);

	// Compute characteristic speed / p.m. solution variable
	static std::vector<real_t> phyCharSpeed(std::vector<real_t>);

	// Compute speed of sound / p.m. solution variable
	static real_t phyA(std::vector<real_t>);

private:
	// ConvPhyFlux variable
	static ConvPhyFlux _phyFlux;
};

// Convective Physical Flux macro
#define PHY_FLUX(u) ConvPhyFlux::phyFlux(u)

#define PHY_CHAR_SPEED(u) ConvPhyFlux::phyCharSpeed(u)

#define PHY_SOUND_SPEED(u) ConvPhyFlux::phyA(u)