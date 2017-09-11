#pragma once
#include "DataType.h"

class InitialCondition
{
public:
	// Constructor / p.m. initial condition type
	InitialCondition(Type);

	// Destructor
	~InitialCondition();

public:
	// Functions
	inline Type getType() { return _type; }

	// calculate initial condition / p.m. x coordinate, solution type
	real_t initializer(real_t, int_t) const;

protected:
	// Variables
	Type _type;

protected:
	// Functions
	// Initial condition functions / p.m. x coordinate, solution type
	real_t smooth(real_t, int_t) const;
	real_t shocktube(real_t, int_t) const;
	real_t expansion(real_t, int_t) const;
	real_t blast(real_t, int_t) const;
	real_t double_shock(real_t, int_t) const;
	real_t contact(real_t, int_t) const;
	real_t blastInteraction(real_t, int_t) const;
	real_t Shu_Osher(real_t, int_t) const;
};