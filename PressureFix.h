#pragma once
#include "DataType.h"
#include "Zone.h"
#include "FRbasis.h"

class PressureFix
{
public:
	// Constructor / p.m. max iteration, Zone(Object)
	PressureFix(int_t, std::vector<std::shared_ptr<Zone>>);

	// Destructor
	~PressureFix();

public:
	// Functions
	// Apply pressure fix / p.m. Zone(Object)
	void apply(std::vector<std::shared_ptr<Zone>>);

protected:
	// Variables
	int_t _maxIter;
	int_t _polyOrder;
	int_t _num_cell;
	std::vector<std::vector<bool>> _marker;

protected:
	// Functions
	// Mark negative pressure / p.m. Zone(Object) / r.t. end of iteration(true = stop)
	bool marking(std::vector<std::shared_ptr<Zone>>);

	// Apply pressure fix once / p.m. Zone(Object), iteration number
	void press(std::vector<std::shared_ptr<Zone>>);

	// Destroy high order term / p.m. Zone(Object)
	void destroy(std::vector<std::shared_ptr<Zone>>);
};