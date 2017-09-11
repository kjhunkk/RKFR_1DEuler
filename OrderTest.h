#pragma once
#include "DataType.h"
#include "Zone.h"
#include "Grid.h"

class OrderTest
{
public:
	// Constructor
	OrderTest();

	// Destructor
	~OrderTest();

public:
	// Functions
	// Set exact solution
	void setExact(std::vector<real_t> exact) { _exact = exact; _num = exact.size(); }

	// Convert Zone to solution array / p.m. Zone
	std::vector<real_t> ZoneToArray(std::shared_ptr<Zone>);

	// Convert Zone to polynomical solution array / p.m. Zone
	std::vector<real_t> ZoneToPoly(std::shared_ptr<Zone>);

	// Calculate error / p.m. computed solution
	real_t L1error(std::vector<real_t>);

	real_t L2error(std::vector<real_t>);

	real_t Linf_error(std::vector<real_t>);

protected:
	// Variables
	std::vector<real_t> _exact;
	int_t _num;
};