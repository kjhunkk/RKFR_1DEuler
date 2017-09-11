#pragma once
#include "DataType.h"
#include "Zone.h"

class Boundary
{
public:
	// Consturctor(should be constructed after zone initilizing) / p.m. boundary condition, Zone(object)
	Boundary(Type, std::shared_ptr<Zone>);

	// Destructor
	~Boundary();

public:
	// Functions
	void apply(std::shared_ptr<Zone>&);

protected:
	// Variables
	int_t _num_cell;
	int_t _polyOrder;
	Type _type;
	real_t _begin;
	real_t _end;
	std::vector<real_t> _beginDOF;
	std::vector<real_t> _endDOF;
};