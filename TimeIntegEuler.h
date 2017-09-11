#pragma once
#include "DataType.h"
#include "TimeInteg.h"

class TimeIntegEuler : public TimeInteg
{
public:
	// Constructor / p.m. limiter type, limiting variable, CFL number, target time, convective flux, Zone(object), Boundary(object)
	TimeIntegEuler(Type, int_t, real_t, real_t, std::shared_ptr<ConvFlux>, std::vector<std::shared_ptr<Zone>>, std::vector<std::shared_ptr<Boundary>>);

	// Destructor
	virtual ~TimeIntegEuler();

public:
	// Functions
	// Compute time integration / p.m. Zone(object) / r.t. go/stop
	virtual bool march(std::vector<std::shared_ptr<Zone>>);	

protected:
	// Variables
	std::vector<std::vector<std::vector<real_t> > > _temp_DOF;
};