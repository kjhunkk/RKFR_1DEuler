#pragma once
#include "DataType.h"
#include "TimeInteg.h"

class TimeIntegRK : public TimeInteg
{
public:
	// Constructor / p.m. limiter type, limiting variable, CFL number, target time, convective flux, Zone(object), Boundary(object), RK order
	TimeIntegRK(Type, int_t, real_t, real_t, std::shared_ptr<ConvFlux>, std::vector<std::shared_ptr<Zone>>, std::vector<std::shared_ptr<Boundary>>, int_t);

	// Destructor
	virtual ~TimeIntegRK();

public:
	// Functions
	// Compute time integration / p.m. Zone(object) / r.t. go/stop
	virtual bool march(std::vector<std::shared_ptr<Zone>>);

protected:
	// Variables
	int_t _RKorder;
	// temporary DOF for TVD-RK / RK order, solution type number, DG degree, cell index
	std::vector<std::vector<std::vector<std::vector<real_t> > > > _temp_DOF;
};