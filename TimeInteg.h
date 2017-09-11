#pragma once
#include "DataType.h"
#include "Zone.h"
#include "ConvFlux.h"
#include "Boundary.h"
#include "Limiter.h"
#include "FRbasis.h"
#include "ConvPhyFlux.h"
#include "PressureFix.h"

class TimeInteg
{
public:
	// Constructor / p.m. limiter type, limiting variable, CFL number, target time, convective flux, Zone(object), Boundary(object)
	TimeInteg(Type, int_t, real_t, real_t, std::shared_ptr<ConvFlux>, std::vector<std::shared_ptr<Zone>>, std::vector<std::shared_ptr<Boundary>>);

	// Destructor
	virtual ~TimeInteg();

public:
	// Functions
	inline void reset() { _currentTime = 0.0; }

	inline void reset(real_t target) { _currentTime = 0.0; _targetTime = target; }

	inline real_t getTime() const { return _currentTime; }

	inline real_t getTimeStep() const { return _timeStep; }

	inline real_t getTargetTime() const { return _targetTime; }

	// Compute time integration / p.m. Zone(object) / r.t. go/stop
	virtual bool march(std::vector<std::shared_ptr<Zone>>) = 0;

protected:
	// Variables
	std::vector<std::vector<real_t> > _temp_solution;
	std::vector<std::vector<std::vector<real_t> > > _prev_DOF;
	std::vector<std::vector<std::vector<real_t> > > _temp_RHS;
	std::vector<std::shared_ptr<Zone>> _zone;
	std::shared_ptr<ConvFlux> _convFlux;
	std::vector<std::shared_ptr<Boundary>> _bdry;
	std::shared_ptr<FRbasis> _basis;
	std::shared_ptr<Limiter> _limiter;
	std::shared_ptr<PressureFix> _pressureFix;
	Type _limiterType;
	int_t _limitVar;
	real_t _CFL;
	real_t _currentTime;
	real_t _timeStep;
	real_t _targetTime;

protected:
	// Functions
	// Compute right hand side / p.m. Zone to compute, solution type to compute
	std::vector<std::vector<real_t> > computeRHS(std::vector<std::shared_ptr<Zone>>, int_t) const;

	// Compute time step / p.m. Zone(object)
	void computeTimeStep(std::vector<std::shared_ptr<Zone>>);

	// Print time variables
	void print() const;
};