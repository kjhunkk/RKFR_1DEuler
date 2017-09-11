#include "TimeIntegEuler.h"

TimeIntegEuler::TimeIntegEuler(Type limiterType, int_t limitVar, real_t CFL, real_t targetTime, std::shared_ptr<ConvFlux> convFlux, std::vector<std::shared_ptr<Zone>> zone, std::vector<std::shared_ptr<Boundary>> bdry)
	:TimeInteg(limiterType, limitVar, CFL, targetTime, convFlux, zone, bdry)
{
	_temp_DOF.resize(zone.size());
	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		_temp_DOF[inum].resize(zone[inum]->getPolyOrder() + 1);
		for (int_t idegree = 0; idegree <= zone[inum]->getPolyOrder(); ++idegree)
			_temp_DOF[inum][idegree].resize(zone[inum]->getGrid()->getNumCell());
	}
}

TimeIntegEuler::~TimeIntegEuler()
{

}

bool TimeIntegEuler::march(std::vector<std::shared_ptr<Zone>> zone)
{
	// Apply boundary condition
	for(int_t inum = 0; inum < zone.size(); ++inum)
		_bdry[inum]->apply(zone[inum]);

	// Marching starts
	bool procedure = true;
	if (abs(_currentTime) < epsilon) MESSAGE("Marching starts.....");

	// Calculate time step
	if ((_currentTime + _timeStep) > _targetTime)
	{
		_timeStep = _targetTime - _currentTime;
		procedure = false;
	}
	else computeTimeStep(zone);

	// Apply hMLP limiter
	_limiter->hMLP_Limiter(zone);

	// Pressure fix
	_pressureFix->apply(zone);

	// Save previous degree of freedom
	for(int_t inum = 0; inum < zone.size(); ++inum)
		_prev_DOF[inum] = zone[inum]->getDOF();

	// Calculate RHS
	for(int_t inum = 0; inum < zone.size(); ++inum)
		_temp_RHS[inum] = computeRHS(zone, inum);

	// Calculate DOF
	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		for (int_t idegree = 0; idegree <= zone[inum]->getPolyOrder(); ++idegree)
		{
			for (int_t icell = 0; icell < zone[inum]->getGrid()->getNumCell(); ++icell)
				_temp_DOF[inum][idegree][icell] = _prev_DOF[inum][idegree][icell] + _timeStep*_temp_RHS[inum][idegree][icell];
		}
	}

	for(int_t inum = 0; inum < zone.size(); ++inum)
		zone[inum]->setDOF(_temp_DOF[inum]);

	// Apply hMLP limiter
	_limiter->hMLP_Limiter(zone);

	// Pressure fix
	_pressureFix->apply(zone);

	// Calculate solution
	for(int_t inum = 0; inum < zone.size(); ++inum)
		zone[inum]->calSolution();

	// Update current time
	_currentTime += _timeStep;

	// Print finish condition
	if (!procedure) print();

	return procedure;
}