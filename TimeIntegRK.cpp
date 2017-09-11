#include "TimeIntegRK.h"

TimeIntegRK::TimeIntegRK(Type limiterType, int_t limitVar, real_t CFL, real_t targetTime, std::shared_ptr<ConvFlux> convFlux, std::vector<std::shared_ptr<Zone>> zone, std::vector<std::shared_ptr<Boundary>> bdry, int_t RKorder)
	:TimeInteg(limiterType, limitVar, CFL, targetTime, convFlux, zone, bdry)
{
	_RKorder = RKorder;
	_temp_DOF.resize(RKorder);
	for (int_t iorder = 0; iorder < RKorder; ++iorder)
	{
		_temp_DOF[iorder].resize(zone.size());
		for (int_t inum = 0; inum < zone.size(); ++inum)
		{
			_temp_DOF[iorder][inum].resize(zone[inum]->getPolyOrder() + 1);
			for (int_t idegree = 0; idegree <= zone[inum]->getPolyOrder(); ++idegree)
				_temp_DOF[iorder][inum][idegree].resize(zone[inum]->getGrid()->getNumCell());
		}
	}
}

TimeIntegRK::~TimeIntegRK()
{

}

bool TimeIntegRK::march(std::vector<std::shared_ptr<Zone>> zone)
{
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

	// Declare temprorary Zone for TVD Runge-Kutta time integration
	std::vector<std::shared_ptr<Zone>> temp_zone;
	temp_zone.resize(zone.size());
	for(int_t i = 0; i < zone.size(); ++i)
		temp_zone[i] = std::make_shared<Zone>(*zone[i]);

	// ----------------------First step--------------------------
	// Apply boundary condition
	for(int_t inum = 0; inum < zone.size(); ++inum)
		_bdry[inum]->apply(temp_zone[inum]);

	// Apply hMLP limiter
	_limiter->hMLP_Limiter(temp_zone);

	// Pressure fix
	_pressureFix->apply(temp_zone);

	// Save previous degree of freedom
	for(int_t inum = 0; inum < zone.size(); ++inum)
		_prev_DOF[inum] = temp_zone[inum]->getDOF();

	// Calculate RHS
	for(int_t inum = 0; inum < zone.size(); ++inum)
		_temp_RHS[inum] = computeRHS(temp_zone, inum);

	// Calculate DOF
	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		for (int_t idegree = 0; idegree <= zone[inum]->getPolyOrder(); ++idegree)
		{
			for (int_t icell = GHOST; icell < zone[inum]->getGrid()->getNumCell() - GHOST; ++icell)
				_temp_DOF[0][inum][idegree][icell] = _prev_DOF[inum][idegree][icell] + _timeStep*_temp_RHS[inum][idegree][icell];
		}
	}	

	// Update temporary Zone object
	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		temp_zone[inum]->setDOF(_temp_DOF[0][inum]);
		temp_zone[inum]->calSolution();
	}

	// ---------------------Second step---------------------------
	// Apply boundary condition
	for (int_t inum = 0; inum < zone.size(); ++inum)
		_bdry[inum]->apply(temp_zone[inum]);

	// Apply hMLP limiter
	_limiter->hMLP_Limiter(temp_zone);

	// Pressure fix
	_pressureFix->apply(temp_zone);

	// Calculate RHS
	for (int_t inum = 0; inum < zone.size(); ++inum)
		_temp_RHS[inum] = computeRHS(temp_zone, inum);

	// Calculate DOF
	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		for (int_t idegree = 0; idegree <= zone[inum]->getPolyOrder(); ++idegree)
		{
			for (int_t icell = GHOST; icell < zone[inum]->getGrid()->getNumCell() - GHOST; ++icell)
				_temp_DOF[1][inum][idegree][icell] = 0.75*_prev_DOF[inum][idegree][icell] + 0.25*(_temp_DOF[0][inum][idegree][icell] + _timeStep*_temp_RHS[inum][idegree][icell]);
		}
	}

	// Update temporary Zone object
	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		temp_zone[inum]->setDOF(_temp_DOF[1][inum]);
		temp_zone[inum]->calSolution();
	}

	// ----------------------Third step---------------------------
	// Apply boundary condition
	for (int_t inum = 0; inum < zone.size(); ++inum)
		_bdry[inum]->apply(temp_zone[inum]);

	// Apply hMLP limiter
	_limiter->hMLP_Limiter(temp_zone);

	// Pressure fix
	_pressureFix->apply(temp_zone);

	// Calculate RHS
	for (int_t inum = 0; inum < zone.size(); ++inum)
		_temp_RHS[inum] = computeRHS(temp_zone, inum);

	// Calculate DOF
	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		for (int_t idegree = 0; idegree <= zone[inum]->getPolyOrder(); ++idegree)
		{
			for (int_t icell = GHOST; icell < zone[inum]->getGrid()->getNumCell() - GHOST; ++icell)
				_temp_DOF[2][inum][idegree][icell] = CONST13*_prev_DOF[inum][idegree][icell] + CONST23*(_temp_DOF[1][inum][idegree][icell] + _timeStep*_temp_RHS[inum][idegree][icell]);
		}
	}

	// Update temporary Zone object
	for (int_t inum = 0; inum < zone.size(); ++inum)
		temp_zone[inum]->setDOF(_temp_DOF[2][inum]);

	// Apply boundary condition
	for (int_t inum = 0; inum < zone.size(); ++inum)
		_bdry[inum]->apply(temp_zone[inum]);

	// Apply hMLP limiter
	_limiter->hMLP_Limiter(temp_zone);

	// Pressure fix
	_pressureFix->apply(temp_zone);

	// Update solution zone
	for (int_t inum = 0; inum < zone.size(); ++inum)
		zone[inum]->setDOF(temp_zone[inum]->getDOF());

	// Calculate solution
	for (int_t inum = 0; inum < zone.size(); ++inum)
		zone[inum]->calSolution();

	// Update current time
	_currentTime += _timeStep;

	// Print finish condition
	if (!procedure) print();

	return procedure;
}