#include "TimeInteg.h"

TimeInteg::TimeInteg(Type limiterType, int_t limitVar, real_t CFL, real_t targetTime, std::shared_ptr<ConvFlux> convFlux, std::vector<std::shared_ptr<Zone>> zone, std::vector<std::shared_ptr<Boundary>> bdry)
{
	// Initializing variables
	_limiterType = limiterType;
	_CFL = CFL; _targetTime = targetTime;
	_limitVar = limitVar;
	_currentTime = 0.0; _timeStep = 0.0;

	// Initializing objects
	_convFlux = convFlux;
	_zone = zone; _bdry = bdry;
	_basis = std::make_shared<FRbasis>(zone[0]->getGrid());
	_limiter = std::make_shared<Limiter>(_limiterType, zone, _limitVar);
	_pressureFix = std::make_shared<PressureFix>(5, zone);

	// Initializing temporary variables
	_temp_solution.resize(3);
	_prev_DOF.resize(3);
	_temp_RHS.resize(3);
	for (int_t isol = 0; isol < zone.size(); ++isol)
	{
		_temp_solution[isol].resize(zone[0]->getGrid()->getNumCell());
		_prev_DOF[isol].resize(zone[0]->getPolyOrder() + 1);
		_temp_RHS[isol].resize(zone[0]->getPolyOrder() + 1);
	}
	for (int_t isol = 0; isol < zone.size(); ++isol)
		for (int_t iorder = 0; iorder <= zone[0]->getPolyOrder(); ++iorder)
		{
			_prev_DOF[isol][iorder].resize(zone[0]->getGrid()->getNumCell());
			_temp_RHS[isol][iorder].resize(zone[0]->getGrid()->getNumCell());
		}
}

TimeInteg::~TimeInteg()
{

}

std::vector<std::vector<real_t> > TimeInteg::computeRHS(std::vector<std::shared_ptr<Zone>> zone, int_t solType) const
{
	real_t sizeX = zone[solType]->getGrid()->getSizeX();
	real_t inv_sizeX = 1.0/(zone[solType]->getGrid()->getSizeX());
	int_t num_cell = zone[solType]->getGrid()->getNumCell();
	int_t polyOrder = zone[solType]->getPolyOrder();

	// Temporary degree of freedom
	std::vector<std::vector<real_t> > RHS;
	RHS.resize(polyOrder + 1);
	for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		RHS[iorder].resize(num_cell);

	// Initializing RHS
	for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		for (int_t icell = 0; icell < num_cell; ++icell)
			RHS[iorder][icell] = 0.0;
	
	// FR flux
	std::vector<std::vector<real_t> > diff_flux;
	diff_flux.resize(polyOrder + 1);
	for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		diff_flux[iorder].resize(num_cell);
	real_t back_flux, for_flux; // computational flux
	real_t left_face, right_face; // face coordinate
	real_t diff_disc_flux; // derivative of discontinuous flux
	std::vector<real_t> left_value, mid_value, right_value;
	left_value.resize(3); mid_value.resize(3), right_value.resize(3);
	std::shared_ptr<FRbasis> basis = std::make_shared<FRbasis>(zone[solType]->getGrid());

	for (int_t icell = GHOST; icell < num_cell - GHOST; ++icell)
	{
		for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		{
			// Compute face coordinates
			left_face = zone[solType]->getGrid()->getCell()[icell]->getPosX() - 0.5*sizeX;
			right_face = zone[solType]->getGrid()->getCell()[icell]->getPosX() + 0.5*sizeX;

			// Computational flux
			// backward flux
			left_value = { zone[0]->getPolySolution(icell - 1, left_face), zone[1]->getPolySolution(icell - 1, left_face), zone[2]->getPolySolution(icell - 1, left_face) };
			right_value = { zone[0]->getPolySolution(icell, left_face), zone[1]->getPolySolution(icell, left_face), zone[2]->getPolySolution(icell, left_face) };
			back_flux = _convFlux->computeFlux(left_value, right_value)[solType];
			
			// forward flux
			left_value = { zone[0]->getPolySolution(icell, right_face), zone[1]->getPolySolution(icell, right_face), zone[2]->getPolySolution(icell, right_face) };
			right_value = { zone[0]->getPolySolution(icell + 1, right_face), zone[1]->getPolySolution(icell + 1, right_face), zone[2]->getPolySolution(icell + 1, right_face) };
			for_flux = _convFlux->computeFlux(left_value, right_value)[solType];
			
			// Compute derivative of discontinuous flux
			diff_disc_flux = 0.0;
			for (int_t iorder2 = 0; iorder2 <= polyOrder; ++iorder2)
			{
				for (int_t i = 0; i < 3; ++i)
					mid_value[i] = zone[i]->getDOF()[iorder2][icell];
				diff_disc_flux += ConvPhyFlux::phyFlux(mid_value)[solType]*basis->diff_LagrangeP(basis->getSolutionPoint(iorder, polyOrder), iorder2, polyOrder);
			}
						
			// Compute derivative of FR flux
			for (int_t i = 0; i < 3; ++i)
			{
				left_value[i] = zone[i]->getPolySolution(icell, left_face);
				right_value[i] = zone[i]->getPolySolution(icell, right_face);
			}

			diff_flux[iorder][icell] = diff_disc_flux
				+ (back_flux - ConvPhyFlux::phyFlux(left_value)[solType])*_convFlux->DcorrectFunc(false, basis->getSolutionPoint(iorder, polyOrder))
				+ (for_flux - ConvPhyFlux::phyFlux(right_value)[solType])*_convFlux->DcorrectFunc(true, basis->getSolutionPoint(iorder, polyOrder));

			// Coordinate conversion from computational to physical
			diff_flux[iorder][icell] *= 2.0*inv_sizeX;
		}
	}

	// Calculate RHS
	for (int_t icell = GHOST; icell < num_cell - GHOST; ++icell)
		for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
			RHS[iorder][icell] = -diff_flux[iorder][icell];
	
	return RHS;
}

void TimeInteg::computeTimeStep(std::vector<std::shared_ptr<Zone>> zone)
{
	std::vector<real_t> max_speed(zone[0]->getGrid()->getNumCell(), 0.0);
	std::vector<real_t> conv_value(3, 0.0);
	std::vector<real_t> char_speed(3, 0.0);
	for (int_t icell = GHOST; icell < max_speed.size() - GHOST; ++icell)
	{
		conv_value = { zone[0]->getDescSolution()[icell], zone[1]->getDescSolution()[icell], zone[2]->getDescSolution()[icell] };
		char_speed = PHY_CHAR_SPEED(conv_value);
		for (int_t i = 0; i < char_speed.size(); ++i)
			char_speed[i] = abs(char_speed[i]);
		max_speed[icell] = *std::max_element(char_speed.begin(), char_speed.end());
		_timeStep = _CFL*zone[0]->getGrid()->getSizeX() / *std::max_element(max_speed.begin(), max_speed.end()) / double(2 * zone[0]->getPolyOrder() + 1);
	}	
}

void TimeInteg::print() const
{
	MESSAGE("Marching finished.....");
	MESSAGE("CFL = " + std::to_string(_CFL));
	MESSAGE("Target time = " + std::to_string(_targetTime));
	MESSAGE("Final time step = " + std::to_string(_timeStep));
	MESSAGE("Current time = " + std::to_string(_currentTime));
}