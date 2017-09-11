#include "Zone.h"

Zone::Zone(std::shared_ptr<Grid> grid, int_t solType)
{
	_grid = grid;
	_solType = solType;
	_polyOrder = 0;
	_solution.resize(grid->getNumCell());
	_DOF.resize(1);
	_DOF[0].resize(grid->getNumCell());

	// DG basis
	_basis = std::make_shared<FRbasis>(_grid);
}

Zone::Zone(std::shared_ptr<Grid> grid, int_t polyOrder, int_t solType)
{
	_grid = grid;
	_solType = solType;
	_polyOrder = polyOrder;

	// Print error message if polynomial order is less than 0
	if (_polyOrder < 0)
		ERROR("Polynomial order");

	_solution.resize(grid->getNumCell());
	_DOF.resize(polyOrder + 1);
	for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		_DOF[iorder].resize(grid->getNumCell());

	// DG basis
	_basis = std::make_shared<FRbasis>(_grid);
}

Zone::~Zone()
{

}

real_t Zone::getPolySolution(int_t icell, real_t x) const
{
	real_t u = 0;
	real_t tester = 0;

	// Calculate polynomial solution at x / beware of time level
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		u += _DOF[idegree][icell] * _basis->LagrangeP(_basis->toCompCoord(icell,x), idegree, _polyOrder);

	return u;
}

real_t Zone::getPolySolution(int_t icell, real_t x, std::shared_ptr<Zone> zone) const
{
	real_t u = 0;
	// Calculate polynomial solution at x / beware of time level
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		u += zone->getDOF()[idegree][icell] * _basis->LagrangeP(_basis->toCompCoord(icell,x), idegree, _polyOrder);

	return u;
}

void Zone::initialize(std::shared_ptr<InitialCondition> initialCondition)
{
	// Temporary cell object
	std::vector<std::shared_ptr<Cell> > temp_cell = _grid->getCell();

	// Grid size
	real_t temp_dx = _grid->getSizeX();
	real_t tiny;

	// Initializing Degree of freedom
	for (int_t icell = 0; icell < _grid->getNumCell(); ++icell)
	{
		real_t temp_x = temp_cell[icell]->getPosX();
		_DOF[0][icell] = initialCondition->initializer(temp_x + 0.5*temp_dx*_basis->getSolutionPoint(0, _polyOrder) + epsilon, _solType);
		_DOF[_polyOrder][icell] = initialCondition->initializer(temp_x + 0.5*temp_dx*_basis->getSolutionPoint(_polyOrder, _polyOrder) - epsilon, _solType);
		if (_polyOrder > 1)
		{
			for (int_t iorder = 1; iorder < _polyOrder; ++iorder)
				_DOF[iorder][icell] = initialCondition->initializer(temp_x + 0.5*temp_dx*_basis->getSolutionPoint(iorder, _polyOrder), _solType);
		}
	}

	// Calculate solution
	calSolution();
}

void Zone::calSolution()
{
	for (int_t icell = 0; icell < _grid->getNumCell(); ++icell)
		_solution[icell] = getAverage(icell);
}

real_t Zone::getAverage(int_t num_cell) const
{
	// Temporary object
	std::vector<real_t> nodal_DOF(_polyOrder + 1, 0.0);
	for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
		nodal_DOF[iorder] = _DOF[iorder][num_cell];

	// Vandermonde matrix
	std::vector<std::vector<real_t>> Vmatrix, inv_Vmatrix;
	std::vector<real_t> lobatto_points(_polyOrder + 1, 0.0);
	switch (_polyOrder)
	{
	case 0:
		lobatto_points[0] = 0.0;
		break;
	case 1:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto2_X(iorder);
		break;
	case 2:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto3_X(iorder);
		break;
	case 3:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto4_X(iorder);
		break;
	case 4:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto5_X(iorder);
		break;
	default: ERROR("exceeds maximum order");
	}
	Vmatrix = Vandermonde(lobatto_points);
	inv_Vmatrix = inv_Vandermonde(lobatto_points);

	// Convert nodal to modal expression
	std::vector<real_t> modal_DOF(_polyOrder + 1, 0.0);
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		for (int_t jdegree = 0; jdegree <= _polyOrder; ++jdegree)
			modal_DOF[idegree] += inv_Vmatrix[idegree][jdegree] * nodal_DOF[jdegree];

	// return mode 0
	return modal_DOF[0];
}

std::vector<std::vector<real_t>> Zone::Vandermonde(std::vector<real_t> r) const
{
	std::vector<std::vector<real_t>> V;
	V.resize(r.size());
	for (int_t ipoint = 0; ipoint < r.size(); ++ipoint)
		V[ipoint].resize(r.size());

	for (int_t iorder = 0; iorder < r.size(); ++iorder)
		for (int_t ipoint = 0; ipoint < r.size(); ++ipoint)
			V[ipoint][iorder] = LegendreP(iorder, r[ipoint]);

	return V;
}

std::vector<std::vector<real_t>> Zone::inv_Vandermonde(std::vector<real_t> r) const
{
	std::vector<std::vector<real_t>> invV;
	invV.resize(r.size());
	for (int_t ipoint = 0; ipoint < r.size(); ++ipoint)
		invV[ipoint].resize(r.size());

	real_t den1, den2, den3, den4;
	switch (r.size())
	{
	case 1:
		invV[0] = { 1.0 };
		break;
	case 2:
		den1 = 1.0 / (r[0] - r[1]);
		invV[0] = { -r[1] * den1, r[0] * den1 };
		invV[1] = { den1, -den1 };
		break;
	case 3:
		den1 = 1.0 / (r[0] * r[1] + r[0] * r[2] - r[1] * r[2] - pow(r[0], 2.0));
		den2 = 1.0 / (r[0] * r[1] - r[0] * r[2] + r[1] * r[2] - pow(r[1], 2.0));
		den3 = 1.0 / (r[0] * r[1] - r[0] * r[2] - r[1] * r[2] + pow(r[2], 2.0));
		invV[0] = { -CONST13*(3.0*r[1] * r[2] + 1.0)*den1, -CONST13*(3.0*r[0] * r[2] + 1.0)*den2, CONST13*(3.0*r[0] * r[1] + 1.0)*den3 };
		invV[1] = { (r[1] + r[2])*den1, (r[0] + r[2])*den2, -(r[0] + r[1])*den3 };
		invV[2] = { -CONST23*den1, -CONST23*den2, CONST23*den3 };
		break;
	case 4:
		den1 = 1.0 / (pow(r[0], 2.0)*(r[1] + r[2] + r[3] - r[0]) - (r[3] * r[0] * r[1]) - (r[3] * r[0] * r[2]) + (r[3] * r[1] * r[2]) - (r[0] * r[1] * r[2]));
		den2 = 1.0 / (pow(r[1], 2.0)*(r[0] + r[2] + r[3] - r[1]) - (r[3] * r[0] * r[1]) + (r[3] * r[0] * r[2]) - (r[3] * r[1] * r[2]) - (r[0] * r[1] * r[2]));
		den3 = 1.0 / (pow(r[2], 2.0)*(r[0] + r[1] + r[3] - r[2]) + (r[3] * r[0] * r[1]) - (r[3] * r[0] * r[2]) - (r[3] * r[1] * r[2]) - (r[0] * r[1] * r[2]));
		den4 = 1.0 / (pow(r[3], 2.0)*(r[0] + r[1] + r[2] - r[3]) - (r[3] * r[0] * r[1]) - (r[3] * r[0] * r[2]) - (r[3] * r[1] * r[2]) + (r[0] * r[1] * r[2]));
		invV[0] = { CONST13*(r[3] + r[1] + r[2] + 3.0*r[3] * r[1] * r[2])*den1, CONST13*(r[3] + r[0] + r[2] + 3.0*r[3] * r[0] * r[2])*den2, CONST13*(r[3] + r[0] + r[1] + 3.0*r[3] * r[0] * r[1])*den3, CONST13*(r[0] + r[1] + r[2] + 3.0*r[0] * r[1] * r[2])*den4 };
		invV[1] = { -0.2*(5.0*r[3] * r[1] + 5.0*r[3] * r[2] + 5.0*r[1] * r[2] + 3.0)*den1, -0.2*(5.0*r[3] * r[0] + 5.0*r[3] * r[2] + 5.0*r[0] * r[2] + 3.0)*den2, -0.2*(5.0*r[3] * r[0] + 5.0*r[3] * r[1] + 5.0*r[0] * r[1] + 3.0)*den3, -0.2*(5.0*r[0] * r[1] + 5.0*r[0] * r[2] + 5.0*r[1] * r[2] + 3.0)*den4 };
		invV[2] = { CONST23*(r[3] + r[1] + r[2])*den1, CONST23*(r[3] + r[0] + r[2])*den2, CONST23*(r[3] + r[0] + r[1])*den3, CONST23*(r[0] + r[1] + r[2])*den4 };
		invV[3] = { -0.4*den1, -0.4*den2, -0.4*den3, -0.4*den4 };
		break;
	default: ERROR("exceeds maximum order");
	}
	return invV;
}

void Zone::print() const
{
	std::cout << "X coordinate" << "\t\t" << "Solution" << "\t\t";
	for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
		std::cout << "DOF(" << iorder << ")\t\t\t";
	std::cout << "\n";

	// Print variables
	for (int_t icell = 0; icell < _grid->getNumCell(); ++icell)
	{
		std::cout << std::to_string(_grid->getCell()[icell]->getPosX()) << "\t\t" << std::to_string(_solution[icell]) << "\t\t";
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			std::cout << std::to_string(_DOF[iorder][icell]) << "\t\t";
		std::cout << "\n";
	}
}