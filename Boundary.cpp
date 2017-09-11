#include "Boundary.h"

Boundary::Boundary(Type type, std::shared_ptr<Zone> zone)
{
	// Initialize
	_num_cell = zone->getGrid()->getNumCell();
	_polyOrder = zone->getPolyOrder();
	_type = type;

	_begin = zone->getDescSolution()[GHOST];
	_end = zone->getDescSolution()[_num_cell - 1 - GHOST];
	_beginDOF.resize(_polyOrder + 1);
	_endDOF.resize(_polyOrder + 1);

	for (int_t idegree = 0; idegree <= zone->getPolyOrder(); ++idegree)
	{
		_beginDOF[idegree] = zone->getDOF()[idegree][GHOST];
		_endDOF[idegree] = zone->getDOF()[idegree][_num_cell - 1 - GHOST];
	}
}

Boundary::~Boundary()
{

}

void Boundary::apply(std::shared_ptr<Zone>& zone)
{
	std::vector<real_t> temp_solution = zone->getDescSolution();
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();

	if (_type == "constant")
	{
		for (int_t icell = 0; icell < GHOST; ++icell)
		{
			temp_solution[icell] = _begin;
			temp_solution[_num_cell - 1 - icell] = _end;
			for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
			{
				temp_DOF[idegree][icell] = _beginDOF[idegree];
				temp_DOF[idegree][_num_cell - 1 - icell] = _endDOF[idegree];
			}
		}
		zone->setDescSolution(temp_solution);
		zone->setDOF(temp_DOF);
	}

	else if (_type == "periodic")
	{
		for (int_t icell = 0; icell < GHOST; ++icell)
		{
			temp_solution[icell] = zone->getDescSolution()[_num_cell - 2*GHOST + icell];
			temp_solution[_num_cell - GHOST + icell] = zone->getDescSolution()[icell + GHOST];
			for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
			{
				temp_DOF[idegree][icell] = zone->getDOF()[idegree][_num_cell - 2*GHOST + icell];
				temp_DOF[idegree][_num_cell - GHOST + icell] = zone->getDOF()[idegree][icell + GHOST];
			}
		}
		zone->setDescSolution(temp_solution);
		zone->setDOF(temp_DOF);
	}

	else ERROR("cannot find proper boundary condition");
}