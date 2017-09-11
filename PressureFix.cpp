#include "PressureFix.h"

PressureFix::PressureFix(int_t maxIter, std::vector<std::shared_ptr<Zone>> zone)
{
	_maxIter = maxIter;
	_polyOrder = zone[0]->getPolyOrder();
	_num_cell = zone[0]->getGrid()->getNumCell();
	_marker.resize(_polyOrder + 1);
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		_marker[idegree].resize(_num_cell);
}

PressureFix::~PressureFix()
{

}

void PressureFix::apply(std::vector<std::shared_ptr<Zone>> zone)
{
	for (int_t iter = 0; iter < _maxIter - 1; ++iter)
	{
		if (marking(zone)) return;
		press(zone);
	}
	marking(zone);
	destroy(zone);
}

bool PressureFix::marking(std::vector<std::shared_ptr<Zone>> zone)
{
	real_t pressure = 0;
	bool stopper = true;
	std::vector<std::shared_ptr<Cell> > cell = zone[0]->getGrid()->getCell();
	real_t u1;
	real_t u2;
	real_t u3;
	real_t sizeX = zone[0]->getGrid()->getSizeX();
	std::shared_ptr<FRbasis> basis = std::make_shared<FRbasis>(zone[0]->getGrid());

	for (int_t icell = 0; icell < _num_cell; ++icell)
	{
		real_t posX = cell[icell]->getPosX();
		for (int_t idegree = 0; idegree <= zone[0]->getPolyOrder(); ++idegree)
		{
			u1 = zone[0]->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone[0]->getPolyOrder()));
			u2 = zone[1]->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone[1]->getPolyOrder()));
			u3 = zone[2]->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone[2]->getPolyOrder()));
			pressure = (GAMMA - 1)*(u3 - 0.5*pow(u2, 2.0) / u1);
			if (pressure > 0.0) _marker[idegree][icell] = true;
			else _marker[idegree][icell] = false;
			stopper = (stopper && _marker[idegree][icell]);
		}
	}

	return stopper;
}

void PressureFix::press(std::vector<std::shared_ptr<Zone>> zone)
{
	std::vector<std::vector<std::vector<real_t>>> temp_DOF;
	temp_DOF.resize(zone.size());
	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		temp_DOF[inum] = zone[inum]->getDOF();
	}	

	for (int_t icell = 0; icell < zone[0]->getGrid()->getNumCell(); ++icell)
	{
		for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		{
			for (int_t inum = 0; inum < zone.size(); ++inum)
				if (_marker[idegree][icell] == false)
					temp_DOF[inum][idegree][icell] = 0.5*(temp_DOF[inum][idegree][icell] + zone[inum]->getAverage(icell));
		}
	}

	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		zone[inum]->setDOF(temp_DOF[inum]);
		zone[inum]->calSolution();
	}
}

void PressureFix::destroy(std::vector<std::shared_ptr<Zone>> zone)
{
	std::vector<std::vector<std::vector<real_t>>> temp_DOF;
	temp_DOF.resize(zone.size());
	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		temp_DOF[inum] = zone[inum]->getDOF();
	}

	for (int_t icell = 0; icell < zone[0]->getGrid()->getNumCell(); ++icell)
	{
		for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		{
			for (int_t inum = 0; inum < zone.size(); ++inum)
				if (_marker[idegree][icell] == false)
					temp_DOF[inum][idegree][icell] = zone[inum]->getAverage(icell);
		}
	}

	for (int_t inum = 0; inum < zone.size(); ++inum)
	{
		zone[inum]->setDOF(temp_DOF[inum]);
		zone[inum]->calSolution();
	}
}