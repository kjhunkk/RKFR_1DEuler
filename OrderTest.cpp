#include "OrderTest.h"

OrderTest::OrderTest()
{
	_num = 0;
}

OrderTest::~OrderTest()
{

}

std::vector<real_t> OrderTest::ZoneToArray(std::shared_ptr<Zone> zone)
{
	std::vector<real_t> solution;
	std::vector<std::shared_ptr<Cell> > cell = zone->getGrid()->getCell();

	for(int_t icell = 0; icell < zone->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true) solution.push_back(zone->getDescSolution()[icell]);
	}

	return solution;
}

std::vector<real_t> OrderTest::ZoneToPoly(std::shared_ptr<Zone> zone)
{
	std::vector<real_t> solution;
	std::vector<std::shared_ptr<Cell> > cell = zone->getGrid()->getCell();
	real_t sizeX = zone->getGrid()->getSizeX();

	for (int_t icell = 0; icell < zone->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			real_t posX = cell[icell]->getPosX();
			for (int_t idegree = 0; idegree < QuadDegree; ++idegree)
			{
				solution.push_back(zone->getPolySolution(icell, posX + 0.5*sizeX*Gauss3_X(idegree)));
			}
		}
	}

	return solution;
}

real_t OrderTest::L1error(std::vector<real_t> computed)
{
	int_t num = _exact.size();
	if (num != computed.size()) ERROR("different number of solutions");

	real_t L1 = 0;
	for (int_t icell = 0; icell < num; ++icell)
	{
		L1 += abs(computed[icell] - _exact[icell]);
	}

	L1 /= double(num);

	return L1;
}

real_t OrderTest::L2error(std::vector<real_t> computed)
{
	if (_num != computed.size()) ERROR("different number of solutions");

	real_t L2 = 0;
	for (int_t icell = 0; icell < _num; ++icell)
	{
		L2 += pow(computed[icell] - _exact[icell], 2.0);
	}
	L2 /= double(_num);
	L2 = sqrt(L2);

	return L2;
}

real_t OrderTest::Linf_error(std::vector<real_t> computed)
{
	if (_num != computed.size()) ERROR("different number of solutions");

	real_t Linf = 0;
	real_t diff = 0;
	for (int_t icell = 0; icell < _num; ++icell)
	{
		Linf =std::max(Linf, abs(computed[icell] - _exact[icell]));
	}

	return Linf;
}