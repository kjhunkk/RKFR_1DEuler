#include "Post.h"

Post::Post(std::shared_ptr<Reader> reader)
{
	_reader = reader;
}

Post::~Post()
{

}

void Post::cons_solution(std::shared_ptr<Zone> zone1, std::shared_ptr<Zone> zone2, std::shared_ptr<Zone> zone3) const
{
	// Determine file name
	std::string fileName = "./output/";
	fileName += "cons_RKFR_1D_P";
	fileName += zone1->getPolyOrder();
	fileName += "_";
	fileName += _reader->getInitial();
	fileName += "_";
	fileName += _reader->getCorrectFunc();
	fileName += "_CFL=";
	fileName += std::to_string(_reader->getCFL());
	fileName += ".plt";

	// Cell size
	real_t dx = zone1->getGrid()->getSizeX();

	// Build solution arrays to post
	std::vector<std::shared_ptr<Cell> > cell = zone1->getGrid()->getCell();
	std::vector<real_t> solution1 = zone1->getDescSolution();
	std::vector<real_t> solution2 = zone2->getDescSolution();
	std::vector<real_t> solution3 = zone3->getDescSolution();
	std::vector<real_t> X;
	std::vector<real_t> U1, U2, U3;
	for (int_t icell = 0; icell < zone1->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			X.push_back(cell[icell]->getPosX());
			U1.push_back(solution1[icell]);
			U2.push_back(solution2[icell]);
			U3.push_back(solution3[icell]);
		}
	}

	// Write solutions
	write(fileName, "X", "Rho", "Rho*U", "E", X, U1, U2, U3);
}

void Post::cons_solution(const std::string& name, std::shared_ptr<Zone> zone1, std::shared_ptr<Zone> zone2, std::shared_ptr<Zone> zone3) const
{
	// Determine file name
	std::string fileName = "./output/";
	fileName += "cons_RKFR_1D_P";
	fileName += std::to_string(zone1->getPolyOrder());
	fileName += "_";
	fileName += _reader->getInitial();
	fileName += "_";
	fileName += _reader->getCorrectFunc();
	fileName += "_CFL=";
	fileName += std::to_string(_reader->getCFL());
	fileName += "_";
	fileName += name;
	fileName += ".plt";

	// Cell size
	real_t dx = zone1->getGrid()->getSizeX();

	// Build solution arrays to post
	std::vector<std::shared_ptr<Cell> > cell = zone1->getGrid()->getCell();
	std::vector<real_t> solution1 = zone1->getDescSolution();
	std::vector<real_t> solution2 = zone2->getDescSolution();
	std::vector<real_t> solution3 = zone3->getDescSolution();
	std::vector<real_t> X;
	std::vector<real_t> U1, U2, U3;
	for (int_t icell = 0; icell < zone1->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			X.push_back(cell[icell]->getPosX());
			U1.push_back(solution1[icell]);
			U2.push_back(solution2[icell]);
			U3.push_back(solution3[icell]);
		}
	}

	// Write solutions
	write(fileName, "X", "Rho", "Rho*U", "E", X, U1, U2, U3);
}

void Post::cons_FRsolution(const std::string& name, std::shared_ptr<Zone> zone1, std::shared_ptr<Zone> zone2, std::shared_ptr<Zone> zone3) const
{
	// Determine file name
	std::string fileName = "./output/";
	fileName += "cons_RKFR_1D_P";
	fileName += std::to_string(zone1->getPolyOrder());
	fileName += "_FRsolution_";
	fileName += _reader->getInitial();
	fileName += "_";
	fileName += _reader->getCorrectFunc();
	fileName += "_CFL=";
	fileName += std::to_string(_reader->getCFL());
	fileName += "_";
	fileName += name;
	fileName += ".plt";

	// Build solution arrays to post
	std::vector<std::shared_ptr<Cell> > cell = zone1->getGrid()->getCell();
	std::vector<real_t> X;
	std::vector<real_t> U1;
	std::vector<real_t> U2;
	std::vector<real_t> U3;
	real_t sizeX = zone1->getGrid()->getSizeX();
	real_t dx = zone1->getGrid()->getSizeX() / double(POST_GRID_NUM);
	std::shared_ptr<FRbasis> basis = std::make_shared<FRbasis>(zone1->getGrid());

	for (int_t icell = 0; icell < zone1->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			real_t posX = cell[icell]->getPosX();
			for (int_t idegree = 0; idegree <= zone1->getPolyOrder(); ++idegree)
			{
				X.push_back(posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone1->getPolyOrder()));
				U1.push_back(zone1->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone1->getPolyOrder())));
				U2.push_back(zone2->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone2->getPolyOrder())));
				U3.push_back(zone3->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone3->getPolyOrder())));
			}
		}
	}

	// Write solutions
	write(fileName, "X", "Rho", "Rho*U", "E", X, U1, U2, U3);
}

void Post::solution(std::shared_ptr<Zone> zone1, std::shared_ptr<Zone> zone2, std::shared_ptr<Zone> zone3) const
{
	// Determine file name
	std::string fileName = "./output/";
	fileName += "RKFR_1D_P";
	fileName += zone1->getPolyOrder();
	fileName += "_";
	fileName += _reader->getInitial();
	fileName += "_";
	fileName += _reader->getCorrectFunc();
	fileName += "_CFL=";
	fileName += std::to_string(_reader->getCFL());
	fileName += ".plt";

	// Cell size
	real_t dx = zone1->getGrid()->getSizeX();

	// Build solution arrays to post
	std::vector<std::shared_ptr<Cell> > cell = zone1->getGrid()->getCell();
	std::vector<real_t> solution1 = zone1->getDescSolution();
	std::vector<real_t> solution2 = zone2->getDescSolution();
	std::vector<real_t> solution3 = zone3->getDescSolution();
	std::vector<real_t> X;
	std::vector<real_t> U1, U2, U3;
	for (int_t icell = 0; icell < zone1->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			X.push_back(cell[icell]->getPosX());
			U1.push_back(solution1[icell]);
			U2.push_back(solution2[icell] / solution1[icell]);
			U3.push_back((GAMMA - 1)*(solution3[icell] - 0.5*pow(solution2[icell], 2.0) / solution1[icell]));
		}
	}

	// Write solutions
	write(fileName, "X", "Density", "Velocity", "Pressure", X, U1, U2, U3);
}

void Post::solution(const std::string& name, std::shared_ptr<Zone> zone1, std::shared_ptr<Zone> zone2, std::shared_ptr<Zone> zone3) const
{
	// Determine file name
	std::string fileName = "./output/";
	fileName += "RKFR_1D_P";
	fileName += std::to_string(zone1->getPolyOrder());
	fileName += "_";
	fileName += _reader->getInitial();
	fileName += "_";
	fileName += _reader->getCorrectFunc();
	fileName += "_CFL=";
	fileName += std::to_string(_reader->getCFL());
	fileName += "_";
	fileName += name;
	fileName += ".plt";

	// Cell size
	real_t dx = zone1->getGrid()->getSizeX();

	// Build solution arrays to post
	std::vector<std::shared_ptr<Cell> > cell = zone1->getGrid()->getCell();
	std::vector<real_t> solution1 = zone1->getDescSolution();
	std::vector<real_t> solution2 = zone2->getDescSolution();
	std::vector<real_t> solution3 = zone3->getDescSolution();
	std::vector<real_t> X;
	std::vector<real_t> U1, U2, U3;
	for (int_t icell = 0; icell < zone1->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			X.push_back(cell[icell]->getPosX());
			U1.push_back(solution1[icell]);
			U2.push_back(solution2[icell] / solution1[icell]);
			U3.push_back((GAMMA - 1)*(solution3[icell] - 0.5*pow(solution2[icell], 2.0) / solution1[icell]));
		}
	}

	// Write solutions
	write(fileName, "X", "Density", "Velocity", "Pressure", X, U1, U2, U3);
}

void Post::FRsolution(const std::string& name, std::shared_ptr<Zone> zone1, std::shared_ptr<Zone> zone2, std::shared_ptr<Zone> zone3) const
{
	// Determine file name
	std::string fileName = "./output/";
	fileName += "RKFR_1D_P";
	fileName += std::to_string(zone1->getPolyOrder());
	fileName += "_FRsolution_";
	fileName += _reader->getInitial();
	fileName += "_";
	fileName += _reader->getCorrectFunc();
	fileName += "_CFL=";
	fileName += std::to_string(_reader->getCFL());
	fileName += "_";
	fileName += name;
	fileName += ".plt";

	// Build solution arrays to post
	std::vector<std::shared_ptr<Cell> > cell = zone1->getGrid()->getCell();
	std::vector<real_t> X;
	std::vector<real_t> U1;
	std::vector<real_t> U2;
	std::vector<real_t> U3;
	real_t u1;
	real_t u2;
	real_t u3;
	real_t sizeX = zone1->getGrid()->getSizeX();
	real_t dx = zone1->getGrid()->getSizeX() / double(POST_GRID_NUM);
	std::shared_ptr<FRbasis> basis = std::make_shared<FRbasis>(zone1->getGrid());

	for (int_t icell = 0; icell < zone1->getGrid()->getNumCell(); ++icell)
	{
		if (cell[icell]->getType() == true)
		{
			real_t posX = cell[icell]->getPosX();
			for (int_t idegree = 0; idegree <= zone1->getPolyOrder(); ++idegree)
			{
				X.push_back(posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone1->getPolyOrder()));
				u1 = zone1->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone1->getPolyOrder()));
				u2 = zone2->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone2->getPolyOrder()));
				u3 = zone3->getPolySolution(icell, posX + 0.5*sizeX*basis->getSolutionPoint(idegree, zone3->getPolyOrder()));
				U1.push_back(u1);
				U2.push_back(u2 / u1);
				U3.push_back((GAMMA - 1)*(u3 - 0.5*pow(u2, 2.0) / u1));
			}
		}
	}

	// Write solutions
	write(fileName, "X", "Density", "Velocity", "Pressure", X, U1, U2, U3);
}

void Post::error(std::shared_ptr<Zone> zone) const
{

}

void Post::write(const std::string& fileName, const std::string& VN1, const std::string& VN2, std::vector<real_t> V1, std::vector<real_t> V2) const
{
	std::ofstream file;
	file.open(fileName, std::ios::trunc);
	int_t _num_element = V1.size();
	if (_num_element != V2.size()) ERROR("Number of elements does not match");

	if (file.is_open())
	{
		MESSAGE("Output file open");
		file << "variables = " << VN1 << ", " << VN2 << "\n";
		file << "zone t = \"RKFR 1D\", i=" << _num_element << ", f=point\n";
		for (int_t ielem = 0; ielem < _num_element; ++ielem)
		{
			file << std::to_string(V1[ielem]) << "\t" << std::to_string(V2[ielem]) << "\n";
		}
		file.close();
	}
	else ERROR("cannot open output file");
}

void Post::write(const std::string& fileName, const std::string& VN1, const std::string& VN2, const std::string& VN3, const std::string& VN4, std::vector<real_t> V1, std::vector<real_t> V2, std::vector<real_t> V3, std::vector<real_t> V4) const
{
	std::ofstream file;
	file.open(fileName, std::ios::trunc);
	int_t _num_element = V1.size();
	if ((_num_element != V2.size()) || (_num_element != V3.size()) || (_num_element != V4.size())) ERROR("Number of elements does not match");

	if (file.is_open())
	{
		MESSAGE("Output file open");
		file << "variables = " << VN1 << ", " << VN2 << ", " << VN3 << ", " << VN4 << "\n";
		file << "zone t = \"RKFR 1D\", i=" << _num_element << ", f=point\n";
		for (int_t ielem = 0; ielem < _num_element; ++ielem)
		{
			file << std::to_string(V1[ielem]) << "\t" << std::to_string(V2[ielem]) << "\t" << std::to_string(V3[ielem]) << "\t" << std::to_string(V4[ielem]) << "\n";
		}
		file.close();
	}
	else ERROR("cannot open output file");
}