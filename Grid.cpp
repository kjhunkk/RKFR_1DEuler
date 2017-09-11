#include "Grid.h"

Grid::Grid(real_t start, real_t area, real_t sizeX)
{
	_area = area;
	_sizeX = sizeX;
	_num_cell = area / sizeX + epsilon;
	_num_cell += 2 * GHOST;
	_cell.resize(_num_cell);
	real_t coord_X = 0.0;
	for (int_t icell = 0; icell < GHOST; ++icell)
	{
		_cell[icell] = std::make_shared<Cell>();
		coord_X = start + sizeX*(0.5 + (icell - GHOST));
		_cell[icell]->initialize(false, sizeX, coord_X);
	}
	for (int_t icell = GHOST; icell < _num_cell - GHOST; ++icell)
	{
		_cell[icell] = std::make_shared<Cell>();
		coord_X = start + sizeX*(0.5 + (icell - GHOST));
		_cell[icell]->initialize(true, sizeX, coord_X);
	}
	for (int_t icell = _num_cell - GHOST; icell < _num_cell; ++icell)
	{
		_cell[icell] = std::make_shared<Cell>();
		coord_X = start + sizeX*(0.5 + (icell - GHOST));
		_cell[icell]->initialize(false, sizeX, coord_X);
	}
}

Grid::~Grid()
{

}