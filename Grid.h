#pragma once
#include "DataType.h"
#include "Cell.h"

class Grid
{
public:
	// Constructor / p.m. start point, area, size dX
	Grid(real_t, real_t, real_t);

	// Destructor
	~Grid();

public:
	// Functions
	inline int_t getNumCell() const { return _num_cell; }

	inline real_t getArea() const { return _area; }

	inline real_t getSizeX() const { return _sizeX; }

	inline std::vector<std::shared_ptr<Cell> > getCell() const { return _cell; }

	inline void setCell(std::vector<std::shared_ptr<Cell>> cell) { _cell = cell; }

protected:
	// Variables
	int_t _num_cell;
	real_t _area;
	real_t _sizeX;
	std::vector<std::shared_ptr<Cell> > _cell;
};