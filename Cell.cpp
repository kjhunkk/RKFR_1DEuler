#include "Cell.h"

Cell::Cell()
{
	_type = 1;
	_sizeX = 0.0;
	_x = 0.0;
}

Cell::~Cell()
{

}

void Cell::initialize(bool type, real_t sizeX, real_t x)
{
	_type = type;
	_sizeX = sizeX;
	_x = x;
}