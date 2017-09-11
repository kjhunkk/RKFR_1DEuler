#pragma once
#include "DataType.h"

class Cell
{
public:
	Cell();
	~Cell();

public:
	// Functions
	inline bool getType() const { return _type; }

	inline real_t getSizeX() const { return _sizeX; }

	inline real_t getPosX() const { return _x; }

	// set cell type and size of cell / p.m. type, sizeX, X coordinate
	void initialize(bool, real_t, real_t);

protected:
	// Variables
	// cell type / real or ghost
	bool _type;
	real_t _sizeX;
	real_t _x;
};