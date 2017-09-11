#pragma once
#include "DataType.h"

class Reader
{
public:
	// Constructor
	Reader();

	// Destructor
	~Reader();

public:
	// Fuctions
	inline Type getFluxScheme() const { return _fluxScheme; }

	inline Type getCorrectFunc() const { return _correctFunc; }

	inline Type getLimiter() const { return _limiter; }

	inline Type getInitial() const { return _initial; }

	inline Type getBoundary() const { return _boundary; }

	inline Type getTimeInteg() const { return _timeInteg; }

	inline int_t getLimitVar() const { return _limitVar; }

	inline int_t getPolyOrder() const { return _polyOrder; }

	inline real_t getStartPoint() const { return _startPoint; }

	inline real_t getArea() const { return _area; }

	inline real_t getSizeX() const { return _sizeX; }

	inline real_t getCFL() const { return _CFL; }

	inline real_t getTargetT() const { return _T; }

	// Read file / p.m. file name / r.t. true/false
	bool readFile(std::string);

protected:
	// Variables
	Type _fluxScheme;
	Type _correctFunc;
	Type _limiter;
	Type _initial;
	Type _boundary;
	Type _timeInteg;
	int_t _limitVar; // converved variable for limiting
	int_t _polyOrder;
	real_t _startPoint;
	real_t _area;
	real_t _sizeX;
	real_t _CFL;
	real_t _T;

protected:
	// Functions
	void print() const;
};