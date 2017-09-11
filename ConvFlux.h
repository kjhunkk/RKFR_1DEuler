#pragma once
#include "DataType.h"
#include "Zone.h"

class ConvFlux
{
public:
	// Constructor / p.m. correction function type
	ConvFlux(int_t, Type);

	// Destructor
	virtual ~ConvFlux();

public:
	// Functions
	// Compute flux / p.m. begin, end
	virtual std::vector<real_t> computeFlux(std::vector<real_t>, std::vector<real_t>) const = 0;

	// Compute correction function / p.m. Right(true)/Left(false), computational coordinate
	real_t correctFunc(bool, real_t) const;

	// Compute derivative of correction function / p.m. Right(true)/Left(false), computational coordinate
	real_t DcorrectFunc(bool, real_t) const;

protected:
	// Variables
	Type _correctFunc;
	int_t _polyOrder;

protected:
	// Functions
	// Left boundary DG correction function / p.m. computational coordinate
	real_t correctDG(bool, real_t) const;

	// Left boundary Derivative of DG correction function / p.m. computational coordinate
	real_t DcorrectDG(bool, real_t) const;

	// Left boundary SG correction function / p.m. computational coordinate
	real_t correctSG(bool, real_t) const;

	// Left boundary Derivative of SG correction function / p.m. computational coordinate
	real_t DcorrectSG(bool, real_t) const;

	// Left boundary Lumping for Lobatto points correction function / p.m. computational coordinate
	real_t correct2(bool, real_t) const;

	// Left boundary Derivative of lumping for Lobatto points correction function / p.m. computational coordinate
	real_t Dcorrect2(bool, real_t) const;

	// Left boundary Legendre-Lobatto correction function / p.m. computational coordinate
	real_t correctLo(bool, real_t) const;

	// Left boundary Derivative of Legendre-Lobatto correction function / p.m. computational coordinate
	real_t DcorrectLo(bool, real_t) const;

	// Left boundary Gauss correction function / p.m. computational coordinate
	real_t correctGa(bool, real_t) const;

	// Left boundary Derivative of Gauss correction function / p.m. computational coordinate
	real_t DcorrectGa(bool, real_t) const;
};