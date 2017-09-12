#pragma once
#include "DataType.h"
#include "Zone.h"
#include "Quadrature.h"

#define CONST_K 1

class Limiter
{
public:
	// Constructor / p.m. limiter type, Zone(object), conserved variable to check troubled cell(0: rho, 1: rho*u, 2: E)
	Limiter(Type, std::vector<std::shared_ptr<Zone>>, int_t);

	// Desturctor
	~Limiter();

public:
	// Functions
	// calculate local projection limiter / p.m. Zone(object)
	void hMLP_Limiter(std::vector<std::shared_ptr<Zone>>);

protected:
	// Variables
	Type _limiter;
	int_t _limitVar;
	int_t _polyOrder;
	int_t _num_cell;
	real_t _size_cell;

protected:
	// Functions
	// Projection n degree DOF to n-1 degree / p.m. Zone(object), current degree, troubled-cell marker
	void troubleCellProject(std::shared_ptr<Zone>, std::vector<int_t>&, const std::vector<bool>&);

	// Troubled-cell marker / p.m. Zone(object), current degree, cell number
	bool troubleCellMarker(std::shared_ptr<Zone>, int_t, int_t) const;

	// MLP-based troubled-cell marker(Augmented MLP condition) / p.m. Zone(object), cell number
	bool augmentMLPmarker(std::shared_ptr<Zone>, int_t) const;

	// MLP-based troubled-cell marker(P1 MLP condition) / p.m. Zone(object), cell number
	bool P1MLPmarker(std::shared_ptr<Zone>, int_t) const;

	// hierarchical troubled-cell marker(Smooth extrema detector) / p.m. Zone(object), cell number
	bool extremaDetector(std::shared_ptr<Zone>, int_t) const;

	// Apply MLP-u limiter / p.m. Zone(object), cell number, MLP limiter function / r.t. DOF
	std::vector<std::vector<real_t>> apply_limiter(std::shared_ptr<Zone>, int_t, real_t) const;

	// Compute MLP limiter / p.m. Zone(object), cell number
	real_t MLP_limit_ftn(std::shared_ptr<Zone>, int_t) const;

	// Compute limiting PI / p.m. vertex difference, linear reconstruction
	real_t limit_PI(real_t, real_t) const;

	// Compute MLP-u2(or MLP-Venkatakrishnan) limiter / p.m. vertex difference, linear reconstruction
	real_t MLP_u2(real_t, real_t) const;

	// Compute projection
	// p.m. degree of projected P, Zone to project, cell number, x coordinate
	// r.t. projected value
	real_t projectionTo(int_t, std::shared_ptr<Zone>, int_t, real_t) const;

	// Compute projection
	// p.m. degree of projected P, Zone to project, cell number
	// r.t. projected DOF
	std::vector<std::vector<real_t>> DOFprojectionTo(int_t, std::shared_ptr<Zone>, int_t) const;

	// Compute cell average
	// p.m. Zone, cell number
	// r.t. cell average
	real_t getAverage(std::shared_ptr<Zone>, int_t) const;

	// Vandermonde matrix / p.m. vector of solution points
	std::vector<std::vector<real_t>> Vandermonde(std::vector<real_t>) const;

	// Inverse Vandermonde matrix / p.m. vector of solution points
	std::vector<std::vector<real_t>> inv_Vandermonde(std::vector<real_t>) const;
};