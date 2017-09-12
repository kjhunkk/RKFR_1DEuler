#include "Limiter.h"

Limiter::Limiter(Type limiter, std::vector<std::shared_ptr<Zone>> zone, int_t limitVar)
{
	_limiter = limiter;
	_limitVar = limitVar;
	_polyOrder = zone[0]->getPolyOrder();
	_num_cell = zone[0]->getGrid()->getNumCell();
	_size_cell = zone[0]->getGrid()->getSizeX();
}

Limiter::~Limiter()
{

}

void Limiter::hMLP_Limiter(std::vector<std::shared_ptr<Zone>> zone)
{
	// No limiter if PO
	if ((_polyOrder == 0) || (_limiter == "none")) return;

	// Temporary variables for hMLP
	// Current projected degree
	std::vector<int_t> projectDegree(_num_cell, _polyOrder);

	// hMLP limiting process
	for (int_t step = 0; step < _polyOrder; ++step)
	{
		// Troubled-cell marker
		std::vector<bool> marker(_num_cell, true);

		// Marking troubled-cell
		for (int_t icell = GHOST; icell < _num_cell - GHOST; ++icell)
		{
			marker[icell] = troubleCellMarker(zone[0], projectDegree[icell], icell);
		}

		// Project troubled-cell
		for(int_t inum = 0; inum < zone.size(); ++inum)
			troubleCellProject(zone[inum], projectDegree, marker);
	}

	// Update zone
	for (int_t inum = 0; inum < zone.size(); ++inum)
		zone[inum]->calSolution();
}

void Limiter::troubleCellProject
(std::shared_ptr<Zone> zone, std::vector<int_t>& degree, const std::vector<bool>& marker)
{
	// Temporary DOF and Zone
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();

	// Projection
	for (int_t icell = 0; icell < _num_cell; ++icell)
	{
		if (marker[icell] == false)
		{
			if (degree[icell] >= 2)
				temp_DOF = DOFprojectionTo(--degree[icell], zone, icell);
			else if (degree[icell] == 1)
			{
				std::vector<std::vector<real_t>> tester = apply_limiter(zone, icell, MLP_limit_ftn(zone, icell));
				temp_DOF = apply_limiter(zone, icell, MLP_limit_ftn(zone, icell));
			}
			else ERROR("step degree");
			zone->setDOF(temp_DOF);
		}
	}

	// Update zone
	zone->setDOF(temp_DOF);
	zone->calSolution();
}

bool Limiter::troubleCellMarker(std::shared_ptr<Zone> zone, int_t degree, int_t icell) const
{
	//bool marker = augmentMLPmarker(zone, icell);
	bool marker = P1MLPmarker(zone, icell);
	if ((marker == false) && (degree > 1)) marker = extremaDetector(zone, icell);
	
	return marker;
}

bool Limiter::augmentMLPmarker(std::shared_ptr<Zone> zone, int_t icell) const
{
	// Augmented MLP condition marker
	bool marker = true;

	// Temporary DOF
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();

	// Variables
	real_t coord_x_left = zone->getGrid()->getCell()[icell]->getPosX() - 0.5*_size_cell;
	real_t coord_x_right = zone->getGrid()->getCell()[icell]->getPosX() + 0.5*_size_cell;

	// Variables to MLP condition
	real_t max_avgQ; /// maximum averaged Q
	real_t min_avgQ; /// minimum averaged Q
	real_t max_appQ; /// maximum approximated Q
	real_t min_appQ; /// minimum approximated Q
	real_t appQ; /// approximated vertex Q

	// Left cell MLP condition
	appQ = zone->getPolySolution(icell, coord_x_left);
	max_appQ = std::max(zone->getPolySolution(icell - 1, coord_x_left), zone->getPolySolution(icell, coord_x_left));
	min_appQ = std::min(zone->getPolySolution(icell - 1, coord_x_left), zone->getPolySolution(icell, coord_x_left));
	max_avgQ = std::max(getAverage(zone, icell - 1), getAverage(zone, icell));
	min_avgQ = std::min(getAverage(zone, icell - 1), getAverage(zone, icell));
	if (!((max_avgQ - max_appQ > -epsilon) && (min_appQ - min_avgQ > -epsilon)))
		marker = false;

	// Right cell MLP condition
	appQ = zone->getPolySolution(icell, coord_x_right);
	max_appQ = std::max(zone->getPolySolution(icell, coord_x_right), zone->getPolySolution(icell + 1, coord_x_right));
	min_appQ = std::min(zone->getPolySolution(icell, coord_x_right), zone->getPolySolution(icell + 1, coord_x_right));
	max_avgQ = std::max(getAverage(zone, icell), getAverage(zone, icell + 1));
	min_avgQ = std::min(getAverage(zone, icell), getAverage(zone, icell + 1));
	if (!((max_avgQ - max_appQ > -epsilon) && (min_appQ - min_avgQ > -epsilon)))
		marker = false;

	return marker;
}

bool Limiter::P1MLPmarker(std::shared_ptr<Zone> zone, int_t icell) const
{
	// P1 MLP condition marker
	bool marker = true;

	// Temporary DOF
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();

	// Variables
	real_t coord_x_left = zone->getGrid()->getCell()[icell]->getPosX() - 0.5*_size_cell;
	real_t coord_x_right = zone->getGrid()->getCell()[icell]->getPosX() + 0.5*_size_cell;

	// Variables to MLP condition
	real_t max_avgQ; /// maximum averaged Q
	real_t min_avgQ; /// minimum averaged Q
	real_t P1_projQ; /// P1 projected Q
	real_t appQ; /// approximated vertex Q

	// Left cell MLP condition
	appQ = zone->getPolySolution(icell, coord_x_left);
	P1_projQ = projectionTo(1, zone, icell, coord_x_left);
	max_avgQ = std::max(getAverage(zone, icell - 1), getAverage(zone, icell));
	min_avgQ = std::min(getAverage(zone, icell - 1), getAverage(zone, icell));
	if (!((max_avgQ - P1_projQ > -epsilon) && (P1_projQ - min_avgQ > -epsilon)))
		marker = false;

	// Right cell MLP condition
	appQ = zone->getPolySolution(icell, coord_x_right);
	P1_projQ = projectionTo(1, zone, icell, coord_x_right);
	max_avgQ = std::max(getAverage(zone, icell), getAverage(zone, icell + 1));
	min_avgQ = std::min(getAverage(zone, icell), getAverage(zone, icell + 1));
	if (!((max_avgQ - P1_projQ > -epsilon) && (P1_projQ - min_avgQ > -epsilon)))
		marker = false;

	return marker;
}

bool Limiter::extremaDetector(std::shared_ptr<Zone> zone, int_t icell) const
{
	// Decomposing the DG-Pn approximation
	real_t avgQ = getAverage(zone, icell);
	real_t max_avgQ;
	real_t min_avgQ;
	real_t P1_projected;
	real_t Pn_projected_slope;
	real_t P1_filtered_Pn;
	
	real_t coord_x_left = zone->getGrid()->getCell()[icell]->getPosX() - 0.5*_size_cell;
	real_t coord_x_right = zone->getGrid()->getCell()[icell]->getPosX() + 0.5*_size_cell;

	real_t leftQ = zone->getPolySolution(icell, coord_x_left);
	real_t rightQ = zone->getPolySolution(icell, coord_x_right);

	// Deactivation threshold
	real_t threshold = std::max(0.001*avgQ, _size_cell);
	if ((abs(leftQ - avgQ) <= threshold) && (abs(rightQ - avgQ) <= threshold))
		return true; /// deactivation

	// Left vertex
	bool leftMarker = false;
	P1_projected = projectionTo(1, zone, icell, coord_x_left);
	Pn_projected_slope = P1_projected - avgQ;
	P1_filtered_Pn = zone->getPolySolution(icell, coord_x_left) - P1_projected;
	max_avgQ = std::max(getAverage(zone, icell - 1), getAverage(zone, icell));
	min_avgQ = std::min(getAverage(zone, icell - 1), getAverage(zone, icell));
	
	// Marking
	if ((Pn_projected_slope > 0.0) && (P1_filtered_Pn < 0.0) && (leftQ > min_avgQ)) leftMarker = true;
	if ((Pn_projected_slope < 0.0) && (P1_filtered_Pn > 0.0) && (leftQ < max_avgQ)) leftMarker = true;

	// Right vertex
	bool rightMarker = false;
	P1_projected = projectionTo(1, zone, icell, coord_x_right);
	Pn_projected_slope = P1_projected - avgQ;
	P1_filtered_Pn = zone->getPolySolution(icell, coord_x_right) - P1_projected;
	max_avgQ = std::max(getAverage(zone, icell), getAverage(zone, icell + 1));
	min_avgQ = std::min(getAverage(zone, icell), getAverage(zone, icell + 1));

	// Marking
	if ((Pn_projected_slope > 0.0) && (P1_filtered_Pn < 0.0) && (rightQ > min_avgQ)) rightMarker = true;
	if ((Pn_projected_slope < 0.0) && (P1_filtered_Pn > 0.0) && (rightQ < max_avgQ)) rightMarker = true;

	// Smooth extrema detect
	return (leftMarker && rightMarker);
}

std::vector<std::vector<real_t>> Limiter::apply_limiter(std::shared_ptr<Zone> zone, int_t num_cell, real_t MLP_ftn) const
{
	// Temporary object
	std::vector<std::vector<real_t>> nodal_DOF = zone->getDOF();

	// Vandermonde matrix
	std::vector<std::vector<real_t>> Vmatrix, inv_Vmatrix;
	std::vector<real_t> lobatto_points(_polyOrder + 1, 0.0);
	switch (_polyOrder)
	{
	case 1:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto2_X(iorder);
		break;
	case 2:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto3_X(iorder);
		break;
	case 3:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto4_X(iorder);
		break;
	case 4:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto5_X(iorder);
		break;
	default: ERROR("exceeds maximum order");
	}
	Vmatrix = Vandermonde(lobatto_points);
	inv_Vmatrix = inv_Vandermonde(lobatto_points);

	// Convert nodal to modal expression
	std::vector<real_t> modal_DOF(_polyOrder + 1, 0.0);
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		for (int_t jdegree = 0; jdegree <= _polyOrder; ++jdegree)
			modal_DOF[idegree] += inv_Vmatrix[idegree][jdegree] * nodal_DOF[jdegree][num_cell];

	// Apply MLP limiter
	modal_DOF[1] = MLP_ftn*modal_DOF[1];

	// Delete ~degree 1
	for (int_t idegree = _polyOrder; idegree > 1; --idegree)
		modal_DOF[idegree] = 0.0;

	// Convert modal to nodal expression
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
	{
		nodal_DOF[idegree][num_cell] = 0.0;
		for (int_t jdegree = 0; jdegree <= _polyOrder; ++jdegree)
			nodal_DOF[idegree][num_cell] += Vmatrix[idegree][jdegree] * modal_DOF[jdegree];
	}

	return nodal_DOF;
}

real_t Limiter::MLP_limit_ftn(std::shared_ptr<Zone> zone, int_t icell) const
{
	// Variables
	real_t limit_ftn_left;
	real_t limit_ftn_right;
	real_t coord_x_left = zone->getGrid()->getCell()[icell]->getPosX() - 0.5*_size_cell;
	real_t coord_x_right = zone->getGrid()->getCell()[icell]->getPosX() + 0.5*_size_cell;
	real_t avgQ = getAverage(zone, icell);
	real_t del_m = projectionTo(1, zone, icell, coord_x_right) - avgQ;

	// Compute MLP function
	if (del_m > epsilon)
	{
		limit_ftn_right = limit_PI(std::max(getAverage(zone, icell), getAverage(zone, icell + 1)) - avgQ, del_m);
		limit_ftn_left = limit_PI(std::min(getAverage(zone, icell - 1), getAverage(zone, icell)) - avgQ, -del_m);
	}
	else if (del_m < -epsilon)
	{
		limit_ftn_right = limit_PI(std::min(getAverage(zone, icell), getAverage(zone, icell + 1)) - avgQ, del_m);
		limit_ftn_left = limit_PI(std::max(getAverage(zone, icell - 1), getAverage(zone, icell)) - avgQ, -del_m);
	}
	else return 1.0;

	return std::min(limit_ftn_right, limit_ftn_left);
}

real_t Limiter::limit_PI(real_t del_p, real_t del_m) const
{
	if (_limiter == "MLP-u1") return std::min(1.0, del_p / del_m);
	if (_limiter == "MLP-u2") return MLP_u2(del_p, del_m);
	ERROR("cannot find limiter");
	return 0;
}

real_t Limiter::MLP_u2(real_t del_p, real_t del_m) const
{
	real_t ep = pow(CONST_K*_size_cell, 1.5);
	real_t num = (pow(del_p, 2.0) + pow(ep, 2.0))*del_m + 2.0*pow(del_m, 2.0)*del_p;
	real_t den = del_m*(pow(del_p, 2.0) + 2.0*pow(del_m, 2.0) + del_m*del_p + pow(ep, 2.0));

	return num / den;
}

real_t Limiter::projectionTo(int_t degree, std::shared_ptr<Zone> zone, int_t icell, real_t coord_x) const
{
	// Temporary object
	std::shared_ptr<Zone> temp_zone = std::make_shared<Zone>(*zone);
	std::vector<std::vector<real_t> > nodal_DOF = zone->getDOF();

	// Vandermonde matrix
	std::vector<std::vector<real_t>> Vmatrix, inv_Vmatrix;
	std::vector<real_t> lobatto_points(_polyOrder + 1, 0.0);
	switch (_polyOrder)
	{
	case 1:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto2_X(iorder);
		break;
	case 2:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto3_X(iorder);
		break;
	case 3:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto4_X(iorder);
		break;
	case 4:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto5_X(iorder);
		break;
	default: ERROR("exceeds maximum order");
	}
	Vmatrix = Vandermonde(lobatto_points);
	inv_Vmatrix = inv_Vandermonde(lobatto_points);

	// Convert nodal to modal expression
	std::vector<real_t> modal_DOF(_polyOrder + 1, 0.0);
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		for (int_t jdegree = 0; jdegree <= _polyOrder; ++jdegree)
			modal_DOF[idegree] += inv_Vmatrix[idegree][jdegree] * nodal_DOF[jdegree][icell];
	
	// Projection to n degree
	for (int_t idegree = _polyOrder; idegree > degree; --idegree)
		modal_DOF[idegree] = 0.0;

	// Convert modal to nodal expression
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
	{
		nodal_DOF[idegree][icell] = 0.0;
		for (int_t jdegree = 0; jdegree <= _polyOrder; ++jdegree)
			nodal_DOF[idegree][icell] += Vmatrix[idegree][jdegree] * modal_DOF[jdegree];
	}
	
	temp_zone->setDOF(nodal_DOF);

	return temp_zone->getPolySolution(icell, coord_x);
}

std::vector<std::vector<real_t>> Limiter::DOFprojectionTo(int_t degree, std::shared_ptr<Zone> zone, int_t icell) const
{
	// Temporary object
	std::shared_ptr<Zone> temp_zone = std::make_shared<Zone>(*zone);
	std::vector<std::vector<real_t> > nodal_DOF = zone->getDOF();

	// Vandermonde matrix
	std::vector<std::vector<real_t>> Vmatrix, inv_Vmatrix;
	std::vector<real_t> lobatto_points(_polyOrder + 1, 0.0);
	switch (_polyOrder)
	{
	case 1:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto2_X(iorder);
		break;
	case 2:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto3_X(iorder);
		break;
	case 3:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto4_X(iorder);
		break;
	case 4:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto5_X(iorder);
		break;
	default: ERROR("exceeds maximum order");
	}
	Vmatrix = Vandermonde(lobatto_points);
	inv_Vmatrix = inv_Vandermonde(lobatto_points);

	// Convert nodal to modal expression
	std::vector<real_t> modal_DOF(_polyOrder + 1, 0.0);
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		for (int_t jdegree = 0; jdegree <= _polyOrder; ++jdegree)
			modal_DOF[idegree] += inv_Vmatrix[idegree][jdegree] * nodal_DOF[jdegree][icell];

	// Projection to n degree
	for (int_t idegree = _polyOrder; idegree > degree; --idegree)
		modal_DOF[idegree] = 0.0;

	// Convert modal to nodal expression
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
	{
		nodal_DOF[idegree][icell] = 0.0;
		for (int_t jdegree = 0; jdegree <= _polyOrder; ++jdegree)
			nodal_DOF[idegree][icell] += Vmatrix[idegree][jdegree] * modal_DOF[jdegree];
	}

	return nodal_DOF;
}

real_t Limiter::getAverage(std::shared_ptr<Zone> zone, int_t num_cell) const
{
	// Temporary object
	std::vector<std::vector<real_t>> temp_nodal_DOF = zone->getDOF();
	std::vector<real_t> nodal_DOF(_polyOrder + 1, 0.0);
	for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
		nodal_DOF[iorder] = temp_nodal_DOF[iorder][num_cell];
	temp_nodal_DOF.clear();

	// Vandermonde matrix
	std::vector<std::vector<real_t>> Vmatrix, inv_Vmatrix;
	std::vector<real_t> lobatto_points(_polyOrder + 1, 0.0);
	switch (_polyOrder)
	{
	case 1:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto2_X(iorder);
		break;
	case 2:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto3_X(iorder);
		break;
	case 3:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto4_X(iorder);
		break;
	case 4:
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			lobatto_points[iorder] = Lobatto5_X(iorder);
		break;
	default: ERROR("exceeds maximum order");
	}
	Vmatrix = Vandermonde(lobatto_points);
	inv_Vmatrix = inv_Vandermonde(lobatto_points);

	// Convert nodal to modal expression
	std::vector<real_t> modal_DOF(_polyOrder + 1, 0.0);
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		for (int_t jdegree = 0; jdegree <= _polyOrder; ++jdegree)
			modal_DOF[idegree] += inv_Vmatrix[idegree][jdegree] * nodal_DOF[jdegree];

	// return mode 0
	return modal_DOF[0];
}

std::vector<std::vector<real_t>> Limiter::Vandermonde(std::vector<real_t> r) const
{
	std::vector<std::vector<real_t>> V;
	V.resize(r.size());
	for (int_t ipoint = 0; ipoint < r.size(); ++ipoint)
		V[ipoint].resize(r.size());

	for (int_t iorder = 0; iorder < r.size(); ++iorder)
		for (int_t ipoint = 0; ipoint < r.size(); ++ipoint)
			V[ipoint][iorder] = LegendreP(iorder, r[ipoint]);

	return V;
}

std::vector<std::vector<real_t>> Limiter::inv_Vandermonde(std::vector<real_t> r) const
{
	std::vector<std::vector<real_t>> invV;
	invV.resize(r.size());
	for (int_t ipoint = 0; ipoint < r.size(); ++ipoint)
		invV[ipoint].resize(r.size());

	real_t den1, den2, den3, den4;
	switch (r.size())
	{
	case 2:
		den1 = 1.0/(r[0] - r[1]);
		invV[0] = { -r[1] * den1, r[0] * den1 };
		invV[1] = { den1, -den1 };
		break;
	case 3:
		den1 = 1.0 / (r[0] * r[1] + r[0] * r[2] - r[1] * r[2] - pow(r[0], 2.0));
		den2 = 1.0 / (r[0] * r[1] - r[0] * r[2] + r[1] * r[2] - pow(r[1], 2.0));
		den3 = 1.0 / (r[0] * r[1] - r[0] * r[2] - r[1] * r[2] + pow(r[2], 2.0));
		invV[0] = { -CONST13*(3.0*r[1] * r[2] + 1.0)*den1, -CONST13*(3.0*r[0] * r[2] + 1.0)*den2, CONST13*(3.0*r[0] * r[1] + 1.0)*den3 };
		invV[1] = { (r[1] + r[2])*den1, (r[0] + r[2])*den2, -(r[0] + r[1])*den3 };
		invV[2] = { -CONST23*den1, -CONST23*den2, CONST23*den3 };
		break;
	case 4:
		den1 = 1.0 / (pow(r[0], 2.0)*(r[1] + r[2] + r[3] - r[0]) - (r[3] * r[0] * r[1]) - (r[3] * r[0] * r[2]) + (r[3] * r[1] * r[2]) - (r[0] * r[1] * r[2]));
		den2 = 1.0 / (pow(r[1], 2.0)*(r[0] + r[2] + r[3] - r[1]) - (r[3] * r[0] * r[1]) + (r[3] * r[0] * r[2]) - (r[3] * r[1] * r[2]) - (r[0] * r[1] * r[2]));
		den3 = 1.0 / (pow(r[2], 2.0)*(r[0] + r[1] + r[3] - r[2]) + (r[3] * r[0] * r[1]) - (r[3] * r[0] * r[2]) - (r[3] * r[1] * r[2]) - (r[0] * r[1] * r[2]));
		den4 = 1.0 / (pow(r[3], 2.0)*(r[0] + r[1] + r[2] - r[3]) - (r[3] * r[0] * r[1]) - (r[3] * r[0] * r[2]) - (r[3] * r[1] * r[2]) + (r[0] * r[1] * r[2]));
		invV[0] = { CONST13*(r[3] + r[1] + r[2] + 3.0*r[3] * r[1] * r[2])*den1, CONST13*(r[3] + r[0] + r[2] + 3.0*r[3] * r[0] * r[2])*den2, CONST13*(r[3] + r[0] + r[1] + 3.0*r[3] * r[0] * r[1])*den3, CONST13*(r[0] + r[1] + r[2] + 3.0*r[0] * r[1] * r[2])*den4 };
		invV[1] = { -0.2*(5.0*r[3] * r[1] + 5.0*r[3] * r[2] + 5.0*r[1] * r[2] + 3.0)*den1, -0.2*(5.0*r[3] * r[0] + 5.0*r[3] * r[2] + 5.0*r[0] * r[2] + 3.0)*den2, -0.2*(5.0*r[3] * r[0] + 5.0*r[3] * r[1] + 5.0*r[0] * r[1] + 3.0)*den3, -0.2*(5.0*r[0] * r[1] + 5.0*r[0] * r[2] + 5.0*r[1] * r[2] + 3.0)*den4 };
		invV[2] = { CONST23*(r[3] + r[1] + r[2])*den1, CONST23*(r[3] + r[0] + r[2])*den2, CONST23*(r[3] + r[0] + r[1])*den3, CONST23*(r[0] + r[1] + r[2])*den4 };
		invV[3] = { -0.4*den1, -0.4*den2, -0.4*den3, -0.4*den4 };
		break;
	default: ERROR("exceeds maximum order");
	}
	return invV;
}