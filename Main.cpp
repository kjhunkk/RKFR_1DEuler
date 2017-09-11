#include "DataType.h"
#include "Reader.h"
#include "Post.h"
#include "Grid.h"
#include "Zone.h"
#include "Boundary.h"
#include "OrderTest.h"
#include "TimeInteg.h"
#include "TimeIntegEuler.h"
#include "TimeIntegRK.h"
#include "ConvFluxSW.h"

// Modified 2017-09-01
// by Juhyeon Kim
// Caution : polynomial order starts from 0

int main()
{
	// Read input file
	std::shared_ptr<Reader> reader = std::make_shared<Reader>();
	if (!reader->readFile("./input.inp")) return 0;

	// Initializing objects
	std::shared_ptr<Post> post = std::make_shared<Post>(reader);
	std::shared_ptr<Grid> grid = std::make_shared<Grid>(reader->getStartPoint(), reader->getArea(), reader->getSizeX());
	std::vector<std::shared_ptr<Zone>> conserv;
	conserv.resize(3);
	for (int_t i = 0; i < 3; ++i)
		conserv[i] = std::make_shared<Zone>(grid, reader->getPolyOrder(), i);
	std::shared_ptr<OrderTest> orderTest = std::make_shared<OrderTest>();

	// Initializing solution domain
	for (int_t i = 0; i < 3; ++i)
		conserv[i]->initialize(std::make_shared<InitialCondition>(reader->getInitial()));

	// initializing boundary condition
	std::vector<std::shared_ptr<Boundary>> bdry;
	bdry.resize(conserv.size());
	for (int_t i = 0; i < conserv.size(); ++i)
		bdry[i] = std::make_shared<Boundary>(reader->getBoundary(), conserv[i]);

	// Save exact solution for order test / using density
	std::vector<real_t> exact = orderTest->ZoneToPoly(conserv[0]);
	orderTest->setExact(exact);

	// Post initial condition
	post->solution("initial", conserv[0], conserv[1], conserv[2]);
	if (reader->getPolyOrder() > 0) post->FRsolution("initial", conserv[0], conserv[1], conserv[2]);

	// Initialzing computational flux
	std::shared_ptr<ConvFlux> convFlux;
	if (reader->getFluxScheme().compare("StegerWarming") == 0) convFlux = std::make_shared<ConvFluxSW>(conserv[0]->getPolyOrder(), reader->getCorrectFunc());
	else ERROR("cannot find computational flux");
	
	// Initialzing time integrator
	std::shared_ptr<TimeInteg> timeInteg;

	if (reader->getTimeInteg() == "Euler")
		timeInteg = std::make_shared<TimeIntegEuler>
		(reader->getLimiter(), reader->getLimitVar(), reader->getCFL(), reader->getTargetT(), convFlux, conserv, bdry);
	else if (reader->getTimeInteg() == "RK3")
		timeInteg = std::make_shared<TimeIntegRK>
		(reader->getLimiter(), reader->getLimitVar(), reader->getCFL(), reader->getTargetT(), convFlux, conserv, bdry, 3);
	else ERROR("cannot find time integrator");

	int_t iter = 0;
	// Time marching
	while (timeInteg->march(conserv))
	{
		iter++;
		if (iter % 100 == 0) MESSAGE("Iteration = " + std::to_string(iter));
		if (iter % 10 == 0) MESSAGE("current time = " + std::to_string(timeInteg->getTime()) + "     time step = " + std::to_string(timeInteg->getTimeStep()));
		if (reader->getPolyOrder() >= 0) post->FRsolution("result" + std::to_string(iter), conserv[0], conserv[1], conserv[2]);
	}
	
	// Computed solution array
	std::vector<real_t> computed = orderTest->ZoneToPoly(conserv[0]);

	// Compute L errors
	real_t L1 = orderTest->L1error(computed);
	real_t L2 = orderTest->L2error(computed);
	real_t Linf = orderTest->Linf_error(computed);

	// Print L errors
	std::cout << "L1 error   = " << L1 << "\n";
	std::cout << "L2 error   = " << L2 << "\n";
	std::cout << "Linf error = " << Linf << "\n";

	// Post solution
	post->solution("result", conserv[0], conserv[1], conserv[2]);
	if (reader->getPolyOrder() > 0) post->FRsolution("result", conserv[0], conserv[1], conserv[2]);

	return 0;
}