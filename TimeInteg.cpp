#include "TimeInteg.h"

TimeInteg::TimeInteg(Type PDEtype, Type fluxType, Type limiterType, real_t CFL, real_t targetTime, std::shared_ptr<Zone> zone, std::shared_ptr<Boundary> bdry)
{
	// Initializing variables
	_PDEtype = PDEtype; _fluxType = fluxType; _limiterType = limiterType;
	_CFL = CFL; _targetTime = targetTime;
	_currentTime = 0.0; _timeStep = 0.0;

	// Initializing objects
	_zone = zone; _bdry = bdry;
	_basis = std::make_shared<DGbasis>(zone->getPolyOrder(), zone->getGrid());

	// Initializing temporary variables
	_godFlux = std::make_shared<ConvFluxGodunov>(PDEtype, zone);
	_temp_solution.resize(zone->getGrid()->getNumCell());
	_prev_DOF.resize(zone->getPolyOrder() + 1);
	_temp_RHS.resize(zone->getPolyOrder() + 1);
	for (int_t iorder = 0; iorder <= zone->getPolyOrder(); ++iorder)
	{
		_prev_DOF[iorder].resize(zone->getGrid()->getNumCell());
		_temp_RHS[iorder].resize(zone->getGrid()->getNumCell());
	}
}

TimeInteg::~TimeInteg()
{

}

std::vector<std::vector<real_t> > TimeInteg::computeRHS(std::shared_ptr<Zone> zone) const
{
	real_t sizeX = zone->getGrid()->getSizeX();
	real_t inv_sizeX = 1.0/(zone->getGrid()->getSizeX());
	int_t num_cell = zone->getGrid()->getNumCell();
	int_t polyOrder = zone->getPolyOrder();

	// Temporary degree of freedom
	std::vector<std::vector<real_t> > DOF;
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();
	DOF.resize(polyOrder + 1);
	for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		DOF[iorder].resize(num_cell);
	
	// DG flux
	std::vector<real_t> flux(num_cell, 0.0);
	real_t for_projec, back_projec, left_u, right_u;
	for (int_t icell = GHOST; icell <= num_cell - GHOST; ++icell)
	{
		// Calculate flux inputs
		// Local projection
		if (polyOrder > 0)
		{
			for_projec = PROJEC_COEFF2*temp_DOF[1][icell - 1];
			back_projec = PROJEC_COEFF2*temp_DOF[1][icell];
		}
		if (polyOrder > 1)
		{
			for_projec += PROJEC_COEFF3*temp_DOF[2][icell - 1];
			back_projec += -PROJEC_COEFF3*temp_DOF[2][icell];
		}

		// Cell quantity
		left_u = PROJEC_COEFF1 * temp_DOF[0][icell - 1] + for_projec;
		right_u = PROJEC_COEFF1 * temp_DOF[0][icell] - back_projec;

		// Calculate flux
		flux[icell] = _godFlux->computeFlux(left_u, right_u);
	}

	// Calculate DOF
	for (int_t icell = GHOST; icell < num_cell - GHOST; ++icell)
	{
		// degree 0
		DOF[0][icell] = -inv_sizeX*(flux[icell + 1] - flux[icell]);

		// degree 1
		if (polyOrder > 0)
		{
			DOF[1][icell] = -0.5*inv_sizeX*(flux[icell + 1] + flux[icell]);
			real_t temp_x = zone->getGrid()->getCell()[icell]->getPosX();
			for (int_t idegree = 0; idegree < QuadDegree; ++idegree)
				DOF[1][icell] += 0.5*inv_sizeX*Gauss3_W(idegree)*PHY_FLUX(_PDEtype, zone->getPolySolution(icell, temp_x + 0.5*sizeX*Gauss3_X(idegree)));
		}

		// degree 2
		if (polyOrder > 1)
		{
			DOF[2][icell] = -CONST16*inv_sizeX*(flux[icell + 1] - flux[icell]);
			real_t temp_x = zone->getGrid()->getCell()[icell]->getPosX();
			for (int_t idegree = 0; idegree < QuadDegree; ++idegree)
				DOF[2][icell] += pow(inv_sizeX, 2.0)*Gauss3_W(idegree)*PHY_FLUX(_PDEtype, zone->getPolySolution(icell, temp_x + 0.5*sizeX*Gauss3_X(idegree)))
				*_basis->basis(1, icell, temp_x + 0.5*sizeX*Gauss3_X(idegree));
		}
	}
	
	return DOF;
}

void TimeInteg::computeTimeStep(std::shared_ptr<Zone> zone)
{
	if (_PDEtype == "advection")
		_timeStep = _CFL*zone->getGrid()->getSizeX() / abs(GET_SPEED) / double(2 * zone->getPolyOrder() + 1);

	else if (_PDEtype == "burgers")
	{
		real_t temp_sol1;
		real_t temp_sol2;
		std::vector<real_t> shockSpeed(zone->getGrid()->getNumCell() - 1, 0.0);
		for (int_t icell = 0; icell < shockSpeed.size(); ++icell)
		{
			temp_sol1 = zone->getDescSolution()[icell];
			temp_sol2 = zone->getDescSolution()[icell + 1];
			if (temp_sol1 >= temp_sol2) shockSpeed[icell] = 0.5*abs(temp_sol1 + temp_sol2);
			else shockSpeed[icell] = std::max(abs(temp_sol1), abs(temp_sol2));
		}
		_timeStep = _CFL*zone->getGrid()->getSizeX() / *std::max_element(shockSpeed.begin(), shockSpeed.end()) / double(2*zone->getPolyOrder() + 1);
	}
}

void TimeInteg::print() const
{
	MESSAGE("Marching finished.....");
	MESSAGE("CFL = " + std::to_string(_CFL));
	MESSAGE("Target time = " + std::to_string(_targetTime));
	MESSAGE("Final time step = " + std::to_string(_timeStep));
	MESSAGE("Current time = " + std::to_string(_currentTime));
}