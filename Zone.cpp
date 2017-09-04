#include "Zone.h"

Zone::Zone(std::shared_ptr<Grid> grid)
{
	_grid = grid;
	_polyOrder = 0;
	_solution.resize(grid->getNumCell());
	_DOF.resize(1);
	_DOF[0].resize(grid->getNumCell());

	// DG basis
	_basis = std::make_shared<DGbasis>(_polyOrder, _grid);
}

Zone::Zone(std::shared_ptr<Grid> grid, int_t polyOrder)
{
	_grid = grid;
	_polyOrder = polyOrder;

	// Print error message if polynomial order is less than 0
	if (_polyOrder < 0)
		ERROR("Polynomial order");

	_solution.resize(grid->getNumCell());
	_DOF.resize(polyOrder + 1);
	for (int_t iorder = 0; iorder <= polyOrder; ++iorder)
		_DOF[iorder].resize(grid->getNumCell());

	// DG basis
	_basis = std::make_shared<DGbasis>(_polyOrder, _grid);
}

Zone::~Zone()
{

}

real_t Zone::getPolySolution(int_t icell, real_t x) const
{
	real_t u = 0;
	
	// Calculate polynomial solution at x / beware of time level
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		u += _basis->getCoeff()[idegree] * _DOF[idegree][icell] * _basis->basis(idegree, icell, x);

	return u;
}

real_t Zone::getPolySolution(int_t icell, real_t x, std::shared_ptr<Zone> zone) const
{
	real_t u = 0;
	// Calculate polynomial solution at x / beware of time level
	for (int_t idegree = 0; idegree <= _polyOrder; ++idegree)
		u += _basis->getCoeff()[idegree] * zone->getDOF()[idegree][icell] * _basis->basis(idegree, icell, x);

	return u;
}

void Zone::initialize(std::shared_ptr<InitialCondition> initialCondition)
{
	// Temporary cell object
	std::vector<std::shared_ptr<Cell> > temp_cell = _grid->getCell();

	// Grid size
	real_t temp_dx = _grid->getSizeX();

	// Initializing Degree of freedom
	for (int_t icell = 0; icell < _grid->getNumCell(); ++icell)
	{
		real_t temp_x = temp_cell[icell]->getPosX();
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
		{
			_DOF[iorder][icell] = 0.0;
			for (int_t idegree = 0; idegree < QuadDegree; ++idegree)
			{
				_DOF[iorder][icell] += Gauss3_W(idegree)*initialCondition->initializer(temp_x + 0.5*temp_dx*Gauss3_X(idegree))
					*_basis->basis(iorder, icell, temp_x + 0.5*temp_dx*Gauss3_X(idegree))*0.5;
			}
			_DOF[iorder][icell] /= pow(temp_dx, iorder);
		}
	}

	// Calculate solution
	calSolution();
}

void Zone::calSolution()
{
	for (int_t icell = 0; icell < _grid->getNumCell(); ++icell)
	{
		_solution[icell] = 0.0;
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
		{
			_solution[icell] += _basis->getCoeff()[iorder] * _DOF[iorder][icell] * _basis->basis(iorder, icell, _grid->getCell()[icell]->getPosX());
		}
	}
}

void Zone::print() const
{
	std::cout << "X coordinate" << "\t\t" << "Solution" << "\t\t";
	for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
		std::cout << "DOF(" << iorder << ")\t\t\t";
	std::cout << "\n";

	// Print variables
	for (int_t icell = 0; icell < _grid->getNumCell(); ++icell)
	{
		std::cout << std::to_string(_grid->getCell()[icell]->getPosX()) << "\t\t" << std::to_string(_solution[icell]) << "\t\t";
		for (int_t iorder = 0; iorder <= _polyOrder; ++iorder)
			std::cout << std::to_string(_DOF[iorder][icell]) << "\t\t";
		std::cout << "\n";
	}
}