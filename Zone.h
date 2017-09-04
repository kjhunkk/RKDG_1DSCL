#pragma once
#include "DataType.h"
#include "Grid.h"
#include "InitialCondition.h"
#include "Quadrature.h"
#include "DGbasis.h"

class Zone
{
public:
	// Constructor / p.m. Grid(object)
	Zone(std::shared_ptr<Grid>);

	// Constructor / p.m. Grid(object), DG polynomial order
	Zone(std::shared_ptr<Grid>, int_t);

	// Destructor
	~Zone();

public:
	// Functions
	inline std::shared_ptr<Grid> getGrid() const { return _grid; }

	// Get Descrete solution
	inline std::vector<real_t> getDescSolution() const { return _solution; }

	inline std::vector<std::vector<real_t> > getDOF() const { return _DOF; }

	inline int_t getPolyOrder() const { return _polyOrder; }

	// Set Descrete solution
	inline void setDescSolution(std::vector<real_t> solution) { _solution = solution; }

	inline void setDOF(std::vector<std::vector<real_t> > DOF) { _DOF = DOF; }

	// Get polynomial solution at coordinate x / p.m. cell index, x coordinate
	real_t getPolySolution(int_t, real_t) const;

	// Get polynomial solution at coordinate x from specific Zone / p.m. cell index, x coordinate, Zone(object)
	real_t getPolySolution(int_t, real_t, std::shared_ptr<Zone>) const;

	// Initialize solution / p.m. Initial condition(object)
	void initialize(std::shared_ptr<InitialCondition>);

	// Calculate Descrete solution from DOF
	void calSolution();

	// Print Solution variables
	void print() const;

protected:
	// Variables
	std::shared_ptr<Grid> _grid;
	std::shared_ptr<DGbasis> _basis;
	std::vector<real_t> _solution;
	std::vector<std::vector<real_t> > _DOF;
	int_t _polyOrder;
};