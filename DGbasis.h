#pragma once
#include "DataType.h"
#include "Grid.h"

class DGbasis
{
public:
	// Constructor / p.m. polynomial order, Gird(object)
	DGbasis(int_t, std::shared_ptr<Grid>);

	// Destructor
	~DGbasis();

public:
	// Functions
	inline vector_d getCoeff() const { return _coeff; }

	// Basis function / p.m. degree, cell index, x coordinate
	real_t basis(int_t, int_t, real_t);

protected:
	// Variables
	std::shared_ptr<Grid> _grid;
	std::vector<std::shared_ptr<Cell> > _cell;
	std::vector<real_t> _coeff;
	real_t _sizeX;
};