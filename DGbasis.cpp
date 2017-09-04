#include "DGbasis.h"

DGbasis::DGbasis(int_t order, std::shared_ptr<Grid> grid)
{
	_grid = grid;
	_cell = grid->getCell();
	_sizeX = grid->getSizeX();
	_coeff.resize(order + 1);
	_coeff[0] = 1.0;
	if (order > 0) _coeff[1] = 12.0 / _sizeX;
	if (order > 1) _coeff[2] = 180.0 / pow(_sizeX, 2.0);
	if (order > 2) ERROR("Exceed maximum polynomial order");
}

DGbasis::~DGbasis()
{

}

real_t DGbasis::basis(int_t degree, int_t index, real_t x)
{
	switch (degree)
	{
	case 0: return 1.0;
	case 1: return (x - _cell[index]->getPosX());
	case 2: return (pow(x - _cell[index]->getPosX(), 2.0) - pow(_sizeX, 2.0) / 12.0);
	default: return 0.0;
	}
}