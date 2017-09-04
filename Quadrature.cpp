#include "Quadrature.h"

Quadrature::Quadrature()
{

}

Quadrature::~Quadrature()
{

}

real_t Quadrature::gauss2X(int_t index)
{
	switch (index)
	{
	case 0: return -CONST13*sqrt(3.0);
	case 1: return CONST13*sqrt(3.0);
	default: return 0.0;
	}
}

real_t Quadrature::gauss3X(int_t index)
{
	switch (index)
	{
	case 0: return -0.2*sqrt(15.0);
	case 1: return 0.0;
	case 2: return 0.2*sqrt(15.0);
	default: return 0.0;
	}
}

real_t Quadrature::gauss2W(int_t index)
{
	return 1.0;
}

real_t Quadrature::gauss3W(int_t index)
{
	switch (index)
	{
	case 0: return CONST59;
	case 1: return CONST89;
	case 2: return CONST59;
	default: return 0.0;
	}
}

Quadrature Quadrature::_quad;