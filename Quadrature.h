#pragma once
#include "DataType.h"

class Quadrature
{
protected:
	Quadrature();
	~Quadrature();

public:
	// Gauss-Legendre guadrature x and weight / p.m. index
	static real_t gauss2X(int_t);
	static real_t gauss3X(int_t);

	static real_t gauss2W(int_t);
	static real_t gauss3W(int_t);

private:
	// Quadrature variable
	static Quadrature _quad;
};

// Gauss-Legendre quadrature macro
// index : 0~
#define Gauss2_X(index) Quadrature::gauss2X(index)
#define Gauss2_W(index) Quadrature::gauss2W(index)

#define Gauss3_X(index) Quadrature::gauss3X(index)
#define Gauss3_W(index) Quadrature::gauss3W(index)