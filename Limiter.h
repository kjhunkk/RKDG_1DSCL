#pragma once
#include "DataType.h"
#include "Zone.h"

class Limiter
{
public:
	// Constructor / p.m. limiter type, polynomial order
	Limiter(Type, int_t);

	// Constructor / p.m. limiter type, polynomial order, grid size, M constant(M = 2/3*M2, M2 = max(abs(u_xx^0)))
	Limiter(Type, int_t, real_t, real_t);

	// Desturctor
	~Limiter();

public:
	// Functions
	// calculate local projection limiter / p.m. Zone(object)
	void projectLimiter(std::shared_ptr<Zone>);

protected:
	// Variables
	Type _limiter;
	int_t _polyOrder;
	real_t _M;
	real_t _h;

protected:
	// Functions
	// Minmod function / p.m. local projection, forward difference, backward difference
	real_t minmod(real_t, real_t, real_t) const;

	// Modified minmod function / p.m. local projection, forward difference, backward difference
	real_t mod_minmod(real_t, real_t, real_t) const;

};