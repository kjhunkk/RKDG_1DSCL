#pragma once
#include "DataType.h"
#include "Zone.h"

class ConvFlux
{
public:
	// Constructor / p.m. physical flux type, Zone(Object)
	ConvFlux(Type, std::shared_ptr<Zone>);

	// Destructor
	virtual ~ConvFlux();

public:
	// Functions
	// Set Zone
	void setZone(std::shared_ptr<Zone> zone) { _zone = zone; }

	// Compute flux / p.m. begin, end
	virtual real_t computeFlux(real_t, real_t) const = 0;

protected:
	// Variables
	std::shared_ptr<Zone> _zone;
	Type _phyFlux;
};