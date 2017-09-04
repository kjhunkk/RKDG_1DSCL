#include "ConvFlux.h"

ConvFlux::ConvFlux(Type phyFlux, std::shared_ptr<Zone> zone)
{
	_zone = zone;
	_phyFlux = phyFlux;
}

ConvFlux::~ConvFlux()
{

}