#include "Limiter.h"

Limiter::Limiter(Type limiter, int_t polyOrder)
{
	_limiter = limiter;
	_polyOrder = polyOrder;
	_h = 0.0;
	_M = CONST23*10.0;
}

Limiter::Limiter(Type limiter, int_t polyOrder, real_t h, real_t M)
{
	_limiter = limiter;
	_polyOrder = polyOrder;
	_h = h;
	_M = M;
}

Limiter::~Limiter()
{

}

void Limiter::projectLimiter(std::shared_ptr<Zone> zone)
{
	// No limiter if PO
	if ((_polyOrder == 0) || (_limiter == "none")) return;

	// Get degree of freedom from Zone object
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();

	// Temporary zone object
	std::shared_ptr<Zone> temp_zone = zone;
	int_t num_cell = zone->getGrid()->getNumCell();

	// Local projection
	std::vector<real_t> for_projec, back_projec;
	if (_polyOrder > 0)
	{
		for (int_t icell = 0; icell < num_cell; ++icell)
		{
			for_projec.push_back(PROJEC_COEFF2*temp_DOF[1][icell]);
			back_projec.push_back(PROJEC_COEFF2*temp_DOF[1][icell]);
		}

	}
	if (_polyOrder > 1)
	{
		for (int_t icell = 0; icell < num_cell; ++icell)
		{
			for_projec[icell] += PROJEC_COEFF3*temp_DOF[2][icell];
			back_projec[icell] -= PROJEC_COEFF3*temp_DOF[2][icell];
		}
	}

	// Modified local projection
	std::vector<real_t> mod_for_projec, mod_back_projec;
	mod_for_projec = for_projec;
	mod_back_projec = back_projec;
	real_t for_diff;
	real_t back_diff;

	if (_limiter == "minmod")
	{
		for (int_t icell = GHOST; icell < num_cell - GHOST; ++icell)
		{
			for_diff = temp_DOF[0][icell + 1] - temp_DOF[0][icell];
			back_diff = temp_DOF[0][icell] - temp_DOF[0][icell - 1];
			mod_for_projec[icell] = minmod(for_projec[icell], for_diff, back_diff);
			mod_back_projec[icell] = minmod(back_projec[icell], for_diff, back_diff);
		}
	}
	else if (_limiter == "modified_minmod")
	{
		for (int_t icell = GHOST; icell < num_cell - GHOST; ++icell)
		{
			for_diff = temp_DOF[0][icell + 1] - temp_DOF[0][icell];
			back_diff = temp_DOF[0][icell] - temp_DOF[0][icell - 1];
			mod_for_projec[icell] = mod_minmod(for_projec[icell], for_diff, back_diff);
			mod_back_projec[icell] = mod_minmod(back_projec[icell], for_diff, back_diff);
		}
	}
	else ERROR("cannot find limiter");

	// Update limited DOF
	if (_polyOrder == 1)
	{
		for (int_t icell = 0; icell < num_cell; ++icell)
			temp_DOF[1][icell] = CONST16*mod_for_projec[icell];
	}
	else if (_polyOrder == 2)
	{
		for (int_t icell = 0; icell < num_cell; ++icell)
		{
			temp_DOF[1][icell] = CONST1_12*(mod_for_projec[icell] + mod_back_projec[icell]);
			temp_DOF[2][icell] = CONST1_60*(mod_for_projec[icell] - mod_back_projec[icell]);
		}
	}
	else ERROR("exceed maximum polynomial order");

	// Set limited DOF and return
	zone->setDOF(temp_DOF);
}

real_t Limiter::minmod(real_t projection, real_t for_diff, real_t back_diff) const
{
	std::vector<real_t> temp = { abs(projection), abs(for_diff), abs(back_diff) };
	real_t min = *min_element(temp.begin(), temp.end());
	
	if (projection > epsilon)
	{
		if ((for_diff > epsilon) && (back_diff > epsilon)) return min;
		else return 0.0;
	}
	else if (projection < -epsilon)
	{
		if ((back_diff < -epsilon) && (back_diff < -epsilon)) return -min;
		else return 0.0;
	}
	else return 0.0;
}

real_t Limiter::mod_minmod(real_t projection, real_t for_diff, real_t back_diff) const
{
	if (_h < epsilon) ERROR("cannot use modified minmod : no grid size");
	if (abs(projection) <= _M*pow(_h, 2.0)) return projection;
	else return minmod(projection, for_diff, back_diff);
}