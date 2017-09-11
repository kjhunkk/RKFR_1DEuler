#include "InitialCondition.h"

InitialCondition::InitialCondition(Type type)
{
	_type = type;
}

InitialCondition::~InitialCondition()
{

}

real_t InitialCondition::initializer(real_t x, int_t solType) const
{
	if (_type == "smooth") // smooth for order test
		return smooth(x, solType);

	if (_type == "shocktube") // Sod's shock tube
		return shocktube(x, solType);

	if (_type == "supersonicexpansion") // supersonic expansion
		return expansion(x, solType);

	if (_type == "blastwave") // Blast wave
		return blast(x, solType);

	if (_type == "twostrongshock") // two strong shock waves
		return double_shock(x, solType);

	if (_type == "contactdiscontinuity") // slowly-moving contact discontinuities
		return contact(x, solType);

	if (_type == "blastwaveinteraction") // blast wave interaction
		return blastInteraction(x, solType);

	if (_type == "shuosher") // Shu-Osher
		return Shu_Osher(x, solType);

	ERROR("cannot find initial condition");

	return 0.0;
}

real_t InitialCondition::smooth(real_t x, int_t solType) const
{
	switch (solType)
	{
	case 0: return 1.0 + 0.2*sin(2.0*M_PI*x);
	case 1: return 1.0 + 0.2*sin(2.0*M_PI*x);
	case 2: return 0.5*(1.0 + 0.2*sin(2.0*M_PI*x)) + 1.0 / (GAMMA - 1);
	default:
		ERROR("solution type number");
		return 0.0;
	}
}

real_t InitialCondition::shocktube(real_t x, int_t solType) const
{
	switch (solType)
	{
	case 0:
		if (x <= 0.3) return 1.0;
		else return 0.125;
	case 1:
		if (x <= 0.3) return 0.75;
		else return 0.0;
	case 2:
		if (x <= 0.3) return 2.78125;
		else return 0.25;
	default:
		ERROR("solution type number");
		return 0.0;
	}
}

real_t InitialCondition::expansion(real_t x, int_t solType) const
{
	switch (solType)
	{
	case 0: return 1.0;
	case 1:
		if (x <= 0.5) return -2.0;
		else return 2.0;
	case 2:
		if (x <= 0.5) return 3.0;
		else return 3.0;
	default:
		ERROR("solution type number");
		return 0.0;
	}
}

real_t InitialCondition::blast(real_t x, int_t solType) const
{
	switch (solType)
	{
	case 0: return 1.0;
	case 1: return 0.0;
	case 2:
		if (x <= 0.5) return 2500.0;
		else return 0.025;
	default:
		ERROR("solution type number");
		return 0.0;
	}
}

real_t InitialCondition::double_shock(real_t x, int_t solType) const
{
	switch (solType)
	{
	case 0:
		if (x <= 0.4) return 5.99924;
		else return 5.99242;
	case 1:
		if (x <= 0.4) return 117.5701059;
		else return -37.13101182;
	case 2:
		if (x <= 0.4) return 2304.275075;
		else return 230.2755012;
	default:
		ERROR("solution type number");
		return 0.0;
	}
}

real_t InitialCondition::contact(real_t x, int_t solType) const
{
	switch (solType)
	{
	case 0: return 1.0;
	case 1: return -19.59745;
	case 2: 
		if (x <= 0.8) return 2692.030023;
		else return 192.0550233;
	default:
		ERROR("solution type number");
		return 0.0;
	}
}

real_t InitialCondition::blastInteraction(real_t x, int_t solType) const
{
	switch (solType)
	{
	case 0: return 1.0;
	case 1: return 0.0;
	case 2:
		if (x <= -4.0) return 2500.0;
		else if (x >= 4.0) return 0.025;
		else return 250.0;
	default:
		ERROR("solution type number");
		return 0.0;
	}
}

real_t InitialCondition::Shu_Osher(real_t x, int_t solType) const
{
	switch (solType)
	{
	case 0:
		if (x < -4.0) return 3.857143;
		else return 1.0 + 0.2*sin(5.0*x);
	case 1:
		if (x < -4.0) return 10.14185223;
		else return 0.0;
	case 2:
		if (x < -4.0) return 39.16666843;
		else return 2.5;
	}
}