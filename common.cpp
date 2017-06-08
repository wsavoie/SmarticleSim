#include "common.h"
using namespace chrono;




double SaturateValue(double val, double low, double high)
{
	return std::max(std::min(high, val), low);
}
//for single input centered around 0
double SaturateValue(double val, double zeroCenteredVal)
{
	//in case inputted negative value
	return SaturateValue(val, -abs(zeroCenteredVal), abs(zeroCenteredVal));
}

//generates random [min,max]
double genRand(double min, double max)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(min, max);
	return dis(gen);
}
//generates random [0,max]
double genRand(double max)
{
	return genRand(0, max);
}
//generates random [0,1]
double genRand()
{
	return genRand(0, 1);
}
int genRandInt(int min, int max)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> intdis(min, max);
	return intdis(gen);
}