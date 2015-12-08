//#include "chfunction_controller.h"
#include "ChFunction_Controller.h"
#include "Smarticle.h"
#include "Controller.h"
using namespace chrono;

double ChFunctionController::Get_y(double t) {
	double output = ComputeOutput(t);
	double output2 = std::max(std::min(controller_->outputLimit, output), -controller_->outputLimit);
	return output2;
}
double ChFunctionController::ComputeOutput(double t) {
	double friction = .2;
	

	static bool runOne = false;
	double p = controller_->p_gain;
	double i = controller_->i_gain;
	double d = controller_->d_gain;

	double curr_angle = controller_->GetCurrAngle(index_, t);
	double desired_angle = controller_->GetDesiredAngle(index_, t); ///get the next angle

	double error = desired_angle - curr_angle;

	cum_error_ += (error)*dT;
	double output = p*error + d*((error-prevError) / dT) + i*cum_error_;
	prevError = error;
	return output;
}
void ChFunctionController::ResetCumulative()	////TODO reset cum_error_ on gui change
{
	cum_error_ = 0;
}