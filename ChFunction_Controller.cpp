//#include "chfunction_controller.h"
#include "ChFunction_Controller.h"
#include "Smarticle.h"
#include "Controller.h"
using namespace chrono;

double ChFunctionController::Get_y(double t) {
	double output = ComputeOutput(t);
	double output2 = std::max(std::min(controller_->outputLimit, output), -controller_->outputLimit);
	//if (t>.5)
	//	GetLog() << "output: " << output << " output2: " << output2 << "\n";
	return output2;
}

// The low level PID controller in motor.
double ChFunctionController::ComputeOutput(double t) {
	double friction = .2;

	//double J = controller_->GetCurrTorque(index_, t) / friction;
	//double tau = J / friction;
	

	static bool runOne = false;
	double p = controller_->p_gain;
	double i = controller_->i_gain;
	double d = controller_->d_gain;

	double curr_angle = controller_->GetCurrAngle(index_, t);
	double desired_angle = controller_->GetDesiredAngle(index_, t); ///get the next angle

	double error = desired_angle - curr_angle;

	cum_error_ += (error)*dT;
	double output = p*error + d*((error-prevError) / dT) + i*cum_error_;
	//if (t > .3 && !runOne){
	//	runOne = true;
	//	GetLog() << "output: " << output << " output2: " << output << "\n";
	//	prevError = 0;
	//}
	prevError = error;
	return output;
}
void ChFunctionController::ResetCumulative()	////TODO reset cum_error_ on gui change
{
	cum_error_ = 0;
}