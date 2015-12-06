//#include "chfunction_controller.h"
#include "ChFunction_Controller.h"
#include "Smarticle.h"
#include "Controller.h"
using namespace chrono;
double ChFunctionController::Get_y(double t) {
	double output = ComputeOutput(t);
	
	//double output = controller_->computeOutput();
	output = std::max(std::min(controller_->outputLimit, output), -controller_->outputLimit);
	//controller_->~Controller();
	return output;
}

// The low level PID controller in motor.
double ChFunctionController::ComputeOutput(double t) {
	double friction = .2;

	//double J = controller_->GetCurrTorque(index_, t) / friction;
	//double tau = J / friction;

	double p = controller_->p_gain;
	double i = controller_->i_gain;
	double d = controller_->d_gain;

	//double curr_angle = controller_->GetCurrAngle(index_, t);
	//double desired_angle = controller_->GetDesiredAngle(index_, t); ///get the next angle
	//double curr_angular_speed = controller_->GetCurrAngularSpeed(index_,t);
	//double output = controller_->GetDesiredAngularSpeed2(index_, t,(desired_angle-curr_angle));


	double curr_angle = controller_->GetCurrAngle(index_, t);
	double desired_angle = controller_->GetDesiredAngle(index_, t); ///get the next angle
	double curr_angular_speed = controller_->GetCurrAngularSpeed(index_, t);
	double desired_angular_speed = controller_->GetDesiredAngularSpeedForFunction(index_, t);
	cum_error_ += (desired_angle - curr_angle)*dT;
	double output = p * (desired_angle - curr_angle) +
		d*  (desired_angular_speed - curr_angular_speed) +
		cum_error_ * i;

	return output;
}
