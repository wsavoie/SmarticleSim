//#include "chfunction_controller.h"
#include "ChFunction_Controller.h"
#include "Smarticle.h"
#include "Controller.h"
using namespace chrono;

//double ChFunctionController::Get_y(double t) {
//	double output = ComputeOutput(t);
//	double a = controller_->GetCurrTorque(index_,t);
//	//if (a > .00001)
//	//	GetLog() << "hightorque on arm: " << index_ << "\n";
//	output = std::max(std::min(controller_->outputLimit, output), -controller_->outputLimit);
//	return output;
//}
//double ChFunctionController::ComputeOutput(double t) {
//	double friction = .2;
//	
//	double p = controller_->p_gain;
//	double i = controller_->i_gain;
//	double d = controller_->d_gain;
//
//	double prev_angle = controller_->GetCurrAngle(index_, t);
//	controller_->prevAngle_.at(index_) = prev_angle;
//	
//	double curr_angle = controller_->SetGetCurrAngle(index_, t);
//	double desired_angle = controller_->GetDesiredAngle(index_, t); ///get the next angle
//	desired_angle = controller_->LinearInterpolate(index_,curr_angle,desired_angle);
//	double error = desired_angle - curr_angle;
//	double prevError = controller_->prevError_.at(index_);
//
//
//
//	if (controller_->resetCumError)
//		ResetCumulative();
//	cum_error_ += (error)*dT;
//	double output = p*error + d*((error - prevError) / dT) + i*cum_error_;
//	controller_->prevError_.at(index_) = error;
//	return output;
//}

double ChFunctionController::Get_y(double t) {
	double output = ComputeOutput(t);
	//double curr_react_torque = controller_->GetCurrTorque(index_, t);
; //add the torque already being place on the body to the torque for the next step
	//double out_torque =output; //add the torque already being place on the body to the torque for the next step
	GetLog() << "Torque: " << output;
	double out_torque = std::max(std::min(controller_->outputLimit, output), -controller_->outputLimit);
	bool o = false;
	if (abs(out_torque) == controller_->outputLimit)
		o = true;
	GetLog() << "   \tLimited: " << out_torque << " " << "\tLIMIT= " << controller_->outputLimit << "\t" << o<< "\n";

	return out_torque;
	//return output;
}



double ChFunctionController::Get_y2(double t) {
	double output = ComputeOutput(t);
	double out_omega = OutputToOmega(t, output);
	out_omega=controller_->OmegaLimiter(index_, out_omega);
	double out_torque = OmegaToTorque(t, out_omega);

	double curr_react_torque = controller_->GetCurrTorque(index_, t);
	out_torque = curr_react_torque + out_torque; //add the torque already being place on the body to the torque for the next step

	out_torque = std::max(std::min(controller_->outputLimit, out_torque), -controller_->outputLimit);
	return out_torque;
	//return output;
}


double ChFunctionController::ComputeOutput(double t) {
	double p =p_gain;
	double i =i_gain;
	double d =d_gain;

	double curr_angle = controller_->GetAngle(index_, t);
	//GetLog() << "Controller" << index_ << ": " <<curr_angle<<"\n";
	double exp_angle = controller_->GetExpAngle(index_, t);
	double desired_angle = controller_->GetDesiredAngle(index_, t); ///get the next angle
	
	//double desired_angle2 = controller_->LinearInterpolate(index_, curr_angle, desired_angle);
	desired_angle = controller_->LinearInterpolate(index_, curr_angle, desired_angle);
	double error = desired_angle - curr_angle;
	//double error2 = desired_angle2- curr_angle;
	double prevError = controller_->prevError_.at(index_);
	
	//if position was changed resetCumError flag is set to true
	if (controller_->resetCumError)
		ResetCumulative(t);

	controller_->cumError_.at(index_) += (error)*dT;

	//calculate controller output in position
	//double deri = d*((error-prevError) / dT); ///fix negative sign vs pos sign error
	//double integ = i*controller_->cumError_.at(index_);
	double output = p*error + d*((error - prevError) / dT) + i*controller_->cumError_.at(index_);

	//set current error to previous error
	controller_->prevError_.at(index_) = error;
	return output;
}
double ChFunctionController::OutputToOmega(double t, double out) {
	//double curr_ang = controller_->GetCurrAngle(index_,t);
	double curr_ang = controller_->GetAngle(index_, t);
	double d_ang = out - curr_ang;
	double omegaOut = d_ang / dT;
	return omegaOut;
}
double ChFunctionController::OmegaToTorque(double t, double out)
{
	double curr_omega = controller_->GetActuatorOmega(index_, t);
	controller_->GetActuatorOmega(index_, t);
	double d_omega = out - curr_omega;
	return d_omega / dT;
}
void ChFunctionController::ResetCumulative(double t)	////TODO reset cum_error_ on gui change
{
	
	static int bothArmsReset = 0;
	controller_->cumError_.at(index_) = 0;
	//GetLog() << "dt="<<t<<"\n";
	bothArmsReset++; //add one to value, if value>1, both arms have been reset thus value can be set to false
	if (bothArmsReset >1)
	{
		//GetLog() << "reset\n";
		controller_->resetCumError = false;
		bothArmsReset = 0;
	}
}