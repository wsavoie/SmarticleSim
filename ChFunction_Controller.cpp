//#include "chfunction_controller.h"
#include "ChFunction_Controller.h"
#include "Smarticle.h"
#include "Controller.h"
using namespace chrono;

double ChFunctionController::Get_y(double t) {
	double curr_react_torque = controller_->GetCurrTorque(index_, t);
  //add the torque already being place on the body to the torque for the next step
	//double out_torque =output; //add the torque already being place on the body to the torque for the next step
	double output = ComputeOutput(t);
	double out_torque = SaturateValue(ComputeOutput(t) - curr_react_torque, controller_->outputLimit);
		
	bool o = false;
	if (abs(out_torque) == controller_->outputLimit)
		o = true;
	return out_torque;
}
double ChFunctionController::ComputeOutput(double t) {
	
	double curr_ang = controller_->GetAngle(index_,t);
	double exp_ang = controller_->GetExpAngle(index_, t);
	double des_ang = controller_->GetDesiredAngle(index_, t); ///get the next angle
	des_ang = controller_->LinearInterpolate(index_, curr_ang, des_ang); //linear interpolate for situations where gui changes so there isn't a major speed increase
	double error = des_ang - curr_ang;
	double prevError = controller_->prevError_.at(index_);

	double K = p_gain;
	double Ti = i_gain;
	double Td = d_gain;
	double Tt = 1*dT;//read about tt
	double N = 10; //N=[8-20] http://www.cds.caltech.edu/~murray/courses/cds101/fa02/caltech/astrom-ch6.pdf
	double b = 1;

	double ulimLow = this->controller_->smarticle_->angLow;
	double ulimHigh = this->controller_->smarticle_->angHigh;
	double vlim = controller_->omegaLimit;
	double tlim = controller_->outputLimit;

	double bi = K*dT / Ti;//integral gain
	double ad = (2 * Td - N*dT) / (2 * Td + N*dT);
	double bd = 2 * K*N*Td / (2 * Td + N*dT); //deriv gain
	double ao = dT / Tt;
	double ysp = des_ang;
	double y = curr_ang;

	//initializes yold to current value for first iteration
	if (controller_->smarticle_->steps == 0)
		controller_->yold[index_] = y;

	double pp = K*(b*ysp - y);
	controller_->DD[index_] = ad*controller_->DD[index_] - bd*(y - controller_->yold[index_]);
	double v = pp + controller_->II[index_] + controller_->DD[index_];
	double u = SaturateValue(v, tlim);
	controller_->II[index_] = controller_->II[index_] + bi*(ysp - y) + ao*(u - v);
	controller_->yold[index_] = y;
	

	//double vel = (controller_->yold[index_] - y) / dT;
	//vel = SaturateValue(vel, vlim);
	//double tor = (controller_->velOld[index_] - vel) / dT;
	//tor = SaturateValue(tor, tlim);
	//controller_->velOld[index_] = vel;



	double out = u;


	if (controller_->resetCumError)
		ResetCumulative(t);

	return out;
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
	controller_->cumOmegError_.at(index_) = 0;
	bothArmsReset++; //add one to value, if value>1, both arms have been reset thus value can be set to false
	controller_->II[index_] = 0;
	if (bothArmsReset >1)
	{
		//GetLog() << "reset\n";
		controller_->resetCumError = false;
		bothArmsReset = 0;
	}
}