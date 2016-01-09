//#include "chfunction_controller.h"
#include "ChFunction_Controller.h"
#include "Smarticle.h"
#include "Controller.h"
using namespace chrono;


double ChFunctionController::Get_y(double t) {
	double output = ComputeOutput(t);
	//double curr_react_torque = controller_->GetCurrTorque(index_, t);
  //add the torque already being place on the body to the torque for the next step
	//double out_torque =output; //add the torque already being place on the body to the torque for the next step
	//GetLog() << "Torque: " << output;
	double out_torque = std::max(std::min(controller_->outputLimit, output), -controller_->outputLimit);
	bool o = false;
	if (abs(out_torque) == controller_->outputLimit)
		o = true;
	//GetLog() << "   \tLimited: " << out_torque << " " << "\tLIMIT= " << controller_->outputLimit; //<< "\t" << o<< "\n";
	
	////////////
	//out_torque = SpeedControl(t, out_torque);
	//double sp = controller_->GetActuatorOmega(index_,t);
	//GetLog() << "\tSpeed" << controller_->GetActuatorOmega(index_, t) << "\n";
	/////////
	return out_torque;
	//return output;
}

double ChFunctionController::ComputeOutput(double t) {
	double p = p_gain;
	double i = i_gain;
	double d = d_gain;

	double curr_ang = controller_->GetAngle(index_, t);
	double exp_ang = controller_->GetExpAngle(index_, t);
	double des_ang = controller_->GetDesiredAngle(index_, t); ///get the next angle

	des_ang = controller_->LinearInterpolate(index_, curr_ang, des_ang); //linear interpolate for situations where gui changes so there isn't a major speed increase

	double error = des_ang - curr_ang;
	double prevError = controller_->prevError_.at(index_);

	//if position was changed resetCumError flag is set to true
	if (controller_->resetCumError)
		ResetCumulative(t);

	controller_->cumError_.at(index_) += (error)*dT;

	double pt = p*error;
	double it = i*controller_->cumError_.at(index_);
	double dt = d*((error - prevError) / dT);
	double output = p*error + d*((error - prevError) / dT) + i*controller_->cumError_.at(index_);

	//set current error to previous error
	controller_->prevError_.at(index_) = error;



	//////////////////
	//Speed control
	double omLim = controller_->omegaLimit;
	double omErrCond = omLim*dT;
	int cond = (abs(error) > omErrCond);

	double curr_omeg = controller_->GetActuatorOmega(index_, t);

	double des_omeg = omLim*cond*sgn(error);
	
	double omError = des_omeg - curr_omeg;
	double prevOmError = controller_->prevOmegError_.at(index_);
	controller_->cumOmegError_.at(index_) += (omError)*dT;
	controller_->prevOmegError_.at(index_) = omError;
	
	double pTerm = p*dT*omLim*omError;
	double iTerm = i*dT*omLim*controller_->cumOmegError_.at(index_);
	double dTerm = d*dT*omLim*(omError - prevOmError / dT);
	double dTerm2 = d*dT*omLim*((omError - prevOmError) / dT);
	//double output2 = 100 * p*dT*omLim*omError;  ///TODO magic number 100?
	double output2 = 100 * p*dT*omLim*omError;  ///TODO magic number 100?
										//+ d*dT*omLim*((omError - prevOmError) / dT);
										//+ i*dT*omLim*controller_->cumOmegError_.at(index_);
	//GetLog() << "output" << output << "\toutput2" << output2;

	
	
	double output3 = output+output2;
	//GetLog() << "1: " << output << " " << "\t2: " << output2 << "\t3: " << output3 << "\tom: " << curr_omeg<<"\n";
	//////////////////







	return output3;
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
	//GetLog() << "dt="<<t<<"\n";
	bothArmsReset++; //add one to value, if value>1, both arms have been reset thus value can be set to false
	if (bothArmsReset >1)
	{
		//GetLog() << "reset\n";
		controller_->resetCumError = false;
		bothArmsReset = 0;
	}
}