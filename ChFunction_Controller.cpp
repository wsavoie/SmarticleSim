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
	double out_torque = SaturateValue(output, controller_->outputLimit);
	bool o = false;
	if (abs(out_torque) == controller_->outputLimit)
		o = true;
	//GetLog() << "   \tLimited: " << out_torque << " " << "\tLIMIT= " << controller_->outputLimit; //<< "\t" << o<< "\n";
	//GetLog() << "\n";
	////////////
	//out_torque = SpeedControl(t, out_torque);
	//double sp = controller_->GetActuatorOmega(index_,t);
	//GetLog() << "\tSpeed" << controller_->GetActuatorOmega(index_, t) << "\n";
	/////////
	return out_torque;
	//return output;
}
double ChFunctionController::ComputeOutput(double t) {
	////////////////////////////////
	//if position was changed resetCumError flag is set to true

	double curr_ang = controller_->GetAngle(index_, t);
	double exp_ang = controller_->GetExpAngle(index_, t);
	double des_ang = controller_->GetDesiredAngle(index_, t); ///get the next angle

	des_ang = controller_->LinearInterpolate(index_, curr_ang, des_ang); //linear interpolate for situations where gui changes so there isn't a major speed increase

	double error = des_ang - curr_ang;
	double prevError = controller_->prevError_.at(index_);
	double s=1;
	double K = p_gain;
	double Ti = i_gain;
	double Td = d_gain;
	double h = dT;
	double Tt = 1*dT;//read about tt
	double N = 10; //N=[8-20] http://www.cds.caltech.edu/~murray/courses/cds101/fa02/caltech/astrom-ch6.pdf
	double b = 1;

	double ulimLow = this->controller_->smarticle_->angLow;
	double ulimHigh = this->controller_->smarticle_->angHigh;
	double vlim = controller_->omegaLimit;
	double tlim = controller_->outputLimit;

	double bi = K*h / Ti;//integral gain
	double ad = (2 * Td - N*h) / (2 * Td + N*h);
	double bd = 2 * K*N*Td / (2 * Td + N*h); //deriv gain
	double ao = h / Tt;
	double ysp = des_ang;
	double y = curr_ang;



	double pp = K*(b*ysp - y);
	controller_->DD[index_] = ad*controller_->DD[index_] - bd*(y - controller_->yold[index_]);
	double v = pp + controller_->II[index_] + controller_->DD[index_];
	double u = SaturateValue(v, tlim);
	//double u = v;
	controller_->II[index_] = controller_->II[index_] + bi*(ysp - y) + ao*(u - v);
	double vel = (controller_->yold[index_] - y) / h;
	
	//GetLog() << "\nvel:" << des_ang;
	vel = SaturateValue(vel, vlim);
	double tor = (controller_->velOld[index_] - vel) / h;
	tor = SaturateValue(tor, tlim);
	double out = u;

	//GetLog() << "\nouts:" << u << "\t" << vel << "\t" << out;
	controller_->yold[index_] = y;
	controller_->velOld[index_] = vel;

	if (controller_->resetCumError)
		ResetCumulative(t);

	////////////////////////////////


	//double	divScale = 20;
	//if (stapleSize)
	//{
	//	divScale = 24; //(mass of servosize/mass of stapleSize with sizescale=5)
	//}
	//double p = p_gain / divScale;
	//double i = i_gain / divScale;
	//double d = d_gain / divScale;

	//double curr_ang = controller_->GetAngle(index_, t);
	//double exp_ang = controller_->GetExpAngle(index_, t);
	//double des_ang = controller_->GetDesiredAngle(index_, t); ///get the next angle

	//des_ang = controller_->LinearInterpolate(index_, curr_ang, des_ang); //linear interpolate for situations where gui changes so there isn't a major speed increase

	//double error = des_ang - curr_ang;
	//double prevError = controller_->prevError_.at(index_);


	////if position was changed resetCumError flag is set to true
	//if (controller_->resetCumError)
	//	ResetCumulative(t);

	//controller_->cumError_.at(index_) += (error)*dT;
	//
	//double pt = p *(error);
	//double it = i*dT*controller_->cumError_.at(index_);
	//
	//
	//double output = p*error + d*((error - prevError) / dT) + i*controller_->cumError_.at(index_);

	////set current error to previous error
	//controller_->prevError_.at(index_) = error;




	////////////////////
	////Speed control
	//double curr_omeg;
	//double des_omeg;
	//double omError;
	//double prevOmError;

	//bool deadBandOff = true;
	//double omLim = controller_->omegaLimit;
	//double omErrCond = omLim*dT * 4;
	//bool cond = (abs(error) > omErrCond); //1 if error is less than 4 omegaLimit timesteps
	//curr_omeg = controller_->GetActuatorOmega(index_, t);
	//des_omeg = omLim*cond*sgn(error);
	//if (!cond && Smarticle::global_GUI_value!=0) //if (abs(curr_ang) < omErrCond*2 || !cond)
	//{
	//	deadBandOff = false;
	//	
	//}

	//omError = des_omeg - curr_omeg;
	//prevOmError = controller_->prevOmegError_.at(index_);
	//controller_->cumOmegError_.at(index_) += (omError)*dT;
	//controller_->prevOmegError_.at(index_) = omError;
	//

	//double pTerm = 35/divScale*omLim*dT*omError; //30
	//double iTerm = 2/divScale*controller_->cumOmegError_.at(index_);//.75
	//double dTerm = 0.00052/divScale*omLim*dT*((omError - prevOmError) / dT);//.0013
	////double dTerm = .01 * dT*(omError - prevOmError / dT);


	//double output2 = pTerm + dTerm +iTerm;
	//
	//double output3 = output + output2*deadBandOff;

	//// GetLog() << "1: " << omError << " \t2: " << output2 << "\t3: " << curr_omeg << "\tom: " << curr_omeg << "\n";
	////////////////////







	//return output3;
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
	//GetLog() << "dt="<<t<<"\n";
	bothArmsReset++; //add one to value, if value>1, both arms have been reset thus value can be set to false
	controller_->II[index_] = 0;
	if (bothArmsReset >1)
	{
		//GetLog() << "reset\n";
		controller_->resetCumError = false;
		bothArmsReset = 0;
	}
}