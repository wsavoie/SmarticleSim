//#include "chfunction_controller.h"
#include "ChFunction_Controller.h"
#include "Smarticle.h"
#include "Controller.h"
using namespace chrono;

double ChFunctionController::Get_y(double t) {
	double curr_react_torque = controller_->GetCurrTorque(index_, t);
	double ct = controller_->smarticle_->getLinkActuator(index_)->Get_mot_torque();
	//add the torque already being place on the body to the torque for the next step
	//double out_torque =output; //add the torque already being placed on the body to the torque for the next step

	double output = 0;
	switch (this->controller_->smarticle_->getLinkActuator(index_)->Get_eng_mode())
	{
		case ChLinkEngine::ENG_MODE_ROTATION:
			break;
		case ChLinkEngine::ENG_MODE_SPEED:
			output = ComputeOutputSpeed(t);
			output = SaturateValue(output, controller_->omegaLimit);
			break;
		case ChLinkEngine::ENG_MODE_TORQUE:
			output = ComputeOutput(t);
			output=SaturateValue(output, controller_->outputLimit);
			
			
			break;
		default:
			output = ComputeOutput(t);
			output = SaturateValue(output, controller_->outputLimit);
				break;
	}
	////////double out_torque = SaturateValue(ComputeOutput(t) - curr_react_torque, controller_->outputLimit);

	//double out_torque = SaturateValue(ComputeOutput(t), controller_->outputLimit);
	//double mT = controller_->smarticle_->getLinkActuator(index_)->Get_mot_rerot_dtdt()*controller_->smarticle_->l*dT;

	//	bool o = false;
	//if (abs(out_torque) == controller_->outputLimit)
	//	o = true;
	return output;
}
//%%%%%%%%%%%%%%%%%%%%TINGNANS CODE%%%%%%%%%%%%%%%%%%%%%%%%%
//double Get_y(double curr_t)
//{
//	double dt = curr_t - mLastCalled;
//	if (dt < mDt)
//		return mLastValue;
//	// we will use the internal information of the mEngine to compute the error;
//	// then we can use the error to compute the torque needed
//
//	double currRotation = mEngine->Get_mot_rot();
//	double desiredRotation = mEngine->Get_rot_funct()->Get_y(curr_t);
//	double currError = desiredRotation - currRotation;
//
//	double currRotation_dt = mEngine->Get_mot_rot_dt();
//	double desiredRotation_dt = mEngine->Get_rot_funct()->Get_y_dx(curr_t);
//	double currError_dt = desiredRotation_dt - currRotation_dt;
//
//	// trapezoidal method
//	mAccuError += 0.5 * (currError + mLastError) * dt;
//
//	double Out = P * currError + I * mAccuError + D * (currError - mLastError) / dt;
//	if (Max > Min)
//	{
//		// set the limit
//		Out = std::max(std::min(Out, Max), Min);
//	}
//
//	mLastError = currError;
//	mLastCalled = curr_t;
//	mLastValue = Out;
//	// std::cout << mAccuError << std::endl;
//	return Out;
//
//}
//%%%%%%%%%
//double ChFunctionController::Get_y(double curr_t)
//{
//	double curr_react_torque = controller_->GetCurrTorque(index_, curr_t);
//	// we will use the internal information of the mEngine to compute the error;
//	// then we can use the error to compute the torque needed
//	double currRotation = controller_->GetAngle(index_, curr_t);
//	double desiredRotation = controller_->GetDesiredAngle(index_, curr_t); ///get the next angle
//	desiredRotation = controller_->LinearInterpolate(index_, currRotation, desiredRotation); //linear interpolate for situations where gui changes so there isn't a major speed increase
//	double currError = desiredRotation - currRotation;
//
//
//	double currRotation_dt = controller_->GetActuatorOmega(index_, curr_t);
//	double desiredRotation_dt = controller_->smarticle_->getLinkActuator(index_)->Get_rot_funct()->Get_y_dx(curr_t);
//	double currError_dt = desiredRotation_dt - currRotation_dt;
//
//	double prevError = controller_->mLastError[index_];
//
//
//	//mLastError.assign(len, 0);
//	//// accumulated input
//	//mAccuError.assign(len, 0);
//	//// last time the function is called
//	//mLastCalled.assign(len, 0);
//	//// last output
//	//mLastValue.assign(len, 0);
//	//// trapezoidal method
//	controller_->mAccuError[index_] +=  .5*(currError + prevError) * dT;
//	
//	double Out = p_gain * currError + i_gain * controller_->mAccuError[index_] + d_gain * (currError - prevError) / dT;
//	Out = SaturateValue(Out, controller_->outputLimit);
//	controller_->mLastError[index_]= currError;
//	controller_->mLastCalled[index_] = curr_t;
//	controller_->mLastValue[index_] = Out;
//	if (controller_->resetCumError)
//		ResetCumulative(curr_t);
//	// std::cout << mAccuError << std::endl;
//	return Out;
//
//}

//%%%%%%%%%%%%%%%%%%%%TINGNANS CODE%%%%%%%%%%%%%%%%%%%%%%%%%
//	
double ChFunctionController::Get_y(double t) const { return 0; }
double ChFunctionController::ComputeOutputSpeed(double t)
{
	double curr_ang = controller_->GetAngle(index_, t);
	double exp_ang = controller_->GetExpAngle(index_, t);
	double des_ang = controller_->GetDesiredAngle(index_, t); ///get the next angle
	des_ang = controller_->LinearInterpolate(index_, curr_ang, des_ang); //linear interpolate for situations where gui changes so there isn't a major speed increase
	double error = des_ang - curr_ang;



	double prevError = controller_->prevError_.at(index_);

	double K = p_gain;
	double Ti = i_gain;
	double Td = d_gain;
	double Tt = 1 * dT;//read about tt
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
double ChFunctionController::ComputeOutput(double t) {
	
	double curr_ang = controller_->GetAngle(index_,t);

	if (controller_->smarticle_->steps == 0)  //*********************
	{
		controller_->ycurr[index_] = curr_ang;
		controller_->yold[index_] = curr_ang;
	}
	else
	{
		controller_->yold[index_] = controller_->ycurr[index_];
		controller_->ycurr[index_] = curr_ang;	
	}

	double exp_ang = controller_->GetExpAngle(index_, t);
	double des_ang = controller_->GetDesiredAngle(index_, t); ///get the next angle
	des_ang = controller_->LinearInterpolate(index_, curr_ang, des_ang); //linear interpolate for situations where gui changes so there isn't a major speed increase
	double error = des_ang - curr_ang;
	double prevError = controller_->prevError_.at(index_);
	double prevSpeedError = (error - prevError) / dT;
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
	//if (controller_->smarticle_->steps == 0)  //*********************
	//	controller_->yold[index_] = y;

	double pp = K*(b*ysp - y);
	controller_->DD[index_] = ad*controller_->DD[index_] - bd*(y - controller_->yold[index_]);
	double v = pp + controller_->II[index_] + controller_->DD[index_];
	double u = SaturateValue(v, tlim);
	controller_->II[index_] = controller_->II[index_] + bi*(ysp - y) + ao*(u - v);
	//controller_->yold[index_] = y;//************
	



	//double vel = (controller_->yold[index_] - y) / dT;
	//vel = SaturateValue(vel, vlim);
	//double tor = (controller_->velOld[index_] - vel) / dT;
	//tor = SaturateValue(tor, tlim);
	//controller_->velOld[index_] = vel;



	double out = u;
	
	if (this->controller_->smarticle_->getLinkActuator(index_)->Get_eng_mode() == ChLinkEngine::ENG_MODE_ROTATION)


	if (controller_->resetCumError)
		ResetCumulative(t);

	/*if (this->controller_->smarticle_->getLinkActuator(index_)->GetForce_D()->Get_modul_R())
		GetLog() << index_<<"modk exists!\n";
	this->controller_->smarticle_->getLinkActuator(index_)->GetForce_D()->updGet_Force(error, prevSpeedError, t);*/
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