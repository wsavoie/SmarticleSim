//#include "chfunction_controller.h"
#include "ChFunction_Controller.h"
#include "Smarticle.h"
#include "Controller.h"
using namespace chrono;
double ChFunctionController::Get_y(double t) const { return 0; }
double ChFunctionController::Get_y(double t) {
	double curr_react_torque = controller_->GetCurrReactTorque(index_, t);
	double ct = controller_->smarticle_->getLinkActuator(index_)->Get_mot_torque();
	//add the torque already being place on the body to the torque for the next step
	//double out_torque =output; //add the torque already being placed on the body to the torque for the next step

	double output = 0;
	switch (this->controller_->smarticle_->getLinkActuator(index_)->Get_eng_mode())
	{
	case ChLinkEngine::ENG_MODE_ROTATION:
		output = ComputeOutputPosition(t);
		output = SaturateValue(output, controller_->outputLimit);
		break;
	case ChLinkEngine::ENG_MODE_SPEED:
		output = ComputeOutputSpeed(t);

		//output = OmegaToTorque(t, controller_->GetDesiredAngle(index_, t));
		
		
		//output = SaturateValue(output, controller_->omegaLimit);
		//GetLog() << this->controller_->smarticle_->getLinkActuator(index_)->Get_react_torque() << nl;
		break;
	case ChLinkEngine::ENG_MODE_TORQUE:
		output = ComputeOutputTorque(t);
		output = SaturateValue(output, controller_->outputLimit);
		break;
	default:
		output = ComputeOutputTorque(t);
		output = SaturateValue(output, controller_->outputLimit);
		break;
	}
	return output;
}

double ChFunctionController::ComputeOutputTorque(double curr_t)//http://robotics.stackexchange.com/questions/1032/how-is-piv-control-performed
{
	double Out = 0;
	controller_->posOld[index_] = controller_->posCur[index_];
	controller_->posCur[index_] = controller_->GetAngle(index_, curr_t);

	controller_->velOld[index_] = controller_->velCur[index_];
	controller_->velCur[index_] = controller_->GetActuatorOmega(index_, curr_t);

	controller_->torOld[index_] = controller_->torCur[index_];
	controller_->torCur[index_] = controller_->GetCurrTorque(index_, curr_t);
	double curr_react_torque = controller_->GetCurrReactTorque(index_, curr_t);
	double l = controller_->smarticle_->l*controller_->smarticle_->GetArm(index_ * 2)->GetMass();
	//double curr_react_torque1 = controller_->smarticle_->getLinkActuator(index_)->Get_mot_retorque()*l;
	
	//chlinkengine
	// double Get_mot_retorque() const { return mot_retorque; }
	// double Get_mot_torque() const { return mot_torque; }

	//linkmarkers
	//const ChVector<>& GetC_torque() const { return C_torque; }
	
	//linkmasked
	//ChLinkForce* GetForce_Rz() { return force_Rz; }

	// we will use the internal information of the mEngine to compute the error;
	// then we can use the error to compute the torque needed
	double desiredRotation = controller_->GetDesiredAngle(index_, curr_t); ///get the next angle

	desiredRotation = controller_->LinearInterpolate(index_, controller_->posCur[index_], desiredRotation); //linear interpolate for situations where gui changes so there isn't a major speed increase
	double currError = desiredRotation - controller_->posCur[index_];

	double desiredRotation_dt = sgn(desiredRotation - controller_->posCur[index_])*controller_->velCur[index_];
	//GetLog() << controller_->velCur[index_] << "\n";
	//double currError_dt = desiredRotation_dt - controller_->posCur[index_];
	double currError_dt = currError - desiredRotation_dt;

	double desiredRotation_dt_dt = curr_react_torque;
	double currError_dt_dt = desiredRotation_dt_dt - controller_->torCur[index_];

	double prevError = controller_->mLastError[index_];

	double p = p_gain*currError;
	double i = controller_->mAccuError[index_] += i_gain*.5*(currError + prevError)*dT;
	double dd = d_gain*(currError - prevError) / dT; //smooth out instabilities moving average of
	double d = derivAvg(dd);
	controller_->mAccuError[index_] += i_gain*.5*(currError + prevError)*dT;
	//double Out = p_gain * currError + i_gain * controller_->mAccuError[index_] + d_gain*(currError-prevError)/dT;

	double kt = .03; //motor's torque constant
	double J = 6e-6;//inertia of unloaded arm
	double gs = 1;
	double b = 1e-9;	
	

	double out = (p + i + d);
	//double b = 5.5e-4;

	out = SaturateValue(out, controller_->outputLimit) - controller_->velCur[index_] * b;
	controller_->mLastError[index_] = currError;
	controller_->mLastCalled[index_] = curr_t;
	controller_->mLastValue[index_] = out;
	CheckReset();

	return out;

}
double ChFunctionController::ComputeOutputPosition(double curr_t)//%%%%%%%%%%%%%%%%%%%%TINGNANS CODE%%%%%%%%%%%%%%%%%%%%%%%%%
{
	controller_->posOld[index_] = controller_->posCur[index_];
	controller_->posCur[index_] = controller_->GetAngle(index_, curr_t);

	controller_->velOld[index_] = controller_->velCur[index_];
	controller_->velCur[index_] = controller_->GetActuatorOmega(index_, curr_t);

	controller_->torOld[index_] = controller_->torCur[index_];
	controller_->torCur[index_] = controller_->GetCurrTorque(index_, curr_t);
	double curr_react_torque = controller_->GetCurrReactTorque(index_, curr_t);


	// we will use the internal information of the mEngine to compute the error;
	// then we can use the error to compute the torque needed
	double desiredRotation = controller_->GetDesiredAngle(index_, curr_t); ///get the next angle

	desiredRotation = controller_->LinearInterpolate(index_, controller_->posCur[index_], desiredRotation); //linear interpolate for situations where gui changes so there isn't a major speed increase
	double currError = desiredRotation - controller_->posCur[index_];

	double desiredRotation_dt = sgn(desiredRotation - controller_->posCur[index_])*controller_->omegaLimit;
	double currError_dt = desiredRotation_dt - controller_->posCur[index_];

	double desiredRotation_dt_dt = curr_react_torque;
	double currError_dt_dt = desiredRotation_dt_dt - controller_->torCur[index_];

	double prevError = controller_->mLastError[index_];


	//mLastError.assign(len, 0);
	//// accumulated input
	//mAccuError.assign(len, 0);
	//// last time the function is called
	//mLastCalled.assign(len, 0);
	//// last output
	//mLastValue.assign(len, 0);
	//// trapezoidal method
	controller_->mAccuError[index_] += .5*(currError + prevError) * dT;

	double Out = p_gain * currError + i_gain * controller_->mAccuError[index_] + d_gain * (currError - prevError) / dT;


	Out = SaturateValue(Out, controller_->outputLimit);
	controller_->mLastError[index_] = currError;
	controller_->mLastCalled[index_] = curr_t;
	controller_->mLastValue[index_] = Out;
	CheckReset();

	controller_->desPrev[index_] = controller_->GetDesiredAngle(index_, curr_t);
	// std::cout << mAccuError << std::endl;
	return Out;
}
//double ChFunctionController::ComputeOutputSpeed(double t)
//{
//	double curr_ang = controller_->GetAngle(index_, t);
//	double exp_ang = controller_->GetExpAngle(index_, t);
//	double des_ang = controller_->GetDesiredAngle(index_, t); ///get the next angle
//	des_ang = controller_->LinearInterpolate(index_, curr_ang, des_ang); //linear interpolate for situations where gui changes so there isn't a major speed increase
//	double error = des_ang - curr_ang;
//
//
//
//	double prevError = controller_->prevError_.at(index_);
//
//	double K = p_gain;
//	double Ti = i_gain;
//	double Td = d_gain;
//	double Tt = 1 * dT;//read about tt
//	double N = 10; //N=[8-20] http://www.cds.caltech.edu/~murray/courses/cds101/fa02/caltech/astrom-ch6.pdf
//	double b = 1;
//
//	double ulimLow = this->controller_->smarticle_->angLow;
//	double ulimHigh = this->controller_->smarticle_->angHigh;
//	double vlim = controller_->omegaLimit;
//
//	double bi = K*dT / Ti;//integral gain
//	double ad = (2 * Td - N*dT) / (2 * Td + N*dT);
//	double bd = 2 * K*N*Td / (2 * Td + N*dT); //deriv gain
//	double ao = dT / Tt;
//	double ysp = des_ang;
//	double y = curr_ang;
//
//	//initializes yold to current value for first iteration
//	if (controller_->smarticle_->steps == 0)
//		controller_->yold[index_] = y;
//
//	double pp = K*(b*ysp - y);
//	controller_->DD[index_] = ad*controller_->DD[index_] - bd*(y - controller_->yold[index_]);
//	double v = pp + controller_->II[index_] + controller_->DD[index_];
//	double u = SaturateValue(v, vlim);
//
//	controller_->II[index_] = controller_->II[index_] + bi*(ysp - y) + ao*(u - v);
//	controller_->yold[index_] = y;
//
//
//	//double vel = (controller_->yold[index_] - y) / dT;
//	//vel = SaturateValue(vel, vlim);
//	//double tor = (controller_->velOld[index_] - vel) / dT;
//	//tor = SaturateValue(tor, tlim);
//	//controller_->velOld[index_] = vel;
//
//	//GetLog() << u<<nl;
//
//	double out = y;
//
//	CheckReset();
//
//	return out;
//
//}
double ChFunctionController::ComputeOutputSpeed(double t)
{
	double curr_ang = controller_->GetAngle(index_, t);
	double exp_ang = controller_->GetExpAngle(index_, t);
	double des_ang = controller_->GetDesiredAngle(index_, t); ///get the next angle
	des_ang = controller_->LinearInterpolate(index_, curr_ang, des_ang); //linear interpolate for situations where gui changes so there isn't a major speed increase

	//strength in num of smarts
	//GetLog() << "&&&&&&&&&&&&&&"<< this->controller_->smartStr << nl;
	//double maxTor = this->controller_->smartStr * this->controller_->smarticle_->GetMass()*this->controller_->smarticle_->w;
	
	double g=abs(this->controller_->ch_system_->Get_G_acc().z());
	double maxTor = this->controller_->smartStr*g*this->controller_->smarticle_->refSmartMass*this->controller_->smarticle_->w; //based torque off of strength variable* weight of 0.7 l/w smarticle *distance of w;
	double currT=abs(this->controller_->smarticle_->getLinkActuator(index_)->Get_mot_torque());
	double error = des_ang - curr_ang;
	int countMax = 20;

	if (currT > maxTor)
	{

		this->controller_->smarticle_->countOT++;
		if (this->controller_->smarticle_->countOT< countMax)
		{
			error = error / 2;
		}
		if (this->controller_->smarticle_->countOT >= countMax)
		{
			error = 0;
			//GetLog() << "Full counter OT!" << nl;
			controller_->smarticle_->textureAssets[index_ * 2]->SetTextureFilename(controller_->smarticle_->mtextureOT->GetTextureFilename());
			controller_->smarticle_->ChangeToOT(true,index_);
		}
	}
	else
	{
		this->controller_->smarticle_->countOT--;
		if (this->controller_->smarticle_->countOT >= countMax || controller_->smarticle_->arm0OT || controller_->smarticle_->arm2OT)
		{

			controller_->smarticle_->textureAssets[index_ * 2]->SetTextureFilename(controller_->smarticle_->mtextureArm->GetTextureFilename());
			controller_->smarticle_->ChangeToOT(false,index_);
		}
		this->controller_->smarticle_->countOT = 0;
	}

	return error/dT;
}

double ChFunctionController::OutputToOmega(double t, double out) {
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
void ChFunctionController::ResetCumulative(double t = 0)
{

	static int bothArmsReset = 0;
	controller_->cumError_.at(index_) = 0;
	controller_->cumOmegError_.at(index_) = 0;
	controller_->mAccuError[index_] = 0;
	controller_->II[index_] = 0;
	bothArmsReset++; //add one to value, if value>1, both arms have been reset thus value can be set to false
	controller_->II[index_] = 0;
	controller_->mAccuError[index_] = 0;


	if (index_ == 0)
		controller_->dAvg0_.clear();
	else
		controller_->dAvg1_.clear();

	if (bothArmsReset > 1)
	{
		controller_->resetCumError = false;
		bothArmsReset = 0;
	}
}
void ChFunctionController::CheckReset()
{
	double nextAng = controller_->smarticle_->GetNextAngle(index_);
	//resetCumError happens on gait change
	if (controller_->resetCumError)
	{
		//GetLog() << "reset after gait change\n";
		ResetCumulative();
		controller_->prevAngle[index_] = nextAng;
		return;
	}
	//if current directed position is different than previous one reset cumulative error as endpoint is now different

	if (abs(controller_->velCur[index_]) > controller_->omegaLimit)
	{
		ResetCumulative();
		controller_->prevAngle[index_] = nextAng;
		return;
	}

	if (nextAng != controller_->prevAngle[index_])
	{
		//GetLog() << "reset after new point in traj\n";
		ResetCumulative();
	}
	controller_->prevAngle[index_] = nextAng;

}

double ChFunctionController::derivAvg(double newD)
{
	double derivAvg = 0;
	int steps = controller_->smarticle_->steps;
	std::deque<double> *b;
	if (index_ == 0)
		b = &controller_->dAvg0_;
	else
		b = &controller_->dAvg1_;
	b->emplace_front(newD);
	int count = (int)b->size();
	if (b->size() > controller_->maxDAvgSize)
	{
		b->pop_back();
	}
	count = std::min(controller_->maxDAvgSize, (int)b->size());
	for (int i = 0; i < count; i++)
	{
		derivAvg += b->at(i) / count;
		//GetLog() << "idx: " << i << " val:"<< b->at(i)<< "avg: "<< derivAvg <<"\n";
	}
	return derivAvg;
}