#include "Controller.h"
//#include "ChFunction_Controller.h"
#include "Smarticle.h"
#include <utility>
#include <tuple>
//#include "include/vector_utility.h"
//#include "include/contact_reporter.h"


using namespace chrono;
Controller::Controller(chrono::ChSystem *ch_system, Smarticle *smarticle)
	: ch_system_(ch_system), smarticle_(smarticle)
{
	int len = 2;
	desiredOmega_.assign(len, 0);
	successfulMove_.assign(len, false);
	desiredAngle_.assign(len, 0);
	currAngle_.assign(len, 0);
	currOmega_.assign(len, 0);
	currTorque_.assign(len, 0);
	engine_funct0 = ChSharedPtr<ChFunctionController>(new ChFunctionController(0, this));
	engine_funct1 = ChSharedPtr<ChFunctionController>(new ChFunctionController(1, this));
	/*contact_reporter_ = new ExtractContactForce(&contact_force_list_);
	for (size_t i = 0; i < amplitudes_.rows(); ++i) {
		amplitudes_(i) = default_amplitude_;
	}*/
	
	
}
Controller::~Controller()
{
	desiredOmega_.clear();
	desiredAngle_.clear();
	currAngle_.clear();
	currOmega_.clear();
	currTorque_.clear();

	engine_funct0->~ChFunctionController();
	engine_funct1->~ChFunctionController();

}
bool Controller::Step(double dt) {
	bool result = false;
	steps_++;


	for (size_t i = 0; i < smarticle_->numEngs; i++)
	{
		double ang = smarticle_->GetCurrAngle(i);
		smarticle_->SetAngle(i, ang);
		bool cannotMoveNext = smarticle_->CanMoveToNextIdx(i, ang); //bad method name 
		successfulMove_.at(i) = cannotMoveNext;
	}

	result = UseForceControl();
	//UseSpeedControl();

	//if (smarticle_->GetArm0OT() || smarticle_->GetArm2OT())
	//{
	//	result = false; //leaving this heere because maybe want to make it true
	//}

	return result;
}

ChSharedPtr<ChLinkEngine> Controller::GetEngine(size_t index) 
{ 
	return smarticle_->getLinkActuator(index); 
}

double Controller::GetCurrTorque(size_t index, double t)
{
	return smarticle_->GetZReactTorque(index);
}


double Controller::GetCurrAngle(size_t index, double t) 
{
	currAngle_.at(index) = smarticle_->GetCurrAngle(index);
	return currAngle_.at(index);
}
double Controller::GetDesiredAngle(size_t index, double t)
{
	desiredAngle_.at(index) = smarticle_->GetNextAngle(index);
	return 	desiredAngle_.at(index);
}
void Controller::SetDesiredAngle(size_t index, double desang)
{
	desiredAngle_.at(index) = desang;
}
void Controller::SetDesiredAngularSpeed(size_t index, double desOmeg)
{
	desiredOmega_.at(index) = desOmeg;
}
double Controller::GetCurrAngularSpeed(size_t index, double t)
{
	currOmega_.at(index) = smarticle_->GetOmega(index);
	return 	currOmega_.at(index);

}
double Controller::LinearInterpolate(size_t idx, double curr, double des)
{
	double err = (des - curr);
	double omega = fabs(err / dT);
	if (omega > omegaLimit)
	{
		double newVal = omegaLimit*dT*sign(err) + curr;
		//smarticle_->SetNextAngle(idx, newVal);
		SetDesiredAngle(idx, newVal);
		//SetDesiredAngularSpeed(idx,omegaLimit);
		successfulMove_.at(idx) = false;
		return newVal;
	}
	successfulMove_.at(idx) = true;
	if (successfulMove_.at(0) && successfulMove_.at(1))
	{
		//complete successfulmove 
		//successfulMove_.at(idx) = true;
	}
	return des;
}
double Controller::GetDesiredAngularSpeedForFunction(size_t index, double t)
{
	return desiredOmega_.at(index);
}
double Controller::GetDesiredAngularSpeed(size_t index, double t)
{
	desiredOmega_.at(index) = smarticle_->GetNextOmega(index); 
	return desiredOmega_.at(index);
}

double Controller::GetDesiredAngularSpeed2(size_t index, double t,double error)
{
	double J = this->GetCurrTorque(index, t);
	double tau = J / friction;
	double s = currOmega_.at(index);
	double noK = error / (s*(tau*s + 1));
	desiredOmega_.at(index) = noK*(p_gain + d_gain*s);
	return desiredOmega_.at(index);
}

bool Controller::OT()
{
	if (smarticle_->GetArm0OT() || smarticle_->GetArm2OT())
		return true;
	else
		return false;
}

bool Controller::UseSpeedControl() {
	bool res = true;
	for (size_t i = 0; i < smarticle_->numEngs; ++i) {
		GetEngine(i)->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		ChSharedPtr<ChFunctionController> engine_funct(new ChFunctionController(i, this));
		GetEngine(i)->Set_spe_funct(engine_funct);
		
		//ChSharedPtr<ChFunction_Const> mfun1 = GetEngine(i)->Get_spe_funct().DynamicCastTo<ChFunction_Const>();
		//mfun1->Set_yconst(desiredOmega_.at(i));
	}
	GetEngine(0)->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	GetEngine(1)->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	GetEngine(0)->Set_spe_funct(engine_funct0);
	GetEngine(1)->Set_spe_funct(engine_funct1);
	if (successfulMove_.at(0) == false || successfulMove_.at(1) == false)
	{
		res = false;
	}
	return res;
}
bool Controller::UseForceControl() {
	bool res = true;
	//for (size_t i = 0; i < smarticle_->numEngs; ++i) {
		
		GetEngine(0)->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
		GetEngine(1)->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
		
		//ChSharedPtr<ChFunctionController> engine_funct0(new ChFunctionController(0, this));
		//ChSharedPtr<ChFunctionController> engine_funct1(new ChFunctionController(1, this));
		GetEngine(0)->Set_tor_funct(engine_funct0);
		GetEngine(1)->Set_tor_funct(engine_funct1);


		if (successfulMove_.at(0) == false || successfulMove_.at(1)== false)
		{res = false;}

	//}
	return res;
}
