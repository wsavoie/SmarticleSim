#include "Controller.h"
#include "ChFunction_Controller.h"
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
	desiredAngle_.assign(len, 0);
	currAngle_.assign(len, 0);
	currOmega_.assign(len, 0);
	currTorque_.assign(len, 0);
	/*contact_reporter_ = new ExtractContactForce(&contact_force_list_);
	for (size_t i = 0; i < amplitudes_.rows(); ++i) {
		amplitudes_(i) = default_amplitude_;
	}*/
	
	
}
Controller::~Controller()
{
	
	smarticle_->~Smarticle();
	ch_system_->~ChSystem();
	//smarticle_->~Smarticle();
}
bool Controller::Step(double dt) {
	//ProcessCommandQueue(dt);
	bool result = false;
	steps_++;
	if (smarticle_->GetArm0OT() || smarticle_->GetArm2OT())
	{
		result = false; //leaving this heere because maybe want to make it true
	}


	//TODO omega1prev==0 successful motion for gui1...?
	double ang01 = smarticle_->GetCurrAngle(0);
	double ang12 = smarticle_->GetCurrAngle(1);

	bool cannotMoveNext0 = smarticle_->CanMoveToNextIdx(0, ang01); //bad method name
	bool cannotMoveNext1 = smarticle_->CanMoveToNextIdx(1, ang12);

	//set desired angle to current angle if not able to move to next idx yet
	if (cannotMoveNext0)
		SetDesiredAngle(0,ang01);
	if (cannotMoveNext1)
		SetDesiredAngle(0, ang12);
	
	if(~cannotMoveNext0 && ~cannotMoveNext1)
	{
		GetDesiredAngularSpeed(0, dt);
		GetDesiredAngularSpeed(1, dt);
		result = true;
	}	

	//TODO fix choose omega amount 
	UseForceControl();
	//UseSpeedControl();

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
double Controller::GetCurrAngularSpeed(size_t index, double t)
{
	currOmega_.at(index) = smarticle_->GetOmega(index);
	return 	currOmega_.at(index);

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

void Controller::UseSpeedControl() {
	for (size_t i = 0; i < smarticle_->numEngs; ++i) {
		GetEngine(i)->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		ChSharedPtr<ChFunctionController> engine_funct(new ChFunctionController(i, this));
		GetEngine(i)->Set_spe_funct(engine_funct);
		
		//ChSharedPtr<ChFunction_Const> mfun1 = GetEngine(i)->Get_spe_funct().DynamicCastTo<ChFunction_Const>();
		//mfun1->Set_yconst(desiredOmega_.at(i));
	}
}
void Controller::UseForceControl() {
	
	for (size_t i = 0; i < smarticle_->numEngs; ++i) {
		GetEngine(i)->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
		ChSharedPtr<ChFunctionController> engine_funct(new ChFunctionController(i, this));
		GetEngine(i)->Set_tor_funct(engine_funct);
	}
}
