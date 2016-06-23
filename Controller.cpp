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
	int len = smarticle_->numEngs;
	//desiredOmega_.assign(len, 0);
	//desiredAngle_.assign(len, 0);
	//currAngle_.assign(len, 0);
	//currOmega_.assign(len, 0);
	//currTorque_.assign(len, 0);
	prevError_.assign(len, 0);
	prevOmegError_.assign(len, 0);
	cumError_.assign(len, 0);
	cumOmegError_.assign(len, 0);
	successfulMove_.assign(len, false);
	prevAngle_.assign(len, 0);

	mLastError.assign(len,0);
	// accumulated input
	mAccuError.assign(len, 0);
	// last time the function is called
	mLastCalled.assign(len, 0);
	// last output
	mLastValue.assign(len, 0);



	//have to define this way because VS is not completely C++11 compliant
	velOld[0] = 0;
	velOld[1] = 0;

	DD[0] = 0;
	DD[1] = 0;

	II[0] = 0;
	II[1] = 0;

	yold[0] = 0;
	yold[1] = 0;


	//engine_funct0 = std::make_shared<ChFunctionController>(0, this);
	//engine_funct1 = std::make_shared<ChFunctionController>(1, this);
	/*contact_reporter_ = new ExtractContactForce(&contact_force_list_);
	for (size_t i = 0; i < amplitudes_.rows(); ++i) {
		amplitudes_(i) = default_amplitude_;
	}*/
	
	
}
Controller::~Controller()
{
	//desiredOmega_.clear();
	//desiredAngle_.clear();
	//currAngle_.clear();
	//currOmega_.clear();
	//currTorque_.clear();

	prevError_.~vector();
	prevOmegError_.~vector();
	cumError_.~vector();
	cumOmegError_.~vector();
	prevAngle_.~vector();
	successfulMove_.~vector();
	//engine_funct0->~ChFunctionController();
	//engine_funct1->~ChFunctionController();

}
bool Controller::Step(double dt) {
	bool result = false;
	double ang = 0;

	for (size_t i = 0; i < smarticle_->numEngs; i++)
	{
		ang = smarticle_->GetCurrAngle(i);
		double exp = smarticle_->GetExpAngle(i);
		smarticle_->SetAngle(i, ang);
		bool desiredPositionStatus = smarticle_->NotAtDesiredPos(i, ang,exp); //true if not at current position being directed to 
		successfulMove_.at(i) = !desiredPositionStatus;
		UseForceControl(i);
		//UseSpeedControl(i);
		
	}
	result = successfulMove_.at(0) || successfulMove_.at(1);


	return result;
}

std::shared_ptr<ChLinkEngine> Controller::GetEngine(size_t index)
{ 
	return smarticle_->getLinkActuator(index); 
}

double Controller::GetCurrTorque(size_t index, double t)
{
	//return smarticle_->GetZReactTorque(index);
	return smarticle_->GetReactTorqueLen(index);
}

double Controller::GetDesiredAngle(size_t index, double t)
{
		return 	smarticle_->GetNextAngle(index);
}
double Controller::GetExpAngle(size_t idx, double t)
{
	return smarticle_->GetExpAngle(idx);
}

double Controller::GetAngle(size_t index, double t)
{
	return smarticle_->GetCurrAngle(index);
}

double Controller::GetActuatorOmega(size_t index, double t)
{
	return smarticle_->GetActuatorOmega(index);
}

double Controller::LinearInterpolate(size_t idx, double curr, double des)
{
	double errLim = 2.75 * D2R; // THIS NUMBER MATCHES EXPERIMENT IF 
	double err = (des - curr);
	err = SaturateValue(err, errLim);

	return err + curr;
}

double Controller::OmegaLimiter(size_t idx, double omega)
{
	return SaturateValue(omega, omegaLimit);	
}


bool Controller::OT()
{
	if (smarticle_->GetArm0OT() || smarticle_->GetArm2OT())
		return true;
	else
		return false;
}

void Controller::UseSpeedControl(size_t id) {
		//GetEngine(id)->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		auto ef = std::make_shared<ChFunctionController>(id, this);
		auto mfun = std::dynamic_pointer_cast<ChFunction_Const>(GetEngine(id)->Get_spe_funct());
		double y = ef->Get_y(ch_system_->GetChTime());
		mfun->Set_yconst(y);
}
void Controller::UseForceControl(size_t id) {
		//GetEngine(id)->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
		auto ef = std::make_shared<ChFunctionController>(id, this);
		auto mfun = std::dynamic_pointer_cast<ChFunction_Const>(GetEngine(id)->Get_tor_funct());
		double y = ef->Get_y(ch_system_->GetChTime());
		mfun->Set_yconst(y);	
}
