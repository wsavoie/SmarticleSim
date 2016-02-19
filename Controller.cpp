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
<<<<<<< HEAD
	//engine_funct0 = std::shared_ptr<ChFunctionController>(new ChFunctionController(0, this));
	//engine_funct1 = std::shared_ptr<ChFunctionController>(new ChFunctionController(1, this));
=======
	//engine_funct0 = std::make_shared<ChFunctionController>(0, this);
	//engine_funct1 = std::make_shared<ChFunctionController>(1, this);
>>>>>>> develop
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
	successfulMove_.~vector();
	prevAngle_.~vector();


	//engine_funct0->~ChFunctionController();
	//engine_funct1->~ChFunctionController();

}
bool Controller::Step(double dt) {
	bool result = false;
	steps_++;


	for (size_t i = 0; i < smarticle_->numEngs; i++)
	{
		double ang = smarticle_->GetCurrAngle(i);
		double exp = smarticle_->GetExpAngle(i);
		//int ind = smarticle_->moveTypeIdxs.at(smarticle_->moveType);
		//GetLog() << " id " << i << ": C:" << ang << " X:"<<exp ;
		smarticle_->SetAngle(i, ang,false);
		
		//GetLog() << "ANG" << i << ":" << smarticle_->GetAngle(i) << "    ";
		bool desiredPositionStatus = smarticle_->NotAtDesiredPos(i, ang,exp); //true if not at current position being directed to 
		//auto a = ~desiredPositionStatus;
		successfulMove_.at(i) = !desiredPositionStatus;
		UseForceControl(i);
	}
	//GetLog() <<" t: " << dt << " " << "\n";
	//exit(-1);
	
	//UseSpeedControl()

	result = successfulMove_.at(0) || successfulMove_.at(1);

	//if (successfulMove_.at(0) == false && successfulMove_.at(1) == false)
	//{
	//	result = false;
	//}
	//else{
	//	result = true;
	//}

	return result;
}

<<<<<<< HEAD
std::shared_ptr<ChLinkEngine> Controller::GetEngine(size_t index) 
=======
std::shared_ptr<ChLinkEngine> Controller::GetEngine(size_t index)
>>>>>>> develop
{ 
	return smarticle_->getLinkActuator(index); 
}

double Controller::GetCurrTorque(size_t index, double t)
{
	return smarticle_->GetZReactTorque(index);
}
//double Controller::SetGetCurrAngle(size_t index, double t)
//{
//	currAngle_.at(index) = smarticle_->GetCurrAngle(index);
//	return currAngle_.at(index);
//}

//double Controller::GetCurrAngle(size_t index, double t) 
//{
//	return currAngle_.at(index);
//}
double Controller::GetDesiredAngle(size_t index, double t)
{
	//desiredAngle_.at(index) = 
		return 	smarticle_->GetNextAngle(index);
}
double Controller::GetExpAngle(size_t idx, double t)
{
	return smarticle_->GetExpAngle(idx);
}
//void Controller::SetDesiredAngle(size_t index, double desang)
//{
//	desiredAngle_.at(index) = desang;
//}
//void Controller::SetDesiredAngularSpeed(size_t index, double desOmeg)
//{
//	desiredOmega_.at(index) = desOmeg;
//}
double Controller::GetAngle(size_t index, double t)
{
	return smarticle_->GetCurrAngle(index);
}
//void Controller::SetCurrAngle(size_t index, double ang)
//{
//	currAngle_.at(index) = ang;
//}
//void Controller::CalcCurr_Omega(size_t index, double t)
//{
//	currOmega_.at(index) = smarticle_->GetActuatorOmega(index);
//}
double Controller::GetActuatorOmega(size_t index, double t)
{
	return smarticle_->GetActuatorOmega(index);
}
//double Controller::GetCurrOmega(size_t index, double t)
//{
//	return 	currOmega_.at(index);
//}
double Controller::LinearInterpolate(size_t idx, double curr, double des)
{
	double errLim = 10 * D2R;
	double err = (des - curr);
	err = std::max(std::min(errLim, err), -errLim);

	return err + curr;

	//double om = fabs(err / dT);
	//if (om > omegaLimit)
	//{
	//	double newVal = omegaLimit*dT*sign(err) + curr;
	//	//SetDesiredAngle(idx, newVal);
	//	//smarticle_->SetNextAngle(idx, newVal);
	//	return newVal;
	//}
	//return des;
}

double Controller::OmegaLimiter(size_t idx, double omega)
{
	return std::max(std::min(omegaLimit, omega), -omegaLimit);
	
}
//double Controller::GetDesiredAngularSpeedForFunction(size_t index, double t)
//{
//	return desiredOmega_.at(index);
//}
//double Controller::GetDesiredAngularSpeed(size_t index, double t)
//{
//	desiredOmega_.at(index) = smarticle_->GetNextOmega(index); 
//	return desiredOmega_.at(index);
//}
//
//double Controller::GetDesiredAngularSpeed2(size_t index, double t,double error)
//{
//	double J = this->GetCurrTorque(index, t);
//	double tau = J / friction;
//	double s = currOmega_.at(index);
//	double noK = error / (s*(tau*s + 1));
//	desiredOmega_.at(index) = noK*(p_gain + d_gain*s);
//	return desiredOmega_.at(index);
//}

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
<<<<<<< HEAD
		std::shared_ptr<ChFunctionController> engine_funct(new ChFunctionController(i, this)); 
		GetEngine(i)->Set_spe_funct(engine_funct);
		
		//std::shared_ptr<ChFunction_Const> mfun1 = GetEngine(i)->Get_spe_funct().DynamicCastTo<ChFunction_Const>();
=======
		auto engine_funct = std::make_shared<ChFunctionController>(i, this);
		GetEngine(i)->Set_spe_funct(engine_funct);
		
		//auto mfun1 = std::dynamic_pointer_cast<ChFunction_Const>(GetEngine(i)->Get_spe_funct());
>>>>>>> develop
		//mfun1->Set_yconst(desiredOmega_.at(i));
	}
	GetEngine(0)->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	GetEngine(1)->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	//GetEngine(0)->Set_spe_funct(engine_funct0);
	//GetEngine(1)->Set_spe_funct(engine_funct1);

}
void Controller::UseForceControl(size_t id) {
	//for (size_t i = 0; i < smarticle_->numEngs; ++i) {
		
		GetEngine(id)->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
<<<<<<< HEAD
		std::shared_ptr<ChFunctionController> ef = std::make_shared<ChFunctionController>(id, this);
		std::shared_ptr<ChFunction_Const> mfun = std::dynamic_pointer_cast<ChFunction_Const>(GetEngine(id)->Get_tor_funct());
=======
		auto ef = std::make_shared<ChFunctionController>(id, this);
		auto mfun = std::dynamic_pointer_cast<ChFunction_Const>(GetEngine(id)->Get_tor_funct());
>>>>>>> develop
		double y = ef->Get_y(ch_system_->GetChTime());
		mfun->Set_yconst(y);
	
		//GetEngine(0)->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
		//GetEngine(1)->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);

<<<<<<< HEAD
		////std::shared_ptr<ChFunctionController> engine_funct0(new ChFunctionController(0, this));
		//std::shared_ptr<ChFunctionController> ef0(new ChFunctionController(0, this));
		//std::shared_ptr<ChFunctionController> ef1(new ChFunctionController(1, this));
=======
		////auto engine_funct0 = std::make_shared<ChFunctionController>(0, this);
		//auto ef0 = std::make_shared<ChFunctionController>(0, this);
		//auto ef1 = std::make_shared<ChFunctionController>(1, this);
>>>>>>> develop
		//
		////causes 2 calls per chfunction
		////GetEngine(0)->Set_tor_funct(ef0);
		////GetEngine(1)->Set_tor_funct(ef1);

		////doing this fixes the multiple chfunction runs
<<<<<<< HEAD
		//std::shared_ptr<ChFunction_Const> mfun0 = GetEngine(0)->Get_tor_funct().DynamicCastTo<ChFunction_Const>();
		//std::shared_ptr<ChFunction_Const> mfun1 = GetEngine(1)->Get_tor_funct().DynamicCastTo<ChFunction_Const>();
=======
		//auto mfun0 = std::dynamic_pointer_cast<ChFunction_Const>(GetEngine(0)->Get_tor_funct());
		//auto mfun1 = std::dynamic_pointer_cast<ChFunction_Const>(GetEngine(1)->Get_tor_funct());
>>>>>>> develop
		//double y0 = ef0->Get_y(ch_system_->GetChTime());
		//double y1 = ef1->Get_y(ch_system_->GetChTime());
		////GetLog() << "y0: " << y0 << "y1: " << y1 << "\n";
		//mfun0->Set_yconst(y0);
		//mfun1->Set_yconst(y1);

	
		//smarticle_->~Smarticle();

	//}

}
