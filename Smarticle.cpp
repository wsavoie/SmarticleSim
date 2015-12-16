/*
 * Smarticle.cpp
 *
 *  Created on: May 22, 2015
 *      Author: Arman Pazouki
 */

#include "Smarticle.h"
#include "utils/ChUtilsGeometry.h"
#include "utils/ChUtilsCreators.h"


#include <math.h>  /* asin */


//#include "physics/ChSystem.h"  // Arman: take care of this later
#include "chrono_parallel/physics/ChSystemParallel.h"
//#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"


using namespace chrono;

Smarticle::Smarticle(
		  ChSystem* otherSystem
		  ) : m_system(otherSystem) {
	smarticleID = -1;
	density = 7800;
	l = 1;
	w = 1;
	r = .05;
	r2 = .01;
	initPos = ChVector<>(0);
	rotation = QUNIT;
	jointClearance = 0;
	volume = GetVolume();
	
}

std::vector<std::pair<double, double>> Smarticle::global;
std::vector<std::pair<double, double>> Smarticle::gui1;
std::vector<std::pair<double, double>> Smarticle::gui2;
std::vector<std::pair<double, double>> Smarticle::gui3;
ChSharedPtr<ChTexture> Smarticle::mtextureOT = ChSharedPtr<ChTexture>(new ChTexture());
ChSharedPtr<ChTexture> Smarticle::mtextureArm = ChSharedPtr<ChTexture>(new ChTexture());
ChSharedPtr<ChTexture> Smarticle::mtextureMid = ChSharedPtr<ChTexture>(new ChTexture());
double Smarticle::pctActive = 1.0;
double Smarticle::distThresh;
unsigned int Smarticle::global_GUI_value;

Smarticle::~Smarticle()
{
	//this->armsController->~Controller();
	m_system->RemoveLink(link_actuator01);
	m_system->RemoveLink(link_actuator12);
	//m_system->RemoveLink(link_revolute01);
	//m_system->RemoveLink(link_revolute12);
	m_system->RemoveBody(arm0);
	m_system->RemoveBody(arm1);
	m_system->RemoveBody(arm2);
	armsController->~Controller();
	GetLog() << "Removing smarticle\n";
}



void Smarticle::Properties(
		int sID,
		double other_density,
		ChSharedPtr<ChMaterialSurface> surfaceMaterial,
		double other_envelop,
		double other_l,
		double other_w,
		double other_r,
		double other_r2,
		ChVector<> pos,
		ChQuaternion<> rot,
		double other_angle,
		double other_angle2){

	smarticleID = sID;
	density = other_density;
	mat_g = surfaceMaterial;
	l = other_l;
	w = other_w;
	r = other_r;
	r2 = other_r2;
	collisionEnvelop = other_envelop;
	initPos = pos;
	rotation = rot;
	angle1 = other_angle;
	angle2 = other_angle2;
	volume = GetVolume();


	mtextureOT->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));

	mtextureArm->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));
	mtextureMid->SetTextureFilename(GetChronoDataFile("cubetexture_blue_bordersBlue.png"));
}

void Smarticle::Properties(
		int sID,
		double other_density,
		ChSharedPtr<ChMaterialSurface> surfaceMaterial,
		double other_envelop,
		double other_l,
		double other_w,
		double other_r,
		double other_r2,
		double other_omega,
		ChVector<> pos,
		ChQuaternion<> rot,
		double other_angle,
		double other_angle2){

	Properties(sID, other_density, surfaceMaterial, other_envelop, other_l, other_w, other_r, other_r2, pos, rot, other_angle, other_angle2);
	SetDefaultOmega(other_omega);
	//SetOmega(other_omega, true);
}
/////Will added Smarticle properties///////
void Smarticle::Properties(
	int sID,
	int mdumID,
	double other_density,
	ChSharedPtr<ChMaterialSurface> surfaceMaterial,
	double other_envelop,
	double other_l,
	double other_w,
	double other_r,
	double other_r2,
	double other_omega,
	bool willVersion,
	ChVector<> pos,
	ChQuaternion<> rot,
	double other_angle,
	double other_angle2,
	double other_torThresh2,
	double other_angLow,
	double other_angHigh){

	Properties(sID, other_density, surfaceMaterial, other_envelop, other_l, other_w, other_r, other_r2, pos, rot, other_angle, other_angle2);
	SetDefaultOmega(other_omega);
	//SetOmega1(other_omega);
	//SetOmega2(other_omega);
	SetOmega(other_omega);
	moveType = MoveType::GLOBAL;
	prevMoveType = MoveType::GLOBAL;
	//moveTypeIdxs.resize(MoveType::OT, 0);
	moveTypeIdxs.resize(MoveType::OT+1, 0);
	arm0OT = false;
	arm2OT = false;
	SetAngles(other_angle, other_angle2);
	torqueThresh2 = 10000;
	angLow = 0;
	angHigh = 120;
	dumID = mdumID;
	armBroken = false;

	std::tuple<double, double,double,double> a (0.0, 0.0,0.0,0.0);
	torques = { a, a, a, a, a, a, a }; //probably a better way to do this....
	nextOmega = { 0, 0 };
	nextAngle = { 0,0};
	currTorque = { 0, 0 };
	torqueAvg = std::make_tuple(0, 0, 0, 0);
	initialAng0 = other_angle;
	initialAng1 = other_angle2;
	//avg1 = 0;
	//avg2 = 0;
	//armsController = (new Controller(m_system, this));
}
void Smarticle::updateTorqueDeque(double mtorque0, double mtorque1,double momega0, double momega1)
{
	//torque1.pop_back();
	//torque1.push_front(mtorque1);
	//
	//torque2.pop_back();
	//torque2.push_front(mtorque2);
	
	std::tuple<double, double,double,double> oldT = torques.back();
	torques.pop_back();
	torques.emplace_front(mtorque0, mtorque1,momega0,momega1);
	torques.emplace_front(mtorque0, mtorque1, momega0, momega1);
	updateTorqueAvg(oldT);
}
void Smarticle::updateTorqueDeque()
{
	std::tuple < double, double, double, double > oldT = torques.back();
	torques.pop_back();
	torques.emplace_front(getLinkActuator(0)->Get_mot_torque(), getLinkActuator(1)->Get_mot_torque(), getLinkActuator(0)->Get_mot_rot_dt(), getLinkActuator(1)->Get_mot_rot_dt());
	updateTorqueAvg(oldT);
}
//void Smarticle::updateTorqueAvg()
////in case I need to recalculate from all values
//{
//	for (std::deque<std::pair<double, double>>::iterator it = torques.begin(); it != torques.end(); ++it)
//	{
//		torqueAvg = std::make_pair(it->first + torqueAvg.first, it->second + torqueAvg.second);
//	}
//	torqueAvg = std::make_pair(torqueAvg.first / torques.size(), torqueAvg.second / torques.size());
//}
void Smarticle::updateTorqueAvg(std::tuple <double,double,double,double > oldT)
{ //already assuming it is an avg:
	size_t len = torques.size();
	std::get<0>(torqueAvg) = std::get<0>(torqueAvg) -std::get<0>(oldT) / len + std::get<0>(torques.front()) / len;
	std::get<1>(torqueAvg) = std::get<1>(torqueAvg) -std::get<1>(oldT) / len + std::get<1>(torques.front()) / len;
	std::get<2>(torqueAvg) = std::get<2>(torqueAvg) -std::get<2>(oldT) / len + std::get<2>(torques.front()) / len;
	std::get<3>(torqueAvg) = std::get<3>(torqueAvg) -std::get<3>(oldT) / len + std::get<3>(torques.front()) / len;

	//torqueAvg = std::make_pair(
	//	torqueAvg.first - oldT.first / len + torques.front().first / len,
	//	torqueAvg.second - oldT.second / len + torques.front().second / len
	//	);
}

//////////////////////////////////////////////
void Smarticle::SetDefaultOmega(double omega) {
	defaultOmega = omega;
}
void Smarticle::SetOmega(int idx, double momega, bool angularFreq)
{
	if (idx == 0)
		SetOmega1(momega, angularFreq);
	else
		SetOmega2(momega, angularFreq);
}
void Smarticle::SetOmega(double momega, bool angularFreq)
{

	if (angularFreq)
	{
		omega1 = momega;
		omega2 = momega;
		return;
	}
	omega1 = (2 * CH_C_PI)*momega;
	omega2 = (2 * CH_C_PI)*momega;
}
//void Smarticle::SetOmega(double momega1, double momega2, bool angularFreq)
//{
//	if (angularFreq)
//	{
//		omega1 = momega1;
//		omega2 = momega2;
//		return;
//	}
//	omega1 = (2 * CH_C_PI)*momega1;
//	omega2 = (2 * CH_C_PI)*momega2;
//}
void Smarticle::SetOmega1(double momega1, bool angularFreq)
{
	if (angularFreq)
	{
		omega1 = momega1;
		return;
	}
	omega1 = (2 * CH_C_PI)*momega1;
}
void Smarticle::SetOmega2(double momega2, bool angularFreq)
{
	if (angularFreq)
	{
		omega2=momega2;
		return;
	}
	omega2 = (2 * CH_C_PI)*momega2;
}
double Smarticle::GetOmega(int id, bool angularFreq)
{
	if (id == 0)
		return GetOmega1(angularFreq);
	else
		return GetOmega2(angularFreq);
}
double Smarticle::GetActuatorOmega(int id)
{
	return getLinkActuator(id)->Get_mot_rot_dt();	
}
double Smarticle::GetNextOmega(int id)
{
	nextOmega.at(id)= this->ChooseOmegaAmount(GetOmega(id), GetCurrAngle(id), GetNextAngle(id));
	return nextOmega.at(id);
}
double Smarticle::GetOmega1(bool angularFreq)
{
	if (angularFreq)
	{
		return omega1;
	}
	return omega1 / (2 * CH_C_PI);
}

double Smarticle::GetOmega2(bool angularFreq)
{
	if (angularFreq)
	{
		return omega2;
	}
	return omega2 / (2 * CH_C_PI);
}

void Smarticle::CreateArm(int armID, double len, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChVector<> gyr;  	// components gyration
	double vol;			// components volume

	ChSharedBodyPtr arm;

	vol = utils::CalcBoxVolume(ChVector<>(len/2.0, r, r2));
	gyr = utils::CalcBoxGyration(ChVector<>(len/2.0, r, r2)).Get_Diag();
	// create body, set position and rotation, add surface property, and clear/make collision model
	if (USE_PARALLEL) {
		arm = ChSharedBodyPtr(new ChBody(new collision::ChCollisionModelParallel));
	} else {
		arm = ChSharedBodyPtr(new ChBody);
	}


	//$$$$$$$$$$$
	//$$$$$$$$$$$

	ChVector<> posArm = rotation.Rotate(posRel) + initPos;
	arm->SetName("smarticle_arm");
	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);
    arm->SetCollide(true);
    arm->SetBodyFixed(false);
    arm->GetPhysicsItem()->SetIdentifier(dumID + armID);
    if (armID == 1) //this was old code from when I was fixing them to fit
			arm->SetBodyFixed(true);
    else
    	arm->SetBodyFixed(false);

	arm->SetMaterialSurface(mat_g);

	double mass = density * vol;
	//double mass = .005;//.043/3.0; //robot weight 43 grams
	arm->GetCollisionModel()->ClearModel();
	//arm->SetLimitSpeed(true);
	//arm->SetMaxSpeed(CH_C_PI * 2 * 5);
	//arm->SetMaxWvel(CH_C_PI * 2 * 5);
	//arm->ClampSpeed();
	
	if (visualize)
	{
	if (armID == 1)
		arm->AddAsset(mtextureMid);
	else
		arm->AddAsset(mtextureArm);

	}
	arm->GetCollisionModel()->SetEnvelope(collisionEnvelop);
	utils::AddBoxGeometry(arm.get_ptr(), ChVector<>(len / 2.0, r, r2), ChVector<>(0, 0, 0),QUNIT,visualize);

	arm->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2

    arm->GetCollisionModel()->BuildModel(); // this function overwrites the intertia

    // change mass and inertia property
    arm->SetMass(mass);
    arm->SetInertiaXX(mass * gyr);
    //arm->SetDensity(density);

    m_system->AddBody(arm);



	switch (armID) {
	case 0: {
		arm0 = arm;
	} break;
	case 1: {
		arm1 = arm;
	} break;
	case 2: {
		arm2 = arm;
	} break;
	default:
		std::cerr << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
		break;
	}
}
void Smarticle::CreateArm2(int armID, double len,double mr, double mr2, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChVector<> gyr;  	// components gyration
	double vol;			// components volume

	ChSharedBodyPtr arm;

	vol = utils::CalcBoxVolume(ChVector<>(len / 2.0, mr, mr2));
	gyr = utils::CalcBoxGyration(ChVector<>(len / 2.0, mr, mr2)).Get_Diag();
	// create body, set position and rotation, add surface property, and clear/make collision model
	if (USE_PARALLEL) {
		arm = ChSharedBodyPtr(new ChBody(new collision::ChCollisionModelParallel));
	}
	else {
		arm = ChSharedBodyPtr(new ChBody);
	}
	ChVector<> posArm = rotation.Rotate(posRel) + initPos;

	arm->SetName("smarticle_arm");
	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);
	arm->SetCollide(true);
	arm->SetBodyFixed(false);
	arm->GetPhysicsItem()->SetIdentifier(dumID + armID);
	if (armID == 1) //this was old code from when I was fixing them to fit
		arm->SetBodyFixed(false);
	else
		arm->SetBodyFixed(false);

	arm->SetMaterialSurface(mat_g);

	double mass = density * vol;
	//double mass = .005;//.043/3.0; //robot weight 43 grams
	arm->GetCollisionModel()->ClearModel();

	if (visualize)
	{
		if (armID == 1)
			arm->AddAsset(mtextureMid);
		else
			arm->AddAsset(mtextureArm);

	}
	arm->GetCollisionModel()->SetEnvelope(collisionEnvelop);
	utils::AddBoxGeometry(arm.get_ptr(), ChVector<>(len / 2.0, mr, mr2), ChVector<>(0, 0, 0), QUNIT, visualize);

	arm->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2

	arm->GetCollisionModel()->BuildModel(); // this function overwrites the intertia

	// change mass and inertia property
	arm->SetMass(mass);
	arm->SetInertiaXX(mass * gyr);
	//arm->SetDensity(density);

	m_system->AddBody(arm);



	switch (armID) {
	case 0: {
		arm0 = arm;
	} break;
	case 1: {
		arm1 = arm;
	} break;
	case 2: {
		arm2 = arm;
	} break;
	default:
		std::cerr << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
		break;
	}
}

ChSharedBodyPtr Smarticle::GetArm(int armID) {
	switch (armID) {
	case 0:
		return arm0;
	case 1:
		return arm1;
	case 2:
		return arm2;
	default:
		std::cerr << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
		break;
	}
	 return ChSharedBodyPtr();
}

ChSharedPtr<ChLinkLockRevolute> Smarticle::GetRevoluteJoint(int jointID) {

	switch (jointID) {
	case 0:
		return link_revolute01;
	case 1:
		return link_revolute12;
	default:
		std::cerr << "Error! smarticle can only have joints with ids from {0, 1}" << std::endl;
		break;
	}
	return ChSharedPtr<ChLinkLockRevolute>();
}
void Smarticle::SetSpeed(ChVector<> newSpeed)
{
	arm0->SetPos_dt(newSpeed);
	arm1->SetPos_dt(newSpeed);
	arm2->SetPos_dt(newSpeed);
}
void Smarticle::TransportSmarticle(ChVector<> newPosition)
{
	arm0->SetPos(arm0->GetPos() - arm1->GetPos() + newPosition);
	arm2->SetPos(arm2->GetPos() - arm1->GetPos() + newPosition);
	arm1->SetPos(newPosition);
}
void Smarticle::CreateJoints() {
	// link 1
	link_revolute01 = ChSharedPtr<ChLinkLockRevolute>(new ChLinkLockRevolute);
	link_revolute12 = ChSharedPtr<ChLinkLockRevolute>(new ChLinkLockRevolute);
	// ChVector<> pR01(-w / 2.0+r2, 0, 0);
	// ChVector<> pR12(w / 2.0-r2, 0, 0);
	ChVector<> pR01(-w / 2.0-r2, 0, 0);
	ChVector<> pR12(w / 2.0+r2, 0, 0);
	ChQuaternion<> qx = Q_from_AngAxis(CH_C_PI / 2.0, VECT_X);

	// link 1
	link_revolute01->Initialize(arm0, arm1, ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx));
	link_revolute01->SetMotion_axis(ChVector<>(0, 0, 1));
	m_system->AddLink(link_revolute01);


	// link 2
	link_revolute12->Initialize(arm1, arm2, ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx));
	link_revolute12->SetMotion_axis(ChVector<>(0, 0, 1));
	m_system->AddLink(link_revolute12);
}


void Smarticle::CreateActuators() {

	link_actuator01 = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	link_actuator12 = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);

	//current sim
	//ChVector<> pR01(-w / 2.0-r2, 0, 0);
	//ChVector<> pR12(w / 2.0+r2, 0, 0);
	
	ChVector<> pR01(-w / 2.0, 0, 0);
	ChVector<> pR12(w / 2.0, 0, 0);


	ChQuaternion<> qx = Q_from_AngAxis(CH_C_PI / 2.0, VECT_X);
	ChQuaternion<> qy1 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, GetAngle1()));
	ChQuaternion<> qy2 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, GetAngle2()));

	// link 1
	link_actuator01->Initialize(arm0, arm1, false, ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx), ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx*qy1));
	//link_actuator01->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	//link_actuator01->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
	//link_actuator01->Set_eng_mode(ChLinkEngine::ENG_MODE_TO_POWERTRAIN_SHAFT);
	link_actuator01->SetMotion_axis(ChVector<>(0, 0, 1));
	m_system->AddLink(link_actuator01);

	link_actuator12->Initialize(arm1, arm2, false, ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx), ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx*qy2));
	//link_actuator12->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	//link_actuator12->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
	link_actuator12->SetMotion_axis(ChVector<>(0, 0, 1));
	m_system->AddLink(link_actuator12);

}

void Smarticle::Create() {
	//jointClearance =1.0/3.0*r2;
	double a1 = GetAngle1();
	double a2 = GetAngle2();
	jointClearance = 0;
	double l_mod;
	//double l_mod = l + 2 * r2 - jointClearance;

	ChQuaternion<> quat0 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0*5*CH_C_PI/180, angle1, 0));
	ChQuaternion<> quat2 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, -angle2, 0));
	quat0.Normalize();
	quat2.Normalize();

	//current sim version
	//CreateArm(0, l_mod, ChVector<>(-w / 2.0 -r2- (l_mod / 2.0-r2)*cos(angle1), 0, -(l_mod / 2.0-r2)*sin(angle1)), quat0);
	//CreateArm(1, w, ChVector<>(0, 0, 0));
	//CreateArm(2, l_mod, ChVector<>( w / 2.0 +r2 + (l_mod / 2.0-r2)*cos(angle2), 0, -(l_mod / 2.0-r2)*sin(angle2)), quat2);

	////new version
	if (stapleSize)
	{
		l_mod = l + 2 * r2 - jointClearance;
		CreateArm(0, l_mod, ChVector<>(-w / 2.0 - (l / 2.0)*cos(angle1), 0, -(l_mod / 2.0 - r2)*sin(angle1)), quat0);
		CreateArm(1, w, ChVector<>(0, 0, 0));
		CreateArm(2, l_mod, ChVector<>(w / 2.0 +  (l / 2.0)*cos(angle2), 0, -(l_mod / 2.0 - r2)*sin(angle2)), quat2);
	}
	else
	{
		double armt = r;
		double armt2 = .0022 / 2 * sizeScale;
		l_mod = l + 2 * armt2 - jointClearance;
		CreateArm2(0, l_mod, armt, armt2, ChVector<>(-w / 2.0 - (l_mod / 2.0 - armt2)*cos(angle1), 0, -(l_mod / 2.0 - armt2)*sin(angle1)), quat0);
		CreateArm2(1, w,r,r2, ChVector<>(0, 0, 0));
		CreateArm2(2, l_mod, armt, armt2, ChVector<>( w / 2.0 + (l_mod / 2.0 - armt2)*cos(angle2), 0, -(l_mod / 2.0 - armt2)*sin(angle2)), quat2);
	}

	CreateActuators();
	//CreateJoints(); //TODO do we need joints?

	// mass property
	mass = arm0->GetMass() + arm1->GetMass() + arm2->GetMass();
	double r = ChRandom();
	//GetLog() << "rand number=" << r << "pctActive="<<pctActive <<"\n";
	if (r <= pctActive)
	{
		active = true;
		//controller(omega,force)
		armsController = (new Controller(m_system,this));
		armsController->outputLimit= torqueThresh2*2;
		armsController->omegaLimit = defaultOmega;

		//armsController->SetCurrAngle(0, this->GetCurrAngle(0));
		//armsController->SetCurrAngle(1, this->GetCurrAngle(1));
		//link_actuator01->GetLimit_Rz()->Set_min(-CH_C_PI);
		//link_actuator01->GetLimit_Rz()->Set_max(CH_C_PI);
		//link_actuator12->GetLimit_Rz()->Set_min(-CH_C_PI);
		//link_actuator12->GetLimit_Rz()->Set_max(CH_C_PI);
		//link_actuator01->GetLimit_Rz()->Set_active(true);
		//link_actuator12->GetLimit_Rz()->Set_active(true);
	}

}

ChSharedPtr<ChFunction> Smarticle::GetActuatorFunction(int actuatorID) {
	if (actuatorID == 0) {
		return function01;
	} else if (actuatorID == 1) {
		return function12;
	} else {
		std::cout << "Error! smarticle can only have actuators with ids from {0, 1}" << std::endl;
	}
	return ChSharedPtr<ChFunction>(NULL);
}

void Smarticle::SetActuatorFunction(int actuatorID, ChSharedPtr<ChFunction> actuatorFunction) {
	if (actuatorID == 0) {
		function01 = actuatorFunction;
		link_actuator01->Set_rot_funct(function01);
	} else if (actuatorID == 1) {
		function12 = actuatorFunction;
		link_actuator12->Set_rot_funct(function12);
	} else {
		std::cout << "Error! smarticle can only have actuators with ids from {0, 1}" << std::endl;
	}
}

void Smarticle::SetActuatorFunction(int actuatorID, double omega, double dT) {
	double diffTheta = dT * omega;
	ChSharedPtr<ChLinkEngine> mlink_actuator;
	if (actuatorID == 0) {
		mlink_actuator = link_actuator01;
	} else {
		mlink_actuator = link_actuator12;
	}
//	ChSharedPtr<ChFunction_Const> mfun1 = mlink_actuator->Get_rot_funct().DynamicCastTo<ChFunction_Const>();
//	mfun1->Set_yconst(diffTheta + mfun1->Get_yconst());
	ChSharedPtr<ChFunction_Const> mfun2 = mlink_actuator->Get_spe_funct().DynamicCastTo<ChFunction_Const>();
	ChSharedPtr<ChFunction_Const> mfun1 = mlink_actuator->Get_tor_funct().DynamicCastTo<ChFunction_Const>();

	
	mfun2->Set_yconst(omega);
	mfun1->Set_yconst(omega);
}

void Smarticle::SetActuatorFunction(int actuatorID, double omega) {
	ChSharedPtr<ChLinkEngine> mlink_actuator;
	if (actuatorID == 0) {
		mlink_actuator = link_actuator01;
	} else {
		mlink_actuator = link_actuator12;
	}
	ChSharedPtr<ChFunction_Const> mfun2 = mlink_actuator->Get_spe_funct().DynamicCastTo<ChFunction_Const>();
	ChSharedPtr<ChFunction_Const> mfun1 = mlink_actuator->Get_tor_funct().DynamicCastTo<ChFunction_Const>();
	mfun2->Set_yconst(omega);
	mfun1->Set_yconst(omega);
}

ChSharedPtr<ChLinkEngine> Smarticle::getLinkActuator(int id)
{
	if (id == 0)
		return link_actuator01;
	else
		return link_actuator12;
}
double Smarticle::GetVolume() {
//	return r * r2 * (w + 2 * (l + jointClearance));
	return (2 * r) * (2 * r2 )* (w + 2 * (l+2*r2));
}
double Smarticle::GetMass() {
	//	return r * r2 * (w + 2 * (l + jointClearance));
	return mass;
}

ChVector<> Smarticle::Get_cm() {
	return (arm0->GetMass() * arm0->GetPos() + arm1->GetMass() * arm1->GetPos() + arm2->GetMass() * arm2->GetPos()) / mass;

}

ChVector<> Smarticle::Get_InitPos() {
	return initPos;
}
void Smarticle::SetAngle(std::pair<double, double> mangles, bool degrees)
{
	SetAngles(mangles.first, mangles.second, degrees);
}
void Smarticle::SetInitialAngles()
{
	initialAng0 = this->GetAngle1();
	initialAng1 = this->GetAngle2();
}
void Smarticle::SetAngles(double mangle1, double mangle2, bool degrees)
{
	if (degrees)
	{
		angle1 = mangle1*CH_C_PI / 180.0;
		angle2 = mangle2*CH_C_PI / 180.0;
		return;
	}
	else
	{
		angle1 = mangle1;
		angle2 = mangle2;
	}
}
void Smarticle::SetAngle(int id, double mangle,bool degrees)
{
	if (id == 0)
		SetAngle1(mangle, degrees);
	else
		SetAngle2(mangle, degrees);
}
void Smarticle::SetAngle(double mangle, bool degrees)
{
	if (degrees)
	{
		angle1 = mangle*CH_C_PI / 180.0;
		angle2 = mangle*CH_C_PI / 180.0;
	}
	else
	{
		angle1 = mangle;
		angle2 = mangle;
	}
}
void Smarticle::SetAngle1(double mangle1, bool degrees)
{
	if (degrees) { angle1 = mangle1*CH_C_PI / 180.0; }
	else{ angle1 = mangle1; }
}
void Smarticle::SetAngle2(double mangle2, bool degrees)
{
	if (degrees) { angle2 = mangle2*CH_C_PI / 180.0; }
	else{ angle2 = mangle2; }
}
double Smarticle::GetAngle(int id, bool degrees)
{
	if (id == 0)
		return GetAngle1(degrees);
	else
		return GetAngle2(degrees);
}
double Smarticle::GetAngle1(bool degrees)
{
	if (degrees)
		return angle1*180.0 / CH_C_PI;
	else
		return angle1;
}
double Smarticle::GetAngle2(bool degrees)
{
	if (degrees)
		return angle2*180.0 / CH_C_PI;
	else
		return angle2;
}
bool Smarticle::NotAtDesiredPos(int id, double ang)//bad method name
{
	//GetLog() << "expAng" << id << ":" << GetExpAngle(id) << "\n    ";
	double x= ChooseOmegaAmount(GetOmega(id), ang, GetExpAngle(id)); //returns true if anything else but 0 is returned from here	
	//returns true if anything else but 0 is returned from here
	return x;
}
double Smarticle::GetExpAngle(int id)
{
	if (id == 0)
		return mv->at(moveTypeIdxs.at(moveType)).first;
	else
		return mv->at(moveTypeIdxs.at(moveType)).second;
}
double Smarticle::GetCurrAngle(int id)
{
	if (id ==0)
		return this->getLinkActuator(id)->Get_mot_rot();
	else 
		return this->getLinkActuator(id)->Get_mot_rot();
}
double Smarticle::GetNextAngle(int id)
{
	if (id==0)
		nextAngle.at(id) = mv->at((moveTypeIdxs.at(moveType) + 1) % mv->size()).first;
	else
		nextAngle.at(id) = mv->at((moveTypeIdxs.at(moveType) + 1) % mv->size()).second;
	
	return nextAngle.at(id);

	//if (id ==0)
	//	return mv->at((moveTypeIdxs.at(moveType) + 1) % mv->size()).first;
	//else
	//	return mv->at((moveTypeIdxs.at(moveType) + 1) % mv->size()).second;
}
void Smarticle::SetNextAngle(int id, double ang)
{
		nextAngle.at(id) = ang;
}
std::pair<double, double> Smarticle::populateMoveVector()
{

	std::ifstream smarticleMoves;
	smarticleMoves.open("smarticleMoves.csv");
	double mdt, momega, mtorqueThresh2, mangLow, mangHigh;
	smarticleMoves >>
		mdt >>
		momega >>
		mtorqueThresh2 >>
		mangLow >>
		mangHigh;
	//printf("dt %f omega %f torqueThresh2 %f angLow %f angHigh %f", mdt, momega, mtorqueThresh2, mangLow, mangHigh);
	SetDefaultOmega(momega);
	SetOmega(momega);
	char ddCh;
	char ddCh1;
	char ddCh2;
	angHigh = mangHigh;
	angLow = mangLow;
	distThresh = mdt*momega;
	torqueThresh2 = mtorqueThresh2;
	ddCh = '!';
	while (ddCh != '#') {
		smarticleMoves >> ddCh;
	}
	std::string ddSt;
	getline(smarticleMoves, ddSt);

	int smarticleCount = 0;
	std::pair<double, double> angPair;
	std::pair<double, double> firstAngPair;
	double ang1;
	double ang2;
	//smarticleMoves >> ang1 >> ddCh >> ang2 >> ddCh >> ddCh;
	//angPair.first = ang1;
	//angPair.second = ang2;

	//this->SetAngle(angPair.first, angPair.second);
	//this->SetAngle1(ang1);
	//this->SetAngle2(ang2);
	//for some reason it only works with vectors?
	ChVector<> angVals;
	smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
	angPair.first = angVals.x;
	angPair.second = angVals.y;
	firstAngPair.first = angVals.x;
	firstAngPair.second = angVals.y;
	//TODO need to rewrite the below way of reading file, very ugly!

	ot.push_back(angPair);
	ot.push_back(angPair);
	if (global.size() < 1)
	{
		ang1 = angPair.first;
		ang2 = angPair.second;
	global.push_back(angPair);

	//SetAngle(ang1, ang2);
	//SetAngles(ang1, ang2);
	//Global


		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
			angPair.first = angVals.x;
			angPair.second = angVals.y;

			global.push_back(angPair);
			//GetLog() << angVals.x << " " << angVals.y << " ddch:" << ddCh << "\n";
			if (ddCh == '#')
				break;
		}
	}

	if (gui1.size() < 1)
	{
		//GUI1
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
			angPair.first = angVals.x;
			angPair.second = angVals.y;

			gui1.push_back(angPair);
			//GetLog() << angVals.x << " " << angVals.y << " ddch:" << ddCh << "\n";
			if (ddCh == '#')
				break;
		}
	}

	if (gui2.size() < 1)
	{
		//GUI2
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
			angPair.first = angVals.x;
			angPair.second = angVals.y;

			gui2.push_back(angPair);
			//GetLog() << angVals.x << " " << angVals.y << " ddch:" << ddCh << "\n";
			//exit(-1);
			if (ddCh == '#')
				break;
		}
	}

	if (gui3.size() < 1)
	{
		//GUI3
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
			angPair.first = angVals.x;
			angPair.second = angVals.y;

			gui3.push_back(angPair);
			//GetLog() << angVals.x << " " << angVals.y << " ddch:" << ddCh << "\n";
			//exit(-1);
			if (ddCh == '#')
				break;
		}
	}


	//SetAngle(firstAngPair);
	//returning first ang pair but can be set here

		//ot.push_back(angPair);
		//ot.push_back(angPair);
	//exit(-1);
	return firstAngPair;

}
bool Smarticle::GetArm0OT()
{
	return this->arm0OT;
}
bool Smarticle::GetArm2OT()
{
	return this->arm2OT;
}
double Smarticle::ChooseOmegaAmount(double momega, double currAng, double destAng)
{
	//since going from -pi to pi:
	currAng = currAng + CH_C_PI;//TODO if starting value is not 0 degs, this can mess thigns up if pi is not added
	destAng = destAng + CH_C_PI;

	if (fabs(destAng - currAng) > distThresh)
	{
		//if destAng is larger, move with positive momega else
		if (destAng > currAng){return momega;}
		else{ return -1 * momega; }
	}
	//if <= distThresh dont move, let omega = 0;
	return 0;
}
bool Smarticle::MoveToAngle2(std::vector<std::pair<double, double>> *v, double momega1, double momega2,MoveType mtype)
{
	if (!active)
		return false;
	if (arm0OT || arm2OT) //turn off motion at limit
	{
		this->SetActuatorFunction(0, 0);
		this->SetActuatorFunction(1, 0);
		return false;
	}

	//real ang01 and ang12
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();

	SetAngles(ang01, ang12);
	
	//next ang01 and ang12
	double nextAng01 = v->at((moveTypeIdxs.at(moveType) + 1) % v->size()).first;
	double nextAng12 = v->at((moveTypeIdxs.at(moveType) + 1) % v->size()).second;


	double expAng01 = v->at(moveTypeIdxs.at(mtype)).first;
	double expAng12 = v->at(moveTypeIdxs.at(mtype)).second;
	
	
	//GetLog() << "(ang1,ang2)=(" << ang01 * 180 / CH_C_PI << "," << ang12 * 180 / CH_C_PI << ") (nextang1,nextang2)=(" << nextAng01 << ","<< nextAng12<<")\n";
	//exit(-1);
	//different real - expected
	//double ang01Diff = ang01-expAng01;
	//double ang12Diff = ang12-expAng12;

	double omega01 = ChooseOmegaAmount(momega1, ang01, expAng01);
	double omega12 = ChooseOmegaAmount(momega2, ang12, expAng12);

	
	//if within some threshold distance to where curr angle is supposed to be:
	if (omega01==0)
	{
		if (omega12 == 0)
		{
			//arm01=right, arm12=right both angles should try to move to the next angle
			this->SetActuatorFunction(0, ChooseOmegaAmount(momega1, ang01, nextAng01));
			this->SetActuatorFunction(1, ChooseOmegaAmount(momega2, ang12, nextAng12));
			return true;//was able to successfully move to next index
		}
		else
		{
			//arm01=right, arm12=wrong, arm12 must catch up
			this->SetActuatorFunction(0, 0); 
			this->SetActuatorFunction(1, omega12);
			return false;
		}
	}
	// from this point we know arm01 is wrong
	if (omega12 == 0)
	{
			//arm01=wrong, arm12=right
		this->SetActuatorFunction(0, omega01);
		this->SetActuatorFunction(1,0);
		return false;
	}
	else
	{
			//arm01=wrong, arm12=wrong
		this->SetActuatorFunction(0, omega01);
		this->SetActuatorFunction(1, omega12);
		return false;
	}
}
void Smarticle::MoveLoop2(int guiState = 0)
{
	//initialize boolean describing same moveType as last step
	bool sameMoveType = false;
	//static bool successfulMotion = false;
	bool prevSucessful = successfulMotion;
	successfulMotion = false;
	//this pointer will point to the correct moveType vector
	std::vector<std::pair<double, double>> *v;
	//set prevMoveType to previous value
	this->prevMoveType = this->moveType;
	//get previous values from last timestep
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();

	double omega1Prev = link_actuator01->Get_mot_rot_dt();
	double omega2Prev = link_actuator12->Get_mot_rot_dt();

	double torque01 = fabs(link_actuator01->Get_react_torque().Length2()); //use length2 to avoid squareroot calculation be aware of blowing up because too high torque overflows double
	double torque12 = fabs(link_actuator12->Get_react_torque().Length2());
	//GetLog() << "\n*********" << torque01 << " " << torque12 << " thresh: " << torqueThresh2;


	//determine moveType
	switch (guiState)
	{
	case 0:
		this->setCurrentMoveType(GLOBAL);
		v = &global;
		break;
	case 1:
		this->setCurrentMoveType(GUI1);
		v = &gui1;
		break;
	case 2:
		this->setCurrentMoveType(GUI2);
		v = &gui2;
		break;
	case 3:
		this->setCurrentMoveType(GUI3);
		v = &gui3;
		break;
	case 4:
		this->setCurrentMoveType(VIB);
		v = &vib;
		break;
	case 5:
		this->setCurrentMoveType(VIB);
		v = &vib;
		break;
	default:
		this->setCurrentMoveType(GLOBAL);
		v = &global;
		break;
	}




	static ChVector<> rel01 = link_actuator01->GetRelAxis();
	static ChVector<> rel12 = link_actuator12->GetRelAxis();
	static double relr01 = link_actuator01->GetDist();
	static double relr12 = link_actuator12->GetDist();
	//GetLog() << "\n" << "rel1:" << rel01 << "rel2:" << rel12 << "\n";


	//if (fabs(rel01.y) > .05 || fabs(rel12.y) > .05)//if angle in x or y is > .1 radians, it is definitely broken
	//{
	//	GetLog() << "angle bad break! \n";
	//	armBroken = true;
	//}
	//if (fabs(relr01) > .075*r2 || fabs(relr12) > .075*r2)//if distance between markers is .025% of thickness, break!
	//{
	//	GetLog() << "distance wrong break! \n";
	//	armBroken = true;
	//}

	//ChVector<> rel01 = link_actuator01->GetRelRotaxis();
	//ChVector<> rel12 = link_actuator12->GetRelRotaxis();
	////GetLog() <<"\n"<< rel01 << "\n";
	//if (fabs(rel01.x) > .1 || fabs(rel01.y) > .1
	//	|| fabs(rel12.x) > .1 || fabs(rel12.x) > .1)//if angle in x or y is > .1 radians, it is definitely broken
	//{
	//	armBroken = true;
	//}

	if (torque01 > torqueThresh2 || torque12 > torqueThresh2)//one arm is OT 
	{
		this->setCurrentMoveType(OT);
		v = &ot;

		if (torque01 > torqueThresh2) // arm 0 is overtorqued
		{
			if (!arm0OT)//if arm was previously not OT, add color asset
			{
				arm0OT = true;
				arm0->AddAsset(mtextureOT);
				this->ot.clear();
				this->ot.emplace_back(GetAngle1(), GetAngle2());
			}
		}
		if (torque12 > torqueThresh2)// arm 2 is overtorqued
		{
			if (!arm2OT)
			{
				arm2OT = true;
				arm2->AddAsset(mtextureOT);
				this->ot.clear();
				this->ot.emplace_back(GetAngle1(), GetAngle2());
			}
		}

	}
	else //arms are not OT 
	{
		if (arm0OT)//if arm was previously OT, pop off previous armOT color asset
		{
			arm0OT = false;
			arm0->GetAssets().pop_back();
		}
		if (arm2OT)
		{
			arm2OT = false;
			arm2->GetAssets().pop_back();
		}

	}

	sameMoveType = !(moveType^prevMoveType); // !(xor) gives true if values are equal, false if not
	
	switch (this->moveType) //have this in case I want to add different action based on move type
	{
		case GLOBAL://TODO implement different case if sameMoveType was wrong

			successfulMotion = MoveToAngle2(v, omega1, omega2,moveType);
			break;
		case OT:
			//if (sameMoveType){}
			successfulMotion = MoveToAngle2(v, 0, 0, moveType);
			break;
		case GUI1:

			if (sameMoveType)
			{
				if (omega1Prev == 0 && omega2Prev == 0)
				{

					successfulMotion = true;
					break;
				}
			}
			successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
			break;
		case GUI2:

			//if (sameMoveType){}
			successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
			break;
		case GUI3:

			//if (sameMoveType){}
			successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
			break;
		case VIB:
			successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
			//GetLog() << "(0,1,2):" << v->at(0).first << "," << v->at(1).first << "," << v->at(2).first;
			//exit(-1);
			break;
		}
	//add 1 to size if move was successful (i.e. can move on to next move index if reached previous one)
	if (successfulMotion&&active&&(!arm0OT&&!arm2OT))
	{
		moveTypeIdxs.at(moveType) = ((moveTypeIdxs.at(moveType) + 1) % v->size());
	}

	return;
}
void Smarticle::ChangeArmColor(double torque01, double torque12)
{
	double TT2 = torqueThresh2*1.99;
	double r0 = fabs(getLinkActuator(0)->Get_mot_rot_dt());
	double r1 = fabs(getLinkActuator(1)->Get_mot_rot_dt());
	double LIM = .1;
	double moveAmt = CH_C_PI / 90; //2 degrees
	//if (r0 || r1)
	//{
	//	GetLog() << "(r0,r1)= (" << r0 << "," <<  r1 << ")\n";
	//}
	if (fabs(torque01) > TT2&& r0<LIM)
	{
		this->setCurrentMoveType(OT);
		mv = &ot;
		if (!arm0OT)//if not previously OT
		{
			arm0OT = true;
			arm0->AddAsset(mtextureOT);
			this->ot.clear();
			this->ot.emplace_back(GetAngle1() + sign(torque01)*moveAmt, GetAngle2() + sign(torque01)*moveAmt);
			//this->ot.emplace_back(GetAngle1(), GetAngle2());
		}
		//nothing needs to be done if prev OT
	}
	else
	{
		if (arm0OT) //it prev OT but currently not
		{
			arm0OT = false;
			arm0->GetAssets().pop_back();
		}
		// nothing needs to be done if not prev OT
	}


	/////////////////////ARM2///////////////////////
	if (fabs(torque12) > TT2&& r1<LIM)
	{
		this->setCurrentMoveType(OT);
		mv = &ot;
		if (!arm2OT)//if not previously OT
		{
			arm2OT = true;
			arm2->AddAsset(mtextureOT);
			this->ot.clear();
			this->ot.emplace_back(GetAngle1() + sign(torque01)*moveAmt, GetAngle2() + sign(torque01)*moveAmt);
			//this->ot.emplace_back(GetAngle1(), GetAngle2());
		}
		//nothing needs to be done if prev OT
	}
	else
	{
		if (arm2OT) //it prev OT but currently not
		{
			arm2OT = false;
			arm2->GetAssets().pop_back();
		}
		// nothing needs to be done if not prev OT
	}
}
void Smarticle::AssignState(int guiState)
{
	switch (guiState)
	{
	case 0:
		this->setCurrentMoveType(GLOBAL);
		mv = &global;
		break;
	case 1:
		this->setCurrentMoveType(GUI1);
		mv = &gui1;
		break;
	case 2:
		this->setCurrentMoveType(GUI2);
		mv = &gui2;
		break;
	case 3:
		this->setCurrentMoveType(GUI3);
		mv = &gui3;
		break;
	case 4:
		this->setCurrentMoveType(VIB);
		mv = &vib;
		break;
	case 5:
		this->setCurrentMoveType(VIB);
		mv = &vib;
		break;
	default:
		this->setCurrentMoveType(GLOBAL);
		mv = &global;
		break;
	}

}
void Smarticle::ControllerMove(int guiState, double torque01, double torque12)
{
	bool sameMoveType = false;
	bool prevSucessful = successfulMotion;
	successfulMotion = false;
	this->prevMoveType = this->moveType;
	
	AssignState(guiState);
	ChangeArmColor(torque01, torque12);

	//!(moveType^prevMoveType)
	sameMoveType = (moveType==prevMoveType); // !(xor) gives true if values are equal, false if not
	if (!sameMoveType)
	{
		if (active)
			this->armsController->resetCumError = true;	
	}
	//successfulMotion = controller->step(sameMoveType,dT);
	
	//used to have switch but prob not necessary can just use if:

	if (active == false)
	{
		successfulMotion = false;
		return;
	}
	successfulMotion = armsController->Step(m_system->GetChTime());
	//if (this->moveType == OT); //if OT stop moving!
	//{
	//	
	//}

	if (successfulMotion)
		moveTypeIdxs.at(moveType) = ((moveTypeIdxs.at(moveType) + 1) % mv->size());
	//make sure controller step returns true false for movement!!

}
void Smarticle::MoveLoop2(int guiState, double torque01, double torque12)
{
	//initialize boolean describing same moveType as last step
	bool sameMoveType = false;
	//static bool successfulMotion = false;
	bool prevSucessful = successfulMotion;
	successfulMotion = false;
	//this pointer will point to the correct moveType vector
	std::vector<std::pair<double, double>> *v;
	//set prevMoveType to previous value
	this->prevMoveType = this->moveType;
	//get previous values from last timestep
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();

	double omega1Prev = link_actuator01->Get_mot_rot_dt();
	double omega2Prev = link_actuator12->Get_mot_rot_dt();

	

	//GetLog() << "\n*********" << torque01 << " " << torque12 << " thresh: " << torqueThresh2;


	//determine moveType
	switch (guiState)
	{
	case 0:
		this->setCurrentMoveType(GLOBAL);
		v = &global;
		break;
	case 1:
		this->setCurrentMoveType(GUI1);
		v = &gui1;
		break;
	case 2:
		this->setCurrentMoveType(GUI2);
		v = &gui2;
		break;
	case 3:
		this->setCurrentMoveType(GUI3);
		v = &gui3;
		break;
	case 4:
		this->setCurrentMoveType(VIB);
		v = &vib;
		break;
	case 5:
		this->setCurrentMoveType(VIB);
		v = &vib;
		break;
	default:
		this->setCurrentMoveType(GLOBAL);
		v = &global;
		break;
	}




	static ChVector<> rel01 = link_actuator01->GetRelAxis();
	static ChVector<> rel12 = link_actuator12->GetRelAxis();
	static double relr01 = link_actuator01->GetDist();
	static double relr12 = link_actuator12->GetDist();
	//GetLog() << "\n" << "rel1:" << rel01 << "rel2:" << rel12 << "\n";


	//if (fabs(rel01.y) > .05 || fabs(rel12.y) > .05)//if angle in x or y is > .1 radians, it is definitely broken
	//{
	//	GetLog() << "angle bad break! \n";
	//	armBroken = true;
	//}
	//if (fabs(relr01) > .075*r2 || fabs(relr12) > .075*r2)//if distance between markers is .025% of thickness, break!
	//{
	//	GetLog() << "distance wrong break! \n";
	//	armBroken = true;
	//}

	//ChVector<> rel01 = link_actuator01->GetRelRotaxis();
	//ChVector<> rel12 = link_actuator12->GetRelRotaxis();
	////GetLog() <<"\n"<< rel01 << "\n";
	//if (fabs(rel01.x) > .1 || fabs(rel01.y) > .1
	//	|| fabs(rel12.x) > .1 || fabs(rel12.x) > .1)//if angle in x or y is > .1 radians, it is definitely broken
	//{
	//	armBroken = true;
	//}

/////////////torque color change was here/////////////
	if (torque01 > torqueThresh2 || torque12 > torqueThresh2)//one arm is OT 
	{
		this->setCurrentMoveType(OT);
		v = &ot;
	}
	ChangeArmColor(torque01, torque12);

//////////////////////////////////////////////////////


	sameMoveType = !(moveType^prevMoveType); // !(xor) gives true if values are equal, false if not

	switch (this->moveType) //have this in case I want to add different action based on move type
	{
	case GLOBAL://TODO implement different case if sameMoveType was wrong

		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
		break;
	case OT:
		//if (sameMoveType){}
		successfulMotion = MoveToAngle2(v, 0, 0, moveType);
		break;
	case GUI1:

		//if (sameMoveType)
		//{
		//	if (omega1Prev == 0 && omega2Prev == 0)
		//	{

		//		successfulMotion = true;
		//		break;
		//	}
		//}
		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
		break;
	case GUI2:

		//if (sameMoveType){}
		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
		break;
	case GUI3:

		//if (sameMoveType){}
		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
		break;
	case VIB:
		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
		//GetLog() << "(0,1,2):" << v->at(0).first << "," << v->at(1).first << "," << v->at(2).first;
		//exit(-1);
		break;
	}
	//add 1 to size if move was successful (i.e. can move on to next move index if reached previous one)
	if (successfulMotion&&active && (!arm0OT&&!arm2OT))
	{
		moveTypeIdxs.at(moveType) = ((moveTypeIdxs.at(moveType) + 1) % v->size());
	}

	return;
}
ChSharedBodyPtr Smarticle::GetSmarticleBodyPointer()
{
	return arm1;
}
int	Smarticle::GetID()
{
	return arm0->GetPhysicsItem()->GetIdentifier();
}
void Smarticle::setCurrentMoveType(MoveType newMoveType)
{
	this->moveType = newMoveType;
	//global_GUI_value = newMoveType;
}

ChVector<> Smarticle::GetReactTorqueVectors01()
{
	return link_actuator01->Get_react_torque();
}
ChVector<> Smarticle::GetReactTorqueVectors12()
{
	return link_actuator12->Get_react_torque();
}
double Smarticle::GetReactTorqueLen01()
{
	//return (link_actuator01->Get_react_torque().Length2());
	return fabs((link_actuator01->Get_react_torque().z));
}
double Smarticle::GetZReactTorque(int id)
{
	return getLinkActuator(id)->Get_mot_torque();
}
double Smarticle::GetReactTorqueLen12()
{
	//return (link_actuator12->Get_react_torque().Length2());
	return fabs((link_actuator12->Get_react_torque().z));
}

void Smarticle::SetBodyFixed(bool mev){
	arm0->SetBodyFixed(mev);
	arm1->SetBodyFixed(mev);
	arm2->SetBodyFixed(mev);
}
