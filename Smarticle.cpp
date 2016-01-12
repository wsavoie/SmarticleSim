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
//#include "chrono_parallel/physics/ChSystemParallel.h"
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
std::vector<std::pair<double, double>> Smarticle::midTorque;

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
	mtextureMid->SetTextureFilename(GetChronoDataFile("cubetexture_blue_bordersBlueOriented.png"));
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
	omega1 = (2 * PI)*momega;
	omega2 = (2 * PI)*momega;
}
//void Smarticle::SetOmega(double momega1, double momega2, bool angularFreq)
//{
//	if (angularFreq)
//	{
//		omega1 = momega1;
//		omega2 = momega2;
//		return;
//	}
//	omega1 = (2 * PI)*momega1;
//	omega2 = (2 * PI)*momega2;
//}
void Smarticle::SetOmega1(double momega1, bool angularFreq)
{
	if (angularFreq)
	{
		omega1 = momega1;
		return;
	}
	omega1 = (2 * PI)*momega1;
}
void Smarticle::SetOmega2(double momega2, bool angularFreq)
{
	if (angularFreq)
	{
		omega2=momega2;
		return;
	}
	omega2 = (2 * PI)*momega2;
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
	return omega1 / (2 * PI);
}

double Smarticle::GetOmega2(bool angularFreq)
{
	if (angularFreq)
	{
		return omega2;
	}
	return omega2 / (2 * PI);
}

void Smarticle::CreateArm(int armID, double len, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChVector<> gyr;  	// components gyration
	double vol;			// components volume

	ChSharedBodyPtr arm;

	vol = utils::CalcBoxVolume(ChVector<>(len/2.0, r, r2));
	gyr = utils::CalcBoxGyration(ChVector<>(len/2.0, r, r2)).Get_Diag();
	// create body, set position and rotation, add surface property, and clear/make collision model
	arm = ChSharedBodyPtr(new ChBody);


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
			arm->SetBodyFixed(false);
    else
    	arm->SetBodyFixed(false);

	arm->SetMaterialSurface(mat_g);

	double mass = density * vol;
	//double mass = .005;//.043/3.0; //robot weight 43 grams
	arm->GetCollisionModel()->ClearModel();
	//arm->SetLimitSpeed(true);
	//arm->SetMaxSpeed(PI * 2 * 5);
	//arm->SetMaxWvel(PI * 2 * 5);
	//arm->ClampSpeed();
	
	if (visualize)
	{
		arm0_texture = ChSharedPtr<ChTexture>(new ChTexture());
		arm1_texture = ChSharedPtr<ChTexture>(new ChTexture());
		arm2_texture = ChSharedPtr<ChTexture>(new ChTexture());
		switch (armID) {
		case 0:
			arm0_texture = ChSharedPtr<ChTexture>(new ChTexture());
			arm0_texture = mtextureArm;
			arm->AddAsset(arm0_texture);
			break;
		case 1:
			arm1_texture = ChSharedPtr<ChTexture>(new ChTexture());
			arm1_texture = mtextureMid;
			arm->AddAsset(arm1_texture);
			break;
		case 2:
			arm2_texture = ChSharedPtr<ChTexture>(new ChTexture());
			arm2_texture = mtextureArm;
			arm->AddAsset(arm2_texture);
			break;
		default:
			std::cerr << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
			break;
		}
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
	arm = ChSharedBodyPtr(new ChBody);

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
		arm0_texture = ChSharedPtr<ChTexture>(new ChTexture());
		arm1_texture = ChSharedPtr<ChTexture>(new ChTexture());
		arm2_texture = ChSharedPtr<ChTexture>(new ChTexture());
		switch (armID) {
		case 0:
			arm0_texture = ChSharedPtr<ChTexture>(new ChTexture());
			arm0_texture = mtextureArm;
			arm->AddAsset(arm0_texture);
			break;
		case 1:
			arm1_texture = ChSharedPtr<ChTexture>(new ChTexture());
			arm1_texture = mtextureMid;
			arm->AddAsset(arm1_texture);
			break;
		case 2:
			arm2_texture = ChSharedPtr<ChTexture>(new ChTexture());
			arm2_texture = mtextureArm;
			arm->AddAsset(arm2_texture);
			break;
		default:
			std::cerr << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
			break;
		}
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
	ChQuaternion<> qx = Q_from_AngAxis(PI_2, VECT_X);

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


	ChQuaternion<> qx = Q_from_AngAxis(PI_2, VECT_X);
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

	ChQuaternion<> quat0 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0 * 5 * D2R, angle1, 0));
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
		double armt2 = .00806 / 2 * sizeScale; //8.06 mm
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
		armsController->outputLimit= torqueThresh2;
		armsController->omegaLimit = omegaLim;

		//armsController->SetCurrAngle(0, this->GetCurrAngle(0));
		//armsController->SetCurrAngle(1, this->GetCurrAngle(1));
		//link_actuator01->GetLimit_Rz()->Set_min(-PI);
		//link_actuator01->GetLimit_Rz()->Set_max(PI);
		//link_actuator12->GetLimit_Rz()->Set_min(-PI);
		//link_actuator12->GetLimit_Rz()->Set_max(PI);
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
		angle1 = mangle1*D2R;
		angle2 = mangle2*D2R;
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
		angle1 = mangle*D2R;
		angle2 = mangle*D2R;
	}
	else
	{
		angle1 = mangle;
		angle2 = mangle;
	}
}
void Smarticle::SetAngle1(double mangle1, bool degrees)
{
	if (degrees) { angle1 = mangle1*D2R; }
	else{ angle1 = mangle1; }
}
void Smarticle::SetAngle2(double mangle2, bool degrees)
{
	if (degrees) { angle2 = mangle2*D2R; }
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
		return angle1*R2D;
	else
		return angle1;
}
double Smarticle::GetAngle2(bool degrees)
{
	if (degrees)
		return angle2*R2D;
	else
		return angle2;
}
void Smarticle::addInterpolatedPathToVector(double a0i, double a2i, double a0f, double a2f)
{
	double dist1 = omega1*dT;
	double dist2 = omega2*dT;
	int n = std::max(abs((a0f - a0i) / dist1), abs((a2f - a2i) / dist2));
	std::vector<double> a0 = linspace(a0i, a0f, n);
	std::vector<double> a2 = linspace(a2i, a2f, n);

	for (int i = 0; i < n; i++)
	{
		mv->emplace_back(a0.at(i), a2.at(i));
	}

}
std::vector<double> Smarticle::linspace(double a, double b, int n) {
	std::vector<double> vec;
	double step = (b - a) / (n - 1);

	for(int i = 0; i < n; ++i)
	{
		vec.push_back(a);
		a += step;           // could recode to better handle rounding errors
	}
	return vec;
}
bool Smarticle::NotAtDesiredPos(int id, double ang,double exp)//bad method name
{
	//GetLog() << "expAng" << id << ":" << GetExpAngle(id) << "\n    ";
	double x= ChooseOmegaAmount(GetOmega(id), ang, exp); //returns true if anything else but 0 is returned from here	
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
	double mdt, momega, mtorqueThresh2, mangLow, mangHigh,mOmegaLim;
	smarticleMoves >>
		mdt >>
		momega >>
		mtorqueThresh2 >>
		mangLow >>
		mangHigh>>
		mOmegaLim;
	//printf("dt %f omega %f torqueThresh2 %f angLow %f angHigh %f", mdt, momega, mtorqueThresh2, mangLow, mangHigh);
	SetDefaultOmega(momega);
	SetOmega(momega);
	omegaLim = mOmegaLim;

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

	if (midTorque.size() < 1)
	{
		//midtorque
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
			angPair.first = angVals.x;
			angPair.second = angVals.y;

			midTorque.push_back(angPair);
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
	//currAng = currAng + PI;
	//destAng = destAng + PI;
	double deltaAng = destAng - currAng;
	if (abs(deltaAng) > 2*distThresh)
	{		//if destAng is larger, move with positive momega
		return sgn(deltaAng)*momega;
	}
	//if <= distThresh dont move, let omega = 0;
	return 0;
}

void Smarticle::ChangeStateBasedOnTorque(double tor0, double tor1,double timeSinceChange)
{
	//torque01 and torque02 are averaged torque over some amount of steps
	//low thresh		= if both torques are < LT*thresh.
	//med thresh		= if both torques are LT<x<HT
	//hi  thresh		= if one torque is > HT
	
	
	//if (timeSinceChange != 0 && timeSinceChange < 10 * dT) //protects situation where torque increases initially causing LT situations to get out of LT immediately upon shape change
	//	return;


	double LT = .85 * torqueThresh2;
	double MT = .95 * torqueThresh2;
	//double HT = 1.99 *torqueThresh2;

	//double t0 = abs(tor0);
	//double t1 = abs(tor1);
	double t0 = abs(tor0);
	double t1 = abs(tor1);
	if (GetArm0OT() || GetArm2OT())
	{//highest torque threshold, stop moving
		//AssignState(OT);
		specialState = -1;
		return;
	}
	else
	{//
		if (t0 > MT && t1 > MT) //if greater than MT, (and less than OT because above if)
		{
			specialState = MIDT;
			//AssignState(MIDT);
			//TODO clear midt and emplace values
			return;
			
		}
		else if (t0 < LT && t1 < LT) //todo abs value of torque
		//if (t0 < LT && t1 < LT)
		{//LOW TORQUE
			if (lowStressChange)		//if time to switch states
			{
				if (specialState == GUI1 || global_GUI_value==GUI1) //was not already in special state
				{
					specialState = GUI2;
				}
				else//if already in special state switch to a different one 
				{
					specialState = GUI1;
				}
				//AssignState(specialState);
				return;
				//ss.clear();
				//addInterpolatedPathToVector()

			}
			//does specialState = -1; go here?
			return; //maybe to a low torque color change?
		}
	}
	specialState = -1;
}
void Smarticle::ChangeArmColor(double torque01, double torque12)
{
	double TT2 = torqueThresh2*.99;
	double r0 = abs(getLinkActuator(0)->Get_mot_rot_dt());
	double r1 = abs(getLinkActuator(1)->Get_mot_rot_dt());
	double LIM = .1;
	double moveAmt = 2*D2R; //2 degrees
	if (abs(torque01) > TT2)
	{
		arm0_texture = mtextureOT;
		//this->setCurrentMoveType(OT);
		//mv = &ot;
		if (!arm0OT)//if not previously OT
		{
			this->ot.clear();
			//this->ot.emplace_back(GetAngle1() + sign(torque01)*moveAmt, GetAngle2() + sign(torque12)*moveAmt);
			this->ot.emplace_back(GetAngle1(), GetAngle2());
			this->ot.emplace_back(GetAngle1(), GetAngle2());
		}
		arm0OT = true;
		//nothing needs to be done if prev OT
	}
	else
	{
		arm0_texture = mtextureArm;
		arm0OT = false;
		// nothing needs to be done if not prev OT
	}


	/////////////////////ARM2///////////////////////
	if (abs(torque12) > TT2)
	{
		arm2_texture = mtextureOT;
		//this->setCurrentMoveType(OT);
		//mv = &ot;
		if (!arm2OT)//if not previously OT
		{
			this->ot.clear();
			//this->ot.emplace_back(GetAngle1() + sign(torque01)*moveAmt, GetAngle2() + sign(torque12)*moveAmt);
			this->ot.emplace_back(GetAngle1(), GetAngle2());
			this->ot.emplace_back(GetAngle1(), GetAngle2());
		}
		arm2OT = true;
		//nothing needs to be done if prev OT
	}
	else
	{
		arm2_texture = mtextureArm;
		arm2OT = false;
		// nothing needs to be done if not prev OT
	}
}
void Smarticle::GenerateVib(double ang1, double ang2)
{
	this->addInterpolatedPathToVector(ang1, ang2, ang1 + vibAmp, ang2 + vibAmp);//curr				->		curr+vib
	this->addInterpolatedPathToVector(ang1 + vibAmp, ang2 + vibAmp, ang1, ang2);//curr+vib		->		curr
	this->addInterpolatedPathToVector(ang1, ang2, ang1 - vibAmp, ang2 - vibAmp);//curr				->		curr-vib
	this->addInterpolatedPathToVector(ang1 - vibAmp, ang2 - vibAmp, ang1, ang2);//curr-vib		->		curr
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
		this->setCurrentMoveType(MIDT);
		mv = &midTorque;
		break;
	case 6:
		this->setCurrentMoveType(SS);
		mv = &ss;
		break;
	case 7:
		this->setCurrentMoveType(OT);
		mv = &ot;
		break;
	default:
		this->setCurrentMoveType(GLOBAL);
		mv = &global;
		break;
	}

}
double Smarticle::CheckLowStressChangeTime()
{
	static double timeSinceLastChange = 0;
	timeSinceLastChange += dT;
	lowStressChange = false;//always reset
	if (timeSinceLastChange > gaitLengthChangeTime)
	{
		lowStressChange = true;
		timeSinceLastChange = 0; //reset value
		return 0;
	}
	return timeSinceLastChange;
}
void Smarticle::ControllerMove(int guiState, double torque01, double torque12)
{
	
	if (active == false)//TODO put this higher
	{
		successfulMotion = false;
		return;
	}

	bool sameMoveType = false;
	bool prevSucessful = successfulMotion;
	successfulMotion = false;
	this->prevMoveType = this->moveType;
	double timeSinceChange = CheckLowStressChangeTime();
	ChangeArmColor(torque01, torque12);
	ChangeStateBasedOnTorque(torque01,torque12,timeSinceChange);
	if (specialState != -1)
		AssignState(specialState);
	else
		AssignState(guiState);

	//!(moveType^prevMoveType)
	sameMoveType = (moveType==prevMoveType); // !(xor) gives true if values are equal, false if not
	if (!sameMoveType)
	{
			this->armsController->resetCumError = true;	
	}
	//successfulMotion = controller->step(sameMoveType,dT);
	
	//used to have switch but prob not necessary can just use if:

	
	successfulMotion = armsController->Step(m_system->GetChTime());
	//if (this->moveType == OT); //if OT stop moving!
	//{
	//	
	//}

	if (successfulMotion)
		moveTypeIdxs.at(moveType) = ((moveTypeIdxs.at(moveType) + 1) % mv->size());
	//make sure controller step returns true false for movement!!

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
