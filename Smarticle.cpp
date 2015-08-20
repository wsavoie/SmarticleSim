/*
 * Smarticle.cpp
 *
 *  Created on: May 22, 2015
 *      Author: Arman Pazouki
 */

#include "Smarticle.h"
#include "chrono_utils/ChUtilsGeometry.h"
#include "chrono_utils/ChUtilsCreators.h"


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
	jointClearance = 1.3 * r2;
	volume = GetVolume();
}

Smarticle::~Smarticle() {}


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

	jointClearance = 1.3 * r2;
	volume = GetVolume();

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
	moveTypeIdxs.resize(MoveType::OT, 0);
	double x = CH_C_PI/2;
	double y = x /10;
	//global.push_back(std::pair<double, double>(x + 0 * y, x + 0 * y));
	//global.push_back(std::pair<double, double>(x + 1 * y, x + 1 * y));
	//global.push_back(std::pair<double, double>(x + 2 * y, x + 2 * y));
	//global.push_back(std::pair<double, double>(x + 3 * y, x + 3 * y));
	//global.push_back(std::pair<double, double>(x + 2 * y, x + 2 * y));
	//global.push_back(std::pair<double, double>(x + 1 * y, x + 1 * y));
	//global.push_back(std::pair<double, double>(x - 0 * y, x - 0 * y));
	//global.push_back(std::pair<double, double>(x - 1 * y, x - 1 * y));
	//global.push_back(std::pair<double, double>(x - 2 * y, x - 2 * y));
	//global.push_back(std::pair<double, double>(x - 3 * y, x - 3 * y));
	//global.push_back(std::pair<double, double>(x - 2 * y, x - 2 * y));
	//global.push_back(std::pair<double, double>(x - 1 * y, x - 1 * y));
	//global.push_back(std::pair<double, double>(x + 0 * y, x + 0 * y));

	torqueThresh2 = 10000;
	angLow = 0; //TODO: should these be in radians?
	angHigh = 120;
	Smarticle::distThresh = .001*CH_C_PI; //TODO this threshold should be calculated from timestep and default omega rather than being a magic number maybe 1/2 momega*dt?
}
//////////////////////////////////////////////	
void Smarticle::SetDefaultOmega(double omega) {
	defaultOmega = omega;
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
void Smarticle::SetOmega(double momega1, double momega2, bool angularFreq)
{
	if (angularFreq)
	{
		omega1 = momega1;
		omega2 = momega2;
		return;
	}
	omega1 = (2 * CH_C_PI)*momega1;
	omega2 = (2 * CH_C_PI)*momega2;
}
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
	ChVector<> posArm = rotation.Rotate(posRel) + initPos;

	arm->SetName("smarticle_arm");
	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);
    arm->SetCollide(true);
    arm->SetBodyFixed(false);
    if (armID == 1)
    	arm->SetBodyFixed(false);
    else
    	arm->SetBodyFixed(false);

	arm->SetMaterialSurface(mat_g);

	double mass = density * vol;

	arm->GetCollisionModel()->ClearModel();

#if irrlichtVisualization
	ChSharedPtr<ChTexture> mtexture(new ChTexture());
	ChSharedPtr<ChBoxShape> box1(new ChBoxShape);
	box1->GetBoxGeometry().Size = ChVector<>(len / 2.0, r, r2);// upper part, max_x plate

	if (armID==1) 
		mtexture->SetTextureFilename(GetChronoDataFile("blu.png"));
	else 
		mtexture->SetTextureFilename(GetChronoDataFile("pink.png"));
	arm->GetCollisionModel()->AddBox(box1->GetBoxGeometry().Size.x, box1->GetBoxGeometry().Size.y, box1->GetBoxGeometry().Size.z, box1->Pos, box1->Rot);
	arm->GetAssets().push_back(box1);
	arm->AddAsset(mtexture);

#else
	utils::AddBoxGeometry(arm.get_ptr(), ChVector<>(len / 2.0, r, r2), ChVector<>(0, 0, 0));
	
#endif
	
	arm->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2
	arm->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelop);
    arm->GetCollisionModel()->BuildModel(); // this function overwrites the intertia

    // change mass and inertia property
    arm->SetMass(mass);
    arm->SetInertiaXX(mass * gyr);
    arm->SetDensity(density);

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
		std::cout << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
		break;
	}
}


ChSharedBodyPtr Smarticle::GetArm(int armID) {
	ChSharedBodyPtr arm;
	switch (armID) {
	case 0: {
		arm = arm0;
	} break;
	case 1: {
		arm = arm1;
	} break;
	case 2: {
		arm = arm2;
	} break;
	default:
		std::cout << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
		break;
	}
	return arm;
}

ChSharedPtr<ChLinkLockRevolute> Smarticle::GetRevoluteJoint(int jointID) {
	ChSharedPtr<ChLinkLockRevolute> link;
	switch (jointID) {
	case 0: {
		link = link_revolute01;
	} break;
	case 1: {
		link = link_revolute12;
	} break;
	default:
		std::cout << "Error! smarticle can only have joints with ids from {0, 1}" << std::endl;
		break;
	}
	return link;
}
void Smarticle::CreateJoints() {
	// link 1
	link_revolute01 = ChSharedPtr<ChLinkLockRevolute>(new ChLinkLockRevolute);
	ChVector<> pR01(-w/2, 0, 0);
	link_revolute01->Initialize(arm0, arm1,
        ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation * Q_from_AngAxis(CH_C_PI / 2, VECT_X)));
	m_system->AddLink(link_revolute01);


	// link 2
	link_revolute12 = ChSharedPtr<ChLinkLockRevolute>(new ChLinkLockRevolute);
	ChVector<> pR12(w/2, 0, 0);
	link_revolute12->Initialize(arm1, arm2,
        ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation * Q_from_AngAxis(CH_C_PI / 2, VECT_X)));
	m_system->AddLink(link_revolute12);
}

void Smarticle::CreateActuators() {
//	function01 = ChSharedPtr<ChFunction>(new ChFunction_Const(0));
//	function12 = ChSharedPtr<ChFunction>(new ChFunction_Const(0));

	// link 1
	link_actuator01 = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	ChVector<> pR01(-w/2, 0, 0);
	link_actuator01->Initialize(arm0, arm1,
		ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation * Q_from_AngAxis(CH_C_PI / 2.0, VECT_X)));
	link_actuator01->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
//	SetActuatorFunction(0, ChSharedPtr<ChFunction>(new ChFunction_Const(0)));
	link_actuator01->Get_mot_rot();
	
	m_system->AddLink(link_actuator01);

	

	// link 2
	link_actuator12 = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	ChVector<> pR12(w/2, 0, 0);
	link_actuator12->Initialize(arm1, arm2,
	        ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation * Q_from_AngAxis(CH_C_PI / 2, VECT_X)));
	link_actuator12->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
//	SetActuatorFunction(1, ChSharedPtr<ChFunction>(new ChFunction_Const(0)));
	m_system->AddLink(link_actuator12);
}

void Smarticle::Create() {
	jointClearance = 1.3*r2;
	double l_mod = l - jointClearance;

	// ** initialize U
//	CreateArm(0, l_mod, ChVector<>(-w/2 + r2, 0, l_mod/2 + r2 + jointClearance), Q_from_AngAxis(CH_C_PI / 2, VECT_Y));
//	CreateArm(1, w, ChVector<>(0, 0, 0));
//	CreateArm(2, l_mod, ChVector<>(w/2 - r2, 0, l_mod/2 + r2 + jointClearance), Q_from_AngAxis(CH_C_PI / 2, VECT_Y));

	
	ChQuaternion<> quat0 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, angle1, 0));
	ChQuaternion<> quat2 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, angle2, 0));	

	CreateArm(1, w, ChVector<>(0, 0, 0));
	CreateArm(0, l_mod, ChVector<>(-w / 2.0 - jointClearance- (r2 + l_mod / 2.0)*cos(angle1), 0, -(l_mod / 2.0 +jointClearance+ r2)*sin(angle1)),quat0);
	CreateArm(2, l_mod, ChVector<>(w / 2.0 + jointClearance + (r2 + l_mod / 2.0)*cos(angle2), 0, -(l_mod / 2.0 +jointClearance+ r2)*sin(angle2)),quat2);
	////////////////////////////////////////////////////////
	//CreateArm(1, w, ChVector<>(0, 0, 0));
	//CreateArm(0, l_mod, ChVector<>(-w / 2 - jointClearance - r2 - l_mod / 2, 0, 0));//original
	//CreateArm(2, l_mod, ChVector<>(w / 2 + jointClearance + r2 + l_mod / 2, 0, 0));//original
	/////////////////////////////////////////////////////////


	//-(w / 2.0 + jointClearance - r2) + (l_mod / 2.0 + r2)*cos(angle1)
	CreateJoints();
	CreateActuators();
	// mass property
	mass = arm0->GetMass() + arm1->GetMass() + arm2->GetMass();
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
	mfun2->Set_yconst(omega);
}

void Smarticle::SetActuatorFunction(int actuatorID, double omega) {
	ChSharedPtr<ChLinkEngine> mlink_actuator;
	if (actuatorID == 0) {
		mlink_actuator = link_actuator01;
	} else {
		mlink_actuator = link_actuator12;
	}
	ChSharedPtr<ChFunction_Const> mfun2 = mlink_actuator->Get_spe_funct().DynamicCastTo<ChFunction_Const>();
	mfun2->Set_yconst(omega);
}

double Smarticle::GetVolume() {
//	return r * r2 * (w + 2 * (l + jointClearance));
	return (2 * r) * (2 * r2 )* (w + 2 * l);
}

ChVector<> Smarticle::Get_cm() {
	return (arm0->GetMass() * arm0->GetPos() + arm1->GetMass() * arm1->GetPos() + arm2->GetMass() * arm2->GetPos()) / mass;

}

ChVector<> Smarticle::Get_InitPos() {
	return initPos;
}
void Smarticle::SetAngle(std::pair<double, double> mangles, bool degrees)
{
	SetAngle(mangles.first, mangles.second, degrees);
}
void Smarticle::SetAngle(double mangle1, double mangle2, bool degrees)
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

void Smarticle::AddMotion(ChSharedPtr<SmarticleMotionPiece> s_motionPiece) {
	int numAddedMotions = motion_vector.size();
	if (numAddedMotions == 0) {
		current_motion = s_motionPiece;
	}
	motion_vector.push_back(s_motionPiece);
}

void Smarticle::MoveLoop() {
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();
	
	double omega1 = current_motion->joint_01.omega;
	double omega2 = current_motion->joint_12.omega;

	if ((ang01 < current_motion->joint_01.theta1) && (omega1 < 0)) {
		omega1 *= -1;
	}
	if ((ang01 > current_motion->joint_01.theta2) && (omega1 > 0)) {
		omega1 *= -1;
	}
	if ((ang12 < current_motion->joint_12.theta1) && (omega2 < 0)) {
		omega2 *= -1;
	}
	if ((ang12 > current_motion->joint_12.theta2) && (omega2 > 0)) {
		omega2 *= -1;
	}

	current_motion->joint_01.omega = omega1;
	current_motion->joint_12.omega = omega2;

	this->SetActuatorFunction(0, omega1);
	this->SetActuatorFunction(1, omega2);

}
std::pair<double, double> Smarticle::populateMoveVector(std::vector<std::pair<double, double>> &mglobal, std::vector<std::pair<double, double>> &mOT, std::vector<std::pair<double, double>> &mGUI1, std::vector<std::pair<double, double>> &mGUI2)
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
		printf("dt %f omega %f torqueThresh2 %f angLow %f angHigh %f",mdt, momega, mtorqueThresh2, mangLow, mangHigh);
		this->SetDefaultOmega(momega);
		this->SetOmega(momega);
		//this->SetOmega(momega);
		char ddCh;
		char ddCh1;
		char ddCh2;
		this->angHigh = mangHigh;
		this->angLow = mangLow;
		this->distThresh = mdt*omega1;
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
	GetLog() << "\n";
	//for some reason it only works with vectors?
	ChVector<> angVals;
	smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
	angPair.first = angVals.x;
	angPair.second = angVals.y;
	firstAngPair.first = angVals.x;
	firstAngPair.second = angVals.y;
	mglobal.push_back(angPair);
	//GetLog() << angVals.x << " " << angVals.y << " ddch:" << ddCh << "\n";
	this->SetAngle(ang1, ang2);
	//Global
	while (smarticleMoves.good()) {
		smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
		angPair.first = angVals.x;
		angPair.second = angVals.y;
		
		mglobal.push_back(angPair);
		//GetLog() << angVals.x << " " << angVals.y << " ddch:" << ddCh << "\n";
		if (ddCh == '#')
			break;
	}
	
	//GUI1
	while (smarticleMoves.good()) {
		smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
		angPair.first = angVals.x;
		angPair.second = angVals.y;
		
		mGUI1.push_back(angPair);
		//GetLog() << angVals.x << " " << angVals.y << " ddch:" << ddCh << "\n";
		if (ddCh == '#')
			break;
	}
	//GUI2
	while (smarticleMoves.good()) {
		smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
		angPair.first = angVals.x;
		angPair.second = angVals.y;

		mGUI2.push_back(angPair);
		//GetLog() << angVals.x << " " << angVals.y << " ddch:" << ddCh << "\n";
		//exit(-1);
		if (ddCh == '#')
			break;
	}

	SetAngle(firstAngPair);
	//returning first ang pair but can be set here


	//exit(-1);
	return firstAngPair;

}

double Smarticle::ChooseOmegaAmount(double momega, double currAng, double destAng)
{//TODO make sure to determine if omega is greater than angHigh or angLow too!!!!
	//since going from -pi to pi:
	currAng = CH_C_PI + currAng;
	destAng = CH_C_PI + destAng;
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
	//real ang01 and ang12
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();

	
	//expected ang01 and ang12
	double expAng01 = v->at(moveTypeIdxs.at(mtype)).first;
	double expAng12 = v->at(moveTypeIdxs.at(mtype)).second;

	//next ang01 and ang12
	double nextAng01 = v->at((moveTypeIdxs.at(moveType) + 1) % v->size()).first;
	double nextAng12 = v->at((moveTypeIdxs.at(moveType) + 1) % v->size()).second;

	//GetLog() << "(ang1,ang2)=(" << ang01 * 180 / CH_C_PI << "," << ang12 * 180 / CH_C_PI << ") (nextang1,nextang2)=(" << nextAng01 << ","<< nextAng12<<")\n";
	//exit(-1);
	//different real - expected
	double ang01Diff = ang01-expAng01;
	double ang12Diff = ang12-expAng12;
	
	double omega01 = ChooseOmegaAmount(momega1, ang01, expAng01);
	double omega12 = ChooseOmegaAmount(momega2, ang12, expAng12);
	//if within some threshold distance to where curr angle is supposed to be:
	if (omega01==0)
	{
		if (omega12 == 0)
		{
			//arm01=right, arm12=right
			this->SetActuatorFunction(0, ChooseOmegaAmount(momega1, ang01, nextAng01));
			this->SetActuatorFunction(1, ChooseOmegaAmount(momega2, ang12, nextAng12));
			return true;//was able to successfully move to next index
		}
		else
		{
			//arm01=right, arm12=wrong
			this->SetActuatorFunction(0, 0);
			this->SetActuatorFunction(1, omega12);
			return false;
		}	
	}
	// from this point we know arm01 is wrong
	if (ChooseOmegaAmount(momega2, ang12, expAng12) == 0)
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
	bool successfulMotion = false;
	//this pointer will point to the correct moveType vector
	std::vector<std::pair<double, double>> *v;
	
	//set prevMoveType to previous value
	this->prevMoveType = this->moveType;

	//get previous values from last timestep
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();
	//double omega01 = li
	double torque01 = link_actuator01->Get_react_torque().Length2(); //use length2 to avoid squareroot calculation be aware of blowing up because too high torque overflows double
	double torque12 = link_actuator12->Get_react_torque().Length2();



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
	}
	//overTorque takes priority!
	if (torque01 > torqueThresh2 || torque12 > torqueThresh2){
		//this->setCurrentMoveType(OT);
		//v = &ot;
		GetLog() << "*************************\n";
		GetLog() << "Torque overload!\n";
		GetLog() << "*************************\n";
	}

	if (this->moveType == this->prevMoveType){
		sameMoveType = true;
	}
	static bool x = false;

	switch (this->moveType) //TODO finish move type case statements
	{
		case GLOBAL://TODO implement different case if sameMoveType was wrong
			//TODO finish global move case

			successfulMotion = MoveToAngle2(v, omega1, omega2,moveType);
			break;
		case OT:
			//TODO finish ot move case
			if (sameMoveType)
			{
			}
			break;
		case GUI1:
			//TODO finish gui1 move case
			if (sameMoveType)
			{
				x = true;
			}
			successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);

			break;
		case GUI2:
			//TODO finish gui2 move case
			if (sameMoveType)
			{
			}
			successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
			break;
		case GUI3:
			//TODO finish gui3 move case
			if (sameMoveType)
			{
			}
			break;
	}
	//add 1 to size if
	if (successfulMotion)
	{
		moveTypeIdxs.at(moveType) = ((moveTypeIdxs.at(moveType) + 1) % v->size());
	}
	
	return;
}
void Smarticle::setCurrentMoveType(MoveType newMoveType)
{
	this->moveType = newMoveType;
}
bool Smarticle::MoveToRange() {
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();
	double omega1 = current_motion->joint_01.omega;
	double omega2 = current_motion->joint_12.omega;

	int low1 = ang01 < current_motion->joint_01.theta1;
	int high1 = ang01 > current_motion->joint_01.theta2;
	int low2 = ang12 < current_motion->joint_12.theta1;
	int high2 = ang12 > current_motion->joint_12.theta2;

	bool isInRange = (!low1) && (!high1) && (!low2) && (!high2);

	if (low1) {
		this->SetActuatorFunction(0, defaultOmega);
	} else if (high1) {
		this->SetActuatorFunction(0, -defaultOmega);
	}

	if (low2) {
		this->SetActuatorFunction(1, defaultOmega);
	} else if (high2) {
		this->SetActuatorFunction(1, -defaultOmega);
	}

	return isInRange;
}

void Smarticle::MoveSquare() {
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();

	int low1 = ang01 < current_motion->joint_01.theta1;
	int high1 = ang01 > current_motion->joint_01.theta2;
	int low2 = ang12 < current_motion->joint_12.theta1;
	int high2 = ang12 > current_motion->joint_12.theta2;

	if (low2 && !low1) {
		current_motion->motionSubSegment = 2;
	} else if (low1 && !high2) {
		current_motion->motionSubSegment = 3;
	} else if (high2 && !high1) {
		current_motion->motionSubSegment = 0;
	} else if (high1 && !low2) {
		current_motion->motionSubSegment = 1;
	} else {
		current_motion->motionSubSegment = 0;
	}

	switch (current_motion->motionSubSegment) {
	case 0:
		current_motion->joint_01.omega = defaultOmega;
		current_motion->joint_12.omega = 0;
		break;
	case 1:
		current_motion->joint_01.omega = 0;
		current_motion->joint_12.omega = -defaultOmega;
		break;
	case 2:
		current_motion->joint_01.omega = -defaultOmega;
		current_motion->joint_12.omega = 0;
		break;
	case 3:
		current_motion->joint_01.omega = 0;
		current_motion->joint_12.omega = defaultOmega;
		break;
	}

	this->SetActuatorFunction(0, current_motion->joint_01.omega);
	this->SetActuatorFunction(1, current_motion->joint_12.omega);
}

void Smarticle::MoveCircle() {
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();

	double th1_low = current_motion->joint_01.theta1;
	double th1_hig = current_motion->joint_01.theta2;
	double th2_low = current_motion->joint_12.theta1;
	double th2_hig = current_motion->joint_12.theta2;

	double r1 = 0.5 * (th1_hig - th1_low);
	double r2 = 0.5 * (th2_hig - th2_low);

	if (ang01 > r1) {
		ang01 = r1;
	} else if (ang01 < -r1) {
		ang01 = -r1;
	}

//	double gamma = asin(ang01 / r1);
	double gamma = atan2(ang01, ang12);

//	while (fabs(gamma) < .01 * CH_C_PI) {
//		gamma += .01 * CH_C_PI;
//	}
	double gammaDot = defaultOmega / r1;

	// make sure the initial configuration is on the sphere (up to a tolerance)

	double omega1 = r1 * cos(gamma) * gammaDot;
//	double omega2 = -ang01 * omega1 / (ang12);
	double omega2 = -r2 * sin(gamma) * gammaDot;

	double smallOmega = .001 * defaultOmega;
	if (fabs(omega1) < smallOmega) {
		if (ang01 < 0) {
			omega1 = smallOmega;
		} else {
			omega1 = -smallOmega;
		}
	}
	if (fabs(omega2) < smallOmega) {
		if (ang12 < 0) {
			omega2 = smallOmega;
		} else {
			omega2 = -smallOmega;
		}
	}

	printf(" omega1 %f omega2 %f\n", omega1, omega2);

	this->SetActuatorFunction(0, omega1);
	this->SetActuatorFunction(1, omega2);
}

void Smarticle::MoveToAngle(double theta1, double theta2) {
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();

	double omega1;
	double omega2;

	if (fabs(ang01 - theta1) < .01 * CH_C_PI) {
		omega1 = 0;
	} else if (ang01 < theta1) {
		omega1 = defaultOmega;
	} else {
		omega1 = -defaultOmega;
	}

	if (fabs(ang12 - theta2) < .01 * CH_C_PI) {
		omega2 = 0;
	} else if (ang12 < theta2) {
		omega2 = defaultOmega;
	} else {
		omega2 = -defaultOmega;
	}

	this->SetActuatorFunction(0, omega1);
	this->SetActuatorFunction(1, omega2);
}

void Smarticle::MoveRelease() {
	this->MoveToAngle(0, 0);
}

void Smarticle::UpdateMySmarticleMotion() {
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();

//	printf("theta1 min %f max %f and theta 2 min %f max %f \n", current_motion->joint_01.theta1, current_motion->joint_01.theta2, current_motion->joint_12.theta1, current_motion->joint_12.theta2);
//
//	bool inRange = true;
//
//	if (ang01 < current_motion->joint_01.theta1) {
//		this->SetActuatorFunction(0, defaultOmega);
////		current_motion->joint_01.omega = defaultOmega;
//		printf("1 ang01 %f joint_01.theta1 %f omega %f\n", ang01, current_motion->joint_01.theta1, defaultOmega);
//		inRange = false;
//	}
//	if (ang01 > current_motion->joint_01.theta2) {
//		this->SetActuatorFunction(0, -defaultOmega);
////		current_motion->joint_01.omega = -defaultOmega;
//		printf("2 ang01 %f joint_01.theta2 %f omega %f\n", ang01, current_motion->joint_01.theta2, -defaultOmega);
//		inRange = false;
//	}
//	if (ang12 < current_motion->joint_12.theta1) {
//		this->SetActuatorFunction(1, defaultOmega);
////		current_motion->joint_12.omega = defaultOmega;
//		printf("3 ang12 %f joint_12.theta1 %f omega %f\n", ang12, current_motion->joint_12.theta1, defaultOmega);
//		inRange = false;
//	}
//	if (ang12 > current_motion->joint_12.theta2) {
//		this->SetActuatorFunction(1, -defaultOmega);
////		current_motion->joint_12.omega = -defaultOmega;
//		printf("4 ang12 %f joint_12.theta1 %f omega %f\n", ang12, current_motion->joint_12.theta2, -defaultOmega);
//
//		inRange = false;
//	}
//
//	if (!inRange) {
//		printf("yoyoyuoyoyuoyoy\n");
//		return;
//	}
//	printf("bw\n");

	static int count_segment = 0;
	current_motion = motion_vector[1];
	if (motion_vector.size() < 2) {
		return;
	} else {
		count_segment = 1;
		current_motion = motion_vector[count_segment]; // 0 is reserved for release
	}
	if (count_segment < motion_vector.size() - 1) {
		if (motion_vector[count_segment + 1]->startTime < m_system->GetChTime()) {
			count_segment ++;
			current_motion = motion_vector[count_segment];
		}
	}



	MotionType sMotion = current_motion->GetMotionType();
	switch (current_motion->GetMotionType()) {
	case SQUARE_G:
		MoveSquare();
		break;
	case CIRCLE_G:
		MoveCircle();
		break;
	case RELEASE_G:
		MoveRelease();
		break;
	case LOOP_G:
		MoveLoop();
		break;
	default:
		break;
	}
}

ChSharedPtr<SmarticleMotionPiece> Smarticle::Get_Current_Motion() {
	return current_motion;
}

void Smarticle::SetBodyFixed(bool mev){
	arm0->SetBodyFixed(mev);
	arm1->SetBodyFixed(mev);
	arm2->SetBodyFixed(mev);
}