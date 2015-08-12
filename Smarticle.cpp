/*
 * Smarticle.cpp
 *
 *  Created on: May 22, 2015
 *      Author: Arman Pazouki
 */

#include "Smarticle.h"
#include "chrono_utils/ChUtilsGeometry.h"
#include "chrono_utils/ChUtilsCreators.h"


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
	jointClearance = .05 * r2;
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

	jointClearance = .05 * r2;
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

}

void Smarticle::SetDefaultOmega(double omega) {
	defaultOmega = omega;
}


void Smarticle::CreateArm(int armID, double len, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChVector<> gyr;  	// components gyration
	double vol;			// components volume
	double jointClearance = 1.3 * r2;

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
//    if (armID == 1)
//    	arm->SetBodyFixed(true);
//    else
//    	arm->SetBodyFixed(false);

	arm->SetMaterialSurface(mat_g);

	double mass = density * vol;

	arm->GetCollisionModel()->ClearModel();
	utils::AddBoxGeometry(arm.get_ptr(), ChVector<>(len/2.0, r, r2), ChVector<>(0, 0, 0));
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
	        ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation * Q_from_AngAxis(CH_C_PI / 2, VECT_X)));
	link_actuator01->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
//	SetActuatorFunction(0, ChSharedPtr<ChFunction>(new ChFunction_Const(0)));
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
	double l_mod = l - jointClearance;

	// ** initialize U
//	CreateArm(0, l_mod, ChVector<>(-w/2 + r2, 0, l_mod/2 + r2 + jointClearance), Q_from_AngAxis(CH_C_PI / 2, VECT_Y));
//	CreateArm(1, w, ChVector<>(0, 0, 0));
//	CreateArm(2, l_mod, ChVector<>(w/2 - r2, 0, l_mod/2 + r2 + jointClearance), Q_from_AngAxis(CH_C_PI / 2, VECT_Y));

	// ** initialize straight
	CreateArm(0, l_mod, ChVector<>(-w/2 - jointClearance - r2 - l_mod/2, 0, 0));
	CreateArm(1, w, ChVector<>(0, 0, 0));
	CreateArm(2, l_mod, ChVector<>(w/2 + jointClearance + r2 + l_mod/2, 0, 0));

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

void Smarticle::SetAngle(double mangle1, double mangle2, bool degrees = false)
{
	if (degrees)
	{
		angle1 = mangle1*CH_C_PI / 180.0;
		angle2 = mangle2*CH_C_PI / 180.0;
	}
	else
	{
		angle1 = mangle1;
		angle2 = mangle2;
	}
}
void Smarticle::SetAngle(double mangle, bool degrees = false)
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
void Smarticle::SetAngle1(double mangle1, bool degrees = false)
{
	if (degrees) { angle1 = mangle1*CH_C_PI / 180.0; }
	else{ angle1 = mangle1; }
}
void Smarticle::SetAngle2(double mangle2, bool degrees = false)
{
	if (degrees) { angle2 = mangle2*CH_C_PI / 180.0; }
	else{ angle2 = mangle2; }
}

double Smarticle::GetAngle1(bool degrees = true)
{
	if (degrees)
		return angle1*180.0 / CH_C_PI;
	else
		return angle1;
}
double Smarticle::GetAngle2(bool degrees = true)
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
	s_motionPiece->SetMotionSegment(numAddedMotions);
	motion_vector.push_back(s_motionPiece);
}

void Smarticle::UpdateSmarticleMotion() {
	double currentTime = m_system->GetChTime();
	int currentSegment = current_motion->GetMotionSegment();
	double t_in = current_motion->startTime - currentTime;
	if (  (t_in >= 0 && t_in < current_motion->timeInterval) &&
			(currentSegment < motion_vector.size() - 1) ) {
		current_motion = motion_vector[currentSegment + 1];
	}
}

void Smarticle::UpdateSmarticleMotionLoop() {
	int currentSegment = current_motion->GetMotionSegment();
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

void Smarticle::MoveSquare() {
	double ang01 = link_actuator01->Get_mot_rot();
	double ang12 = link_actuator12->Get_mot_rot();
	double omega1 = current_motion->joint_01.omega;
	double omega2 = current_motion->joint_12.omega;


	int low1 = ang01 < current_motion->joint_01.theta1;
	int high1 = ang01 > current_motion->joint_01.theta2;
	int low2 = ang12 < current_motion->joint_12.theta1;
	int high2 = ang12 > current_motion->joint_12.theta2;

	if (low1 && !(low2) && !(high2)) {
		current_motion->motionSubSegment = 3;
	} else if (high1 && !(low2) && !(high2)) {
		current_motion->motionSubSegment = 1;
	} else if (low2 && !(low1) && !(high1)) {
		current_motion->motionSubSegment = 0;
	} else if (high2 && !(low1) && !(high1)) {
		current_motion->motionSubSegment = 2;
	} else if (!(low1) && !(high1) && !(low2) && !(high2)) {
		current_motion->motionSubSegment = 0;
	} else {
		current_motion->motionSubSegment = (current_motion->motionSubSegment + 1) % 4;
	}

	switch (current_motion->motionSubSegment) {
	case 0:
		current_motion->joint_01.omega = defaultOmega;
		current_motion->joint_12.omega = 0;
		break;
	case 1:
		current_motion->joint_01.omega = 0;
		current_motion->joint_12.omega = defaultOmega;
		break;
	case 2:
		current_motion->joint_01.omega = -defaultOmega;
		current_motion->joint_12.omega = 0;
		break;
	case 3:
		current_motion->joint_01.omega = 0;
		current_motion->joint_12.omega = -defaultOmega;
		break;
	}

//	if (
//			ang01 > current_motion->joint_01.theta1 &&
//			ang01 < current_motion->joint_01.theta2 &&
//			ang12 > current_motion->joint_12.theta1 &&
//			ang12 < current_motion->joint_12.theta2)

//	if ((ang01 < current_motion->joint_01.theta1) && (omega1 < 0)) {
//		omega1 *= -1;
//	}
//	if ((ang01 > current_motion->joint_01.theta2) && (omega1 > 0)) {
//		omega1 *= -1;
//	}
//	if ((ang12 < current_motion->joint_12.theta1) && (omega2 < 0)) {
//		omega2 *= -1;
//	}
//	if ((ang12 > current_motion->joint_12.theta2) && (omega2 > 0)) {
//		omega2 *= -1;
//	}
//
//	current_motion->joint_01.omega = omega1;
//	current_motion->joint_12.omega = omega2;
//
//	this->SetActuatorFunction(0, current_motion->joint_01.omega);
//	this->SetActuatorFunction(1, current_motion->joint_12.omega);
//
//
//	printf("move square, omegas %f %f\n", current_motion->joint_01.omega, current_motion->joint_12.omega);

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

	current_motion = motion_vector[1];
	MotionType sMotion = current_motion->GetMotionType();
	switch (current_motion->GetMotionType()) {
	case SQUARE_G:
		MoveSquare();
		break;
	case CIRCLE_G:
//		MoveCircle();
		break;
	case RELEASE_G:
//		MoveRelease();
		break;
	default:
		break;
	}
}

ChSharedPtr<SmarticleMotionPiece> Smarticle::Get_Current_Motion() {
	return current_motion;
}

void SetActuatorFunction(int actuatorID, ChSharedPtr<ChFunction> actuatorFunction);

void SetActuatorFunction(int actuatorID, ChSharedPtr<ChFunction> actuatorFunction);

void Smarticle::SetBodyFixed(bool mev){
	arm0->SetBodyFixed(mev);
	arm1->SetBodyFixed(mev);
	arm2->SetBodyFixed(mev);
}
