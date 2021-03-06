/*
 * Smarticle.cpp
 *
 *  Created on: May 22, 2015
 *      Author: Arman Pazouki
 */

#include "Smarticle.h"
#include "utils/ChUtilsGeometry.h"
#include "utils/ChUtilsCreators.h"
#include <math.h>       /* round, floor, ceil, trunc */



 //#include "physics/ChSystem.h"  // Arman: take care of this later
 //#include "chrono_parallel/physics/ChSystemParallel.h"
 //#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"


using namespace chrono;

double armt2 = .0032 / 2 * sizeScale; //8.06 mm with solar 3.2 without
Smarticle::Smarticle(
	std::shared_ptr<CH_SYSTEM> otherSystem
) : m_system(otherSystem) {
	smarticleID = -1;
	arm_density = 7800;
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
std::vector<std::pair<double, double>> Smarticle::extra1;
std::vector<std::pair<double, double>> Smarticle::extra2;
std::vector<std::pair<double, double>> Smarticle::midTorque;


std::shared_ptr<ChTexture> Smarticle::mtextureOT = std::make_shared<ChTexture>();
std::shared_ptr<ChTexture> Smarticle::mtextureArm = std::make_shared<ChTexture>();
std::shared_ptr<ChTexture> Smarticle::mtextureMid = std::make_shared<ChTexture>();
double Smarticle::pctActive = 1.0;
double Smarticle::distThresh;
unsigned int Smarticle::global_GUI_value;
Smarticle::~Smarticle()
{
	for (size_t i = 0; i < numEngs; i++)
	{
		m_system->RemoveLink(this->link_actuators[i]);
	}
	for (size_t i = 0; i < numSegs; i++)
	{
		m_system->RemoveBody(this->arms[i]);
	}
	if (armsController)
		armsController->~Controller();

	GetLog() << "Deleting Smarticle\n";
}



void Smarticle::Properties(
	int sID,
	double other_density,
	std::shared_ptr<ChMaterialSurface> surfaceMaterial,
	double other_envelope,
	double other_l,
	double other_w,
	double other_r,
	double other_r2,
	ChVector<> pos,
	ChQuaternion<> rot,
	double other_angle,
	double other_angle2) {

	smarticleID = sID;
	density = other_density;
	mat_smarts = surfaceMaterial;
	l = other_l;
	w = other_w;
	r = other_r;
	r2 = other_r2;
	collisionEnvelope = other_envelope;
	initPos = pos;
	rotation = rot;
	angles[0] = other_angle;
	angles[1] = other_angle2;
	volume = GetVolume();



	mtextureOT->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));
	mtextureArm->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));
	mtextureMid->SetTextureFilename(GetChronoDataFile("cubetexture_blue_bordersBlueOrientedARROWS.png"));
}

void Smarticle::Properties(
	int sID,
	double other_density,
	std::shared_ptr<ChMaterialSurface> surfaceMaterial,
	double other_envelope,
	double other_l,
	double other_w,
	double other_r,
	double other_r2,
	double other_omega,
	ChVector<> pos,
	ChQuaternion<> rot,
	double other_angle,
	double other_angle2) {

	Properties(sID, other_density, surfaceMaterial, other_envelope, other_l, other_w, other_r, other_r2, pos, rot, other_angle, other_angle2);
	SetDefaultOmega(other_omega);
	//SetOmega(other_omega, true);
}
/////Will added Smarticle properties///////
void Smarticle::Properties(
	int sID,
	int mdumID,
	double otherArm_density,
	double otherMid_density,
	std::shared_ptr<ChMaterialSurface> surfaceMaterial,
	double other_envelope,
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
	double other_angHigh) {

	Properties(sID, otherArm_density, surfaceMaterial, other_envelope, other_l, other_w, other_r, other_r2, pos, rot, other_angle, other_angle2);
	SetDefaultOmega(other_omega);
	//SetOmega1(other_omega);
	//SetOmega2(other_omega);
	SetOmega(other_omega);
	moveType = (MoveType)global_GUI_value;
	prevMoveType = (MoveType)global_GUI_value;
	//moveTypeIdxs.resize(MoveType::OT, 0);
	moveTypeIdxs.resize(MoveType::OT + 1, 0);
	arm0OT = false;
	arm2OT = false;
	torqueLimit = other_torThresh2;
	angLow = 0;
	angHigh = 120;
	dumID = mdumID;
	armBroken = false;
	arm_density = otherArm_density;
	mid_density = otherMid_density;
	std::tuple<double, double, double, double> a(0.0, 0.0, 0.0, 0.0);
	//torques = { a, a, a, a, a, a, a }; //probably a better way to do this....
	torques = { a }; //probably a better way to do this....
	nextOmega = { 0, 0 };
	nextAngle = { 0,0 };
	currTorque = { 0, 0 };
	torqueAvg = std::make_tuple(0, 0, 0, 0);
	initialAngs[0] = other_angle;
	initialAngs[1] = other_angle2;

	activateStress = .05; //.05 //reassign during constructor

	OTThresh = 1.01;//for single arm 
	MTThresh = OTThresh / 3.0;//for both arms?
	LTThresh = .01;//for both arms?

	//initialize OT timer params
	this->OTMaxTime = .30;
	//give random time so no all change at same time
	this->OTTimer = genRand(OTMaxTime);

	this->OTRunning = false;
	this->OTVal.emplace_back(GUI1);
	this->OTVal.emplace_back(EXTRA1);
	this->OTVal.emplace_back(EXTRA2);
	this->OTValIdx = genRandInt(0, (int)OTVal.size() - 1);
	//GetLog() << "smartRandidx:" << this->OTValIdx << "\n";

	//initialize LT timer params

	this->LTMaxTime = .10;
	this->LTTimer = genRand(LTMaxTime);;
	this->LTRunning = false;
	this->LTVal.emplace_back(GUI2);
	//this->LTVal.emplace_back(GUI3);
	this->LTValIdx = genRandInt(0, (int)LTVal.size() - 1);
	for (int i = 0; i < numSegs; i++)
	{
		textureAssets[i] = std::make_shared<ChTexture>();
		//if (i == 1)
		//{
		//	if (stapleSize) { textureAssets[i]->SetTextureFilename(GetChronoDataFile("cubetexture_blue_bordersBlueOrientedARROWS.png")); }
		//	else { textureAssets[i]->SetTextureFilename(GetChronoDataFile("cubetexture_SmarticlePicture.png")); }
		//}
		//else
		//{
		//	if (stapleSize) { textureAssets[i]->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png")); }
		//	else { textureAssets[i]->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png")); }
		//}
	}
	setTextures();

	for (int i = 0; i < numEngs; i++)
	{
		armOT[i] = false;
	}
	//create texture loop
}

void Smarticle::setTextures()
{
	
	for (int i = 0; i < numSegs; i++)
	{
		if (active)
		{
			if (i == 1)
			{
				if (stapleSize) { textureAssets[i]->SetTextureFilename(GetChronoDataFile("cubetexture_blue_bordersBlueOrientedARROWS.png")); }
				else { textureAssets[i]->SetTextureFilename(GetChronoDataFile("cubetexture_SmarticlePicture.png")); }
			}
			else
			{
				if (stapleSize) { textureAssets[i]->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png")); }
				else { textureAssets[i]->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png")); }
			}
		}
		else
		{	textureAssets[i]->SetTextureFilename(GetChronoDataFile("cubetexture_white_bordersBlack.png"));}
	}
	
}

void Smarticle::updateTorqueDeque()
{
	std::tuple < double, double, double, double > oldT = torques.back();
	torques.pop_back();
	torques.emplace_front(getLinkActuator(0)->Get_mot_torque(), getLinkActuator(1)->Get_mot_torque(), getLinkActuator(0)->Get_mot_rot_dt(), getLinkActuator(1)->Get_mot_rot_dt());
	updateTorqueAvg(oldT);


}
void Smarticle::updateTorqueAvg(std::tuple <double, double, double, double > oldT)
{ //already assuming it is an avg:
	size_t len = torques.size();
	//if statement fills up torque list if not enough frames have passed
	if (steps < torques.size()) //since first frame this method runs on = 1 use <=
	{
		std::get<0>(torqueAvg) = (std::get<0>(torques.front()) + std::get<0>(torqueAvg) * steps) / (steps + 1);
		std::get<1>(torqueAvg) = (std::get<1>(torques.front()) + std::get<1>(torqueAvg) * steps) / (steps + 1);
		std::get<2>(torqueAvg) = (std::get<2>(torques.front()) + std::get<2>(torqueAvg) * steps) / (steps + 1);
		std::get<3>(torqueAvg) = (std::get<3>(torques.front()) + std::get<3>(torqueAvg) * steps) / (steps + 1);
	}
	else
	{
		std::get<0>(torqueAvg) = std::get<0>(torqueAvg) + (std::get<0>(torques.front()) - std::get<0>(oldT)) / len;
		std::get<1>(torqueAvg) = std::get<1>(torqueAvg) + (std::get<1>(torques.front()) - std::get<1>(oldT)) / len;
		std::get<2>(torqueAvg) = std::get<2>(torqueAvg) + (std::get<2>(torques.front()) - std::get<2>(oldT)) / len;
		std::get<3>(torqueAvg) = std::get<3>(torqueAvg) + (std::get<3>(torques.front()) - std::get<3>(oldT)) / len;
	}
}

//////////////////////////////////////////////
void Smarticle::SetDefaultOmega(double omega) {
	defaultOmega = omega;
}

void Smarticle::SetOmega(int idx, double momega, bool angularFreq)
{
	if (angularFreq)
		omegas[idx] = momega;
	else
		omegas[idx] = momega * 2 * PI;
}
void Smarticle::SetOmega(double momega, bool angularFreq)
{
	for (int i = 0; i < numEngs; i++)
	{
		SetOmega(i, momega, angularFreq);
	}
}

double Smarticle::GetOmega(int id, bool angularFreq)
{
	if (angularFreq)
		return omegas[id];
	else
		return omegas[id] * 2 * PI;

}
double Smarticle::GetActuatorOmega(int id)
{
	return getLinkActuator(id)->Get_mot_rot_dt();
}
double Smarticle::GetNextOmega(int id)
{
	nextOmega.at(id) = this->ChooseOmegaAmount(GetOmega(id), GetCurrAngle(id), GetNextAngle(id));
	return nextOmega.at(id);
}
void Smarticle::CreateArms(int armID, double len, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChMatrix33<> gyr;  	// components gyration
	double vol = 0;			// components volume
	double dens = 0;
	double m = 0;

	ChVector<> armDims = ChVector<>(len / 2.0, r, r2); //all funcs take half lengths
	vol = utils::CalcBoxVolume(armDims);
	gyr = utils::CalcBoxGyration(armDims);
	
	// create body, set position and rotation, add surface property, and clear/make collision model
	auto arm = std::make_shared<ChBody>();

	//$$$$$$$$$$$
	//$$$$$$$$$$$

	ChVector<> posArm = rotation.Rotate(posRel) + initPos;
	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);

	arm->SetBodyFixed(false);
	arm->GetPhysicsItem()->SetIdentifier(dumID + armID);
	if (armID == 1)
	{
		arm->SetBodyFixed(false);
		dens = mid_density;
	}
	else
	{
		arm->SetBodyFixed(false);
		dens = arm_density;
	}
	m = dens*vol;
	arm->SetMass(m);
	arm->SetDensity((float)dens);
	arm->SetMaterialSurface(mat_smarts);


	//	double mass = density * vol;
	//double mass = .005;//.043/3.0; //robot weight 43 grams
	arm->GetCollisionModel()->ClearModel();
	//arm->SetLimitSpeed(true);
	//arm->SetMaxSpeed(PI * 2 * 5);
	//arm->SetMaxWvel(PI * 2 * 5);
	//arm->ClampSpeed();

	if (visualize)
	{
		arm->AddAsset(textureAssets[armID]);
		if(armID==1)
			arm->SetName("smarticle_cent");
		else
			arm->SetName("smarticle_arm");
	}

	arm->GetCollisionModel()->SetEnvelope(collisionEnvelope*.1);
	arm->GetCollisionModel()->SetSafeMargin(arm->GetCollisionModel()->GetEnvelope()/2.0);
	utils::AddBoxGeometry(arm.get(), armDims, ChVector<>(0, 0, 0), QUNIT, visualize);

	arm->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2

	arm->GetCollisionModel()->BuildModel(); // this function overwrites the intertia
	arm->SetCollide(true);
	// change mass and inertia property

	arm->SetInertia(gyr*m);
	//arm->SetDensity(density);

	m_system->AddBody(arm);


	arms[armID] = arm;
}
void Smarticle::CreateArm(int armID, double len, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChVector<> gyr;  	// components gyration
	double vol;			// components volume

	vol = utils::CalcBoxVolume(ChVector<>(len / 2.0, r, r2));
	gyr = utils::CalcBoxGyration(ChVector<>(len / 2.0, r, r2)).Get_Diag();
	// create body, set position and rotation, add surface property, and clear/make collision model
	auto arm = std::make_shared<ChBody>();

	//$$$$$$$$$$$
	//$$$$$$$$$$$

	ChVector<> posArm = rotation.Rotate(posRel) + initPos;
	arm->SetName("smarticle_arm");
	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);
	arm->SetCollide(true);
	arm->SetBodyFixed(false);
	arm->GetPhysicsItem()->SetIdentifier(dumID + armID);
	double m = 0;
	if (armID == 1) //this was old code from when I was fixing them to fit
	{
		arm->SetBodyFixed(false);
		m = vol*mid_density;
		arm->SetDensity((float)mid_density);
	}
	else
	{
		arm->SetBodyFixed(false);
		m = vol*arm_density;
		arm->SetDensity((float)arm_density);
	}
	arm->SetMaterialSurface(mat_smarts);


	//	double mass = density * vol;
	//double mass = .005;//.043/3.0; //robot weight 43 grams
	arm->GetCollisionModel()->ClearModel();
	//arm->SetLimitSpeed(true);
	//arm->SetMaxSpeed(PI * 2 * 5);
	//arm->SetMaxWvel(PI * 2 * 5);
	//arm->ClampSpeed();


	if (visualize)
	{
		arm->AddAsset(textureAssets[armID]);
		if (armID == 1)
		{
			arm->SetName("smarticle_cent");
		}
		else { arm->SetName("smarticle_arm"); }
	}

	arm->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(arm.get(), ChVector<>(len / 2.0, r, r2), ChVector<>(0, 0, 0), QUNIT, visualize);

	arm->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2

	arm->GetCollisionModel()->BuildModel(); // this function overwrites the intertia

																					// change mass and inertia property
	arm->SetMass(m);
	arm->SetInertiaXX(m * gyr);
	//arm->SetDensity(density);

	m_system->AddBody(arm);
	arms[armID] = arm;
}

void Smarticle::CreateArmsPill(int armID, double len, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChMatrix33<> gyr;  	// components gyration
	double vol=0;			// components volume
	double dens=0;
	double m = 0;

	//gyr = utils::CalcCapsuleGyration(r, len / 2.0,ChVector<>(0),Angle_to_Quat(ANGLE, ChVector<>(0, 0, PI_2)));
	//gyr = utils::CalcCapsuleGyration(r, len / 2.0, ChVector<>(0) );

	ChVector<> armDims = ChVector<>(len / 2.0, r, r2); //all funcs take half lengths
	
	vol = utils::CalcCylinderVolume(r, len / 2.0);
	gyr = utils::CalcCylinderGyration(r,len/2.0, ChVector<>(0),Angle_to_Quat(ANGLE, ChVector<>(0, 0, PI_2)));

	// create body, set position and rotation, add surface property, and clear/make collision model
	auto arm = std::make_shared<ChBody>();
	ChVector<> posArm = rotation.Rotate(posRel) + initPos;

	arm->SetCollide(true);
	arm->SetBodyFixed(false);
	arm->GetPhysicsItem()->SetIdentifier(dumID + armID);
	if (armID == 1)
	{
		arm->SetBodyFixed(false);
		dens = mid_density;
	}
	else
	{
		arm->SetBodyFixed(false);
		dens = arm_density;
	}
	m = dens*vol;
	arm->SetMass(m);
	arm->SetDensity((float)dens);
	arm->SetMaterialSurface(mat_smarts);


	ChVector<> endPts = ChVector<>(len / 2.0,0,0 );
	arm->GetCollisionModel()->ClearModel();
	if (visualize)
	{
		if (armID == 1)
			arm->SetName("smarticle_cent");
		else
			arm->SetName("smarticle_arm");
	}
	arm->GetCollisionModel()->SetEnvelope(collisionEnvelope);

	//arm->GetCollisionModel()->AddCylinder(r, r, len / 2.0,ChVector<>(0), Angle_to_Quat(ANGLE, ChVector<>(0, 0, PI_2)));  // radius x, radius z, height on y
	utils::AddCylinderGeometry(arm.get(), r, len / 2.0, ChVector<>(0), Angle_to_Quat(ANGLE, ChVector<>(0, 0, PI_2)), true);
	//arm->GetCollisionModel()->AddSphere(r, endPts);
	//arm->GetCollisionModel()->AddSphere(r, endPts*-1);


	arm->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2
	arm->GetCollisionModel()->BuildModel(); // this function overwrites the intertia
	
	arm->SetInertia(gyr*m);
	m_system->AddBody(arm);
	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);

	arms[armID] = arm;
}
double Smarticle::calcOverlapSphereVol(int N)//assuming overlap is equal to radius
{
		double v = 5.0 / 12.0 * PI*(r*r*r); //overlap volume of 2 sphere
		double V = 4.0 / 3.0 * PI*(r*r*r) - v; //non - overlap volume of 2 spheres
		double VV = V - v; //non overlap volume for a double overlapped volume
		return (2.0 * V + (N - 1.0)*v + (N - 2.0)*VV);
}
ChMatrix33<> Smarticle::calcOverlapSphereGyr(int N)//gyration of N spheres glued with overlap equal to radius origin is in center of sphere link length 
{
	double l = N*r+r;
	double pos = 0;
	double el = (2.0 / 5.0) * r*r; //gyr for single sphere at origin
	ChMatrix33<> posMat = ChMatrix33<>(0);
	ChMatrix33<> gyr = ChMatrix33<>(el*N);
	for (int i = 1; i <= N; i++)
	{
		pos = i*r - l / 2.0; //origin is at center of link length 
		posMat.Set33Element(0, 0, pos*pos);
		gyr = gyr + posMat;
	}

	return gyr;
}
void Smarticle::CreateArmsSphere(int armID, double len, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChMatrix33<> gyr= ChMatrix33<>(0);  	// components gyration
	double vol = 0;			// components volume
	double dens = 0;
	double m = 0;
	int N = round(len / r);
	vol = calcOverlapSphereVol(N);
	gyr = calcOverlapSphereGyr(N);
	// create body, set position and rotation, add surface property, and clear/make collision model
	auto arm = std::make_shared<ChBody>();
	ChVector<> posArm = rotation.Rotate(posRel) + initPos;

	arm->SetCollide(true);
	arm->SetBodyFixed(false);
	arm->GetPhysicsItem()->SetIdentifier(dumID + armID);
	if (armID == 1)
	{
		arm->SetBodyFixed(false);
		dens = mid_density;
	}
	else
	{
		arm->SetBodyFixed(false);
		dens = arm_density;
	}
	m = dens*vol;
	arm->SetMass(m);


	arm->SetDensity((float)dens);
	arm->SetMaterialSurface(mat_smarts);

	arm->GetCollisionModel()->ClearModel();
	if (visualize)
	{
		
		if (armID == 1)
			arm->SetName("smarticle_cent");
		else
			arm->SetName("smarticle_arm");
	}
	arm->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	double pos=0;
	for (int i = 1; i <= N; i++)
	{
		pos = i*r - len / 2.0; //origin is at center of link length 
		utils::AddSphereGeometry(arm.get(), r, ChVector<>(pos, 0, 0),QUNIT, visualize);
	}
	arm->AddAsset(textureAssets[armID]);
	arm->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2
	arm->GetCollisionModel()->BuildModel(); // this function overwrites the intertia

	arm->SetInertia(gyr*m);
	m_system->AddBody(arm);
	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);

	arms[armID] = arm;
}

void Smarticle::CreateArm2(int armID, double len, double mr, double mr2, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChVector<> gyr;  	// components gyration
	double vol;			// components volume
	vol = utils::CalcBoxVolume(ChVector<>(len / 2.0, mr, mr2));
	gyr = utils::CalcBoxGyration(ChVector<>(len / 2.0, mr, mr2)).Get_Diag();
	// create body, set position and rotation, add surface property, and clear/make collision model
	auto arm = std::make_shared<ChBody>();

	ChVector<> posArm = rotation.Rotate(posRel) + initPos;

	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);
	arm->SetCollide(true);
	arm->SetBodyFixed(false);
	arm->GetPhysicsItem()->SetIdentifier(dumID + armID);

	double m = 0;
	if (armID == 1) //this was old code from when I was fixing them to fit
	{
		m = vol*mid_density;
		arm->SetDensity((float)mid_density);
		arm->SetBodyFixed(false);
	}
	else
	{
		m = vol*arm_density;
		arm->SetDensity((float)arm_density);
		arm->SetBodyFixed(false);
	}
	//mat_g->SetFriction(.05);
	arm->SetMaterialSurface(mat_smarts);
	//double mass = .005;//.043/3.0; //robot weight 43 grams
	arm->GetCollisionModel()->ClearModel();

	if (visualize)
	{
		arm->AddAsset(textureAssets[armID]);
		if(armID==1)
			arm->SetName("smarticle_cent");
		else
			arm->SetName("smarticle_arm");
	}
	arm->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(arm.get(), ChVector<>(len / 2.0, mr,  mr2), ChVector<>(0, 0, 0), QUNIT, visualize);

	arm->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2
	arm->GetCollisionModel()->BuildModel(); // this function overwrites the intertia

	// change mass and inertia property
	arm->SetMass(m);
	arm->SetInertiaXX(m * gyr);
	//arm->SetDensity(density);

	m_system->AddBody(arm);


	arms[armID] = arm;
}
void Smarticle::CreateArm3(int armID, double len, double mr, double mr2, ChVector<> posRel, ChQuaternion<> armRelativeRot) {
	ChVector<> gyr;  	// components gyration
	double vol;			// components volume
	vol = utils::CalcBoxVolume(ChVector<>(len / 2.0, mr, mr2));
	gyr = utils::CalcBoxGyration(ChVector<>(len / 2.0, mr, mr2)).Get_Diag();
	// create body, set position and rotation, add surface property, and clear/make collision model
	auto arm = std::make_shared<ChBody>();

	ChVector<> posArm = rotation.Rotate(posRel) + initPos;

	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);
	arm->SetCollide(true);
	arm->SetBodyFixed(false);
	arm->GetPhysicsItem()->SetIdentifier(dumID + armID);

	double m = 0;
	if (armID == 1) //this was old code from when I was fixing them to fit
	{
		m = vol*mid_density;
		arm->SetDensity((float)mid_density);
		arm->SetBodyFixed(false);
	}
	else
	{
		m = vol*arm_density;
		arm->SetDensity((float)arm_density);
		arm->SetBodyFixed(false);
	}
	//mat_g->SetFriction(.05);
	arm->SetMaterialSurface(mat_smarts);
	//double mass = .005;//.043/3.0; //robot weight 43 grams
	arm->GetCollisionModel()->ClearModel();

	if (visualize)
	{
		setTextures();
		arm->AddAsset(textureAssets[armID]);
		if (armID == 1)
			arm->SetName("smarticle_cent");
		else
			arm->SetName("smarticle_arm");
		arm->GetCollisionModel()->SetEnvelope(collisionEnvelope);

		switch (armID) {
		case 0:
		{
			utils::AddBoxGeometry(arm.get(), ChVector<>(len / 2.0, mr, mr2), ChVector<>(0, 0, 0), QUNIT, visualize);
			break;
		}
		case 1:
		{
			utils::AddBoxGeometry(arm.get(), ChVector<>(len / 2.0, mr, mr2), ChVector<>(0, 0, 0), QUNIT, visualize);
			//radius,h,pos,rot,vis
			auto a = Angle_to_Quat(ANGLE, ChVector<>(0, PI / 2, PI / 2));
			//utils::AddCylinderGeometry(arm.get(), .004, 0.005, ChVector<>(-len / 2.0 + .01131, -mr + .0085, mr2), a, true);
			//utils::AddCylinderGeometry(arm.get(), .004, 0.005, ChVector<>(-len / 2.0 + .03406, -mr + .0055, -mr2), a, true);
			break;
		}
		case 2:
		{
			utils::AddBoxGeometry(arm.get(), ChVector<>(len / 2.0, mr, mr2), ChVector<>(0, 0, 0), QUNIT, visualize);
			break;
		}
		default:
		{
			std::cerr << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
			break;
		}
		}
	}


	arm->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2
	arm->GetCollisionModel()->BuildModel(); // this function overwrites the intertia

																					// change mass and inertia property
	arm->SetMass(m);
	arm->SetInertiaXX(m * gyr);
	//arm->SetDensity(density);

	m_system->AddBody(arm);


	arms[armID] = arm;
	
}
std::shared_ptr<ChBody> Smarticle::GetArm(int armID) {
	return arms[armID];
}

std::shared_ptr<ChLinkLockRevolute> Smarticle::GetRevoluteJoint(int jointID) {
	return link_revolutes[jointID];
}
void Smarticle::SetSpeed(ChVector<> newSpeed)
{
	for (size_t i = 0; i < numSegs; i++)
	{
		arms[i]->SetPos_dt(newSpeed);
	}
	//arm0->SetPos_dt(newSpeed);
	//arm1->SetPos_dt(newSpeed);
	//arm2->SetPos_dt(newSpeed);
}
void Smarticle::SetAcc(ChVector<> newAcc)
{
	for (size_t i = 0; i < numSegs; i++)
	{
		arms[i]->SetPos_dtdt(newAcc);
	}
}
void Smarticle::SetRotSpeed(ChQuaternion<> newSpeed)
{
	for (size_t i = 0; i < numSegs; i++)
	{
		arms[i]->SetRot_dt(newSpeed);
	}
	//arms[1]->SetRot_dt(newSpeed);
}
void Smarticle::RotateSmarticle(ChQuaternion<> newRotation)
{
	GetLog() << arms[1]->GetRotAxis();
	for (size_t i = 0; i < numSegs; i++)
	{
		arms[i]->SetRot(newRotation);
	}
}

void Smarticle::UpdateState()
{
	for (int i = 0; i < numSegs; i++)
	{
		GetArm(i)->Update();
	}
}
void Smarticle::SetEdges()
{
	UpdateState();
	//ChVector<>pos(mySmarticlesVec[0]->GetArm(1)->TransformPointLocalToParent(ChVector<>(w_smarticle / 2.0, 0, 0)));
	//ChVector<>pos2(mySmarticlesVec[0]->GetArm(2)->TransformPointLocalToParent(ChVector<>(l_smarticle / 2, 0, 0)));

	//vector3df armpos(pos.x(), pos.y(), pos.z());
	//vector3df armpos2(pos2.x(), pos2.y(), pos2.z());
	//application.GetVideoDriver()->setTransform(irr::video::ETS_WORLD, core::IdentityMatrix);
	//application.GetVideoDriver()->draw3DLine(armpos - vector3df(0, 10, 0),
	//	armpos + vector3df(0, 10, 0), irr::video::SColor(70, 255, 0, 0));

	//CreateArm2(0, l, armt, armt2, ChVector<>((-w / 2.0 + (jointClearance)-cos(-angle1)*l / 2), 0, -(l / 2.0)*sin(-angle1) - offPlaneoffset), quat0);
	//CreateArm2(1, w, r, r2, ChVector<>(0, 0, 0));
	//CreateArm2(2, l, armt, armt2, ChVector<>((w / 2.0 - (jointClearance)+cos(-angle2)*l / 2), 0, -(l / 2.0)*sin(-angle2) - offPlaneoffset), quat2);


	double armt = r;
	//double armt2 = .0032 / 2 * sizeScale; //8.06 mm with solar 3.2 without
	if (!stapleSize)
	{
		//verts

		armVerts[1][0] = GetArm(1)->TransformPointLocalToParent(ChVector<>(-w / 2.0, 0, r2));  //check
		armVerts[1][1] = GetArm(1)->TransformPointLocalToParent(ChVector<>(-w / 2.0, 0, -r2));
		armVerts[1][2] = GetArm(1)->TransformPointLocalToParent(ChVector<>(w / 2.0, 0, -r2));
		armVerts[1][3] = GetArm(1)->TransformPointLocalToParent(ChVector<>(w / 2.0, 0, r2));


		armVerts[0][0] = GetArm(0)->TransformPointLocalToParent(ChVector<>(-l / 2.0, 0, r2)); //armt2
		armVerts[0][1] = GetArm(0)->TransformPointLocalToParent(ChVector<>(-l / 2.0, 0, -r2));
		armVerts[0][2] = armVerts[1][1];
		armVerts[0][3] = armVerts[1][0];
		//armVerts[0][2] = GetArm(0)->TransformPointLocalToParent(ChVector<>(l / 2.0 + jointClearance, 0, -r2));
		//armVerts[0][3] = GetArm(0)->TransformPointLocalToParent(ChVector<>(l / 2.0 + jointClearance, 0, r2));



		//armVerts[2][0] = GetArm(2)->TransformPointLocalToParent(ChVector<>(-l / 2.0 + jointClearance, 0, r2));//armt2
		//armVerts[2][1] = GetArm(2)->TransformPointLocalToParent(ChVector<>(-l / 2.0 + jointClearance, 0, -r2));
		armVerts[2][0] = armVerts[1][3];
		armVerts[2][1] = armVerts[1][2];
		armVerts[2][2] = GetArm(2)->TransformPointLocalToParent(ChVector<>(l / 2.0, 0, -r2)); //checked
		armVerts[2][3] = GetArm(2)->TransformPointLocalToParent(ChVector<>(l / 2.0, 0, r2));

		//going clockwise a->b = b-a

		//arm0
		arm0Front = (armVerts[0][0]) - (armVerts[0][3]);
		arm0OuterEdge = (armVerts[0][1]) - (armVerts[0][0]);
		arm0Back = (armVerts[0][2]) - (armVerts[0][1]);
		arm0Edge = (armVerts[0][3]) - (armVerts[0][2]);

		arm0Front = ChVector<>(-arm0Front.y(), arm0Front.x(), arm0Front.z());
		arm0OuterEdge = ChVector<>(-arm0OuterEdge.y(), arm0OuterEdge.x(), arm0OuterEdge.z());
		arm0Back = ChVector<>(-arm0Back.y(), arm0Back.x(), arm0Back.z());
		arm0Edge = ChVector<>(-arm0Edge.y(), arm0Edge.x(), arm0Edge.z());

		//arm1
		arm1Front = (armVerts[1][0]) - (armVerts[1][3]);
		arm10Shared = (armVerts[1][1]) - (armVerts[1][0]);
		arm1Back = (armVerts[1][2]) - (armVerts[1][1]);
		arm12Shared = (armVerts[1][3]) - (armVerts[1][2]);

		arm1Front = ChVector<>(-arm1Front.y(), arm1Front.x(), arm1Front.z());
		arm10Shared = ChVector<>(-arm10Shared.y(), arm10Shared.x(), arm10Shared.z());
		arm1Back = ChVector<>(-arm1Back.y(), arm1Back.x(), arm1Back.z());
		arm12Shared = ChVector<>(-arm12Shared.y(), arm12Shared.x(), arm12Shared.z());



		//arm2
		arm2Front = (armVerts[2][0]) - (armVerts[2][3]);
		arm2Edge = (armVerts[2][1]) - (armVerts[2][0]);
		arm2Back = (armVerts[2][2]) - (armVerts[2][1]);
		arm2OuterEdge = (armVerts[2][3]) - (armVerts[2][2]);

		arm2Front = ChVector<>(-arm2Front.y(), arm2Front.x(), arm2Front.z());
		arm2Edge = ChVector<>(-arm2Edge.y(), arm2Edge.x(), arm2Edge.z());
		arm2Back = ChVector<>(-arm2Back.y(), arm2Back.x(), arm2Back.z());
		arm2OuterEdge = ChVector<>(-arm2OuterEdge.y(), arm2OuterEdge.x(), arm2OuterEdge.z());


		armAxes[0][0] = arm0Front;
		armAxes[0][1] = arm0OuterEdge;
		armAxes[0][2] = arm0Back;
		armAxes[0][3] = arm0Edge;

		armAxes[1][0] = arm1Front;
		armAxes[1][1] = arm10Shared;
		armAxes[1][2] = arm1Back;
		armAxes[1][3] = arm12Shared;

		armAxes[2][0] = arm2Front;
		armAxes[2][1] = arm2Edge;
		armAxes[2][2] = arm2Back;
		armAxes[2][3] = arm2OuterEdge;

	}

}

void Smarticle::RotateSmarticleBy(ChQuaternion<> newAng)
{
	ChQuaternion<> quat0 = Angle_to_Quat(ANGLE, ChVector<>(0, -angles[0], 0));
	ChQuaternion<> quat2 = Angle_to_Quat(ANGLE, ChVector<>(0, angles[1], 0));
	quat0.Normalize();
	quat2.Normalize();
	newAng.Normalize();

	auto c = arms[1]->GetCoord();
	c.pos.z() = c.pos.z() + 1;
	//ChTransform<>::TransformLocalToParent()

	auto ca = arms[1]->GetRotAxis();
	rotation = rotation*newAng;
	arms[1]->SetRot(rotation);
	UpdateState();
	arms[1]->GetCollisionModel()->SyncPosition();

	auto a = arms[0]->GetCoord();
	a.ConcatenatePreTransformation(arms[1]->GetCoord());
	arms[0]->SetCoord(a);
	arms[2]->SetCoord(a);

	UpdateState();

	arms[0]->GetCollisionModel()->SyncPosition();
	arms[2]->GetCollisionModel()->SyncPosition();

	//rotation = newAng;
	//ChVector<>posRel = ChVector<>((-w / 2.0 + (jointClearance)-cos(-angle1)*l / 2), 0, -(l / 2.0)*sin(-angle1) - offPlaneoffset);
	//ChQuaternion<> armRelativeRot = quat0;
	//arm0->SetPos(rotation.Rotate(posRel) + arm0->GetPos());
	//arm0->SetRot(rotation*armRelativeRot);

	//posRel = ChVector<>(0, 0, 0);
	//armRelativeRot = QUNIT;
	//arm1->SetPos(rotation.Rotate(posRel) + arm1->GetPos());
	//arm1->SetRot(rotation*armRelativeRot);

	//
	//posRel = ChVector<>((w / 2.0 - (jointClearance)+cos(-angle2)*l / 2), 0, -(l / 2.0)*sin(-angle2) - offPlaneoffset);
	//armRelativeRot = quat2;
	//arm2->SetPos(rotation.Rotate(posRel) + arm2->GetPos());
	//arm2->SetRot(rotation*armRelativeRot);
	//
	//arm0->SetRot_dt(QUNIT);
	//arm1->SetRot_dt(QUNIT);
	//arm2->SetRot_dt(QUNIT);
	//UpdateState();
	//arm0->GetCollisionModel()->SyncPosition();
	//arm1->GetCollisionModel()->SyncPosition();
	//arm2->GetCollisionModel()->SyncPosition();
	//link_actuator01->SyncCollisionModels();
	//link_actuator12->SyncCollisionModels();


	/*CreateArm2(0, l, armt, armt2, ChVector<>((-w / 2.0 + (jointClearance)-cos(-angle1)*l / 2), 0, -(l / 2.0)*sin(-angle1) - offPlaneoffset), quat0
	GetLog() << "arm1:" << arm1->GetPos();
	GetLog() << "arm2:" << arm2->GetPos();
	double l_mod = l + 2 * r2 - jointClearance;
	GetLog() << "\nl-mod" << l_mod << " w:" << w << " l:"<<l<<"\n";
	GetLog() << "local to parent" << arm1->TransformPointLocalToParent(ChVector<>(w / 2.0 - (jointClearance)+cos(-angle2)*l / 2, 0, -(l / 2.0)*sin(-angle2) - offPlaneoffset));
	GetLog() << "local to parent" << arm1->TransformPointLocalToParent(ChVector<>(w / 2.0+l, 0,0));
	GetLog() << "local to parent" << arm2->TransformPointLocalToParent(ChVector<>(l, 0, 0));*/
}
void Smarticle::TransportSmarticle(ChVector<> newPosition)
{
	arms[0]->SetPos(arms[0]->GetPos() - arms[1]->GetPos() + newPosition);
	arms[2]->SetPos(arms[2]->GetPos() - arms[1]->GetPos() + newPosition);
	arms[1]->SetPos(newPosition);
	
	for (size_t i = 0; i < numSegs; i++)
	{arms[i]->GetCollisionModel()->SyncPosition();}
	
	//arm0->SetPos(arm0->GetPos() - arm1->GetPos() + newPosition);
	//arm2->SetPos(arm2->GetPos() - arm1->GetPos() + newPosition);
	//arm1->SetPos(newPosition);

	//arm0->GetCollisionModel()->SyncPosition();
	//arm1->GetCollisionModel()->SyncPosition();
	//arm2->GetCollisionModel()->SyncPosition();
}
void Smarticle::CreateJoints() {
	// link 1
	std::shared_ptr<ChLinkLockRevolute> link_revolute01 = std::make_shared<ChLinkLockRevolute>();
	std::shared_ptr<ChLinkLockRevolute> link_revolute12 = std::make_shared<ChLinkLockRevolute>();
	// ChVector<> pR01(-w / 2.0+r2, 0, 0);
	// ChVector<> pR12(w / 2.0-r2, 0, 0);
	ChVector<> pR01(-w / 2.0 - r2, 0, 0);
	ChVector<> pR12(w / 2.0 + r2, 0, 0);
	ChQuaternion<> qx = Q_from_AngAxis(PI_2, VECT_X);


		// link 1
		link_revolute01->Initialize(arms[0], arms[1], ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx));
		link_revolute01->SetMotion_axis(ChVector<>(0, 0, 1));


		// link 2
		link_revolute12->Initialize(arms[1], arms[2], ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx));
		link_revolute12->SetMotion_axis(ChVector<>(0, 0, 1));


	m_system->AddLink(link_revolute01);
	m_system->AddLink(link_revolute12);

	link_revolutes[0] = link_revolute01;
	link_revolutes[1] = link_revolute12;
}

void Smarticle::CreateActuators() {
	std::shared_ptr<ChLinkEngine> link_actuator = std::make_shared<ChLinkEngine>();

	std::shared_ptr<ChLinkEngine> link_actuator01 = std::make_shared<ChLinkEngine>();
	std::shared_ptr<ChLinkEngine> link_actuator12 = std::make_shared<ChLinkEngine>();
	ChVector<> pR01(-w / 2.0, 0, 0);
	ChVector<> pR12(w / 2.0, 0, 0);
	if (!stapleSize)
	{
		pR01 = ChVector<>(-(w / 2 - jointClearance), 0, -offPlaneoffset);
		pR12 = ChVector<>((w / 2 - jointClearance), 0, -offPlaneoffset);
	}

	ChQuaternion<> qx1 = Q_from_AngAxis(-PI_2, VECT_X);
	ChQuaternion<> qx2 = Q_from_AngAxis(PI_2, VECT_X);
	ChQuaternion<> qy1 = Q_from_AngAxis(GetAngle(0), VECT_Z);
	ChQuaternion<> qy2 = Q_from_AngAxis(GetAngle(1), VECT_Z);

	qx1.Normalize();
	qx2.Normalize();
	qy1.Normalize();
	qy2.Normalize();

	for (int i = 0; i < numEngs; i++)
	{
		std::shared_ptr<ChLinkEngine> link_actuator = std::make_shared<ChLinkEngine>();
		link_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		if (i == 0)
		{
			link_actuator->Initialize(arms[0], arms[1], false, ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx1*qy1), ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx1));
		}
		else
		{
			link_actuator->Initialize(arms[2], arms[1], false, ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx2*qy2), ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx2));
		}
		link_actuator->GetLimit_Rz()->Set_min(D2R*angLow);
		link_actuator->GetLimit_Rz()->Set_max(D2R*angHigh);
		link_actuator->GetLimit_Rz()->Set_active(true);
		link_actuator->SetMotion_axis(ChVector<>(0, 0, 1));
		link_actuator->Set_mot_inertia(arms[0]->GetInertiaXX().z());
		link_actuators[i] = link_actuator;
		m_system->AddLink(link_actuators[i]);

	}

	//link_actuator01->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	//link_actuator12->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	//	//link_actuator01->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
	//	//link_actuator12->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);

	//	

	//	//link_actuator01->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
	//	//link_actuator12->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
	//	link_actuator01->Initialize(arms[0], arms[1], false, ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx1*qy1), ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx1));
	//	link_actuator01->GetLimit_Rz()->Set_min(D2R*angLow);
	//	link_actuator01->GetLimit_Rz()->Set_max(D2R*angHigh);
	//	link_actuator01->GetLimit_Rz()->Set_active(true);
	//	link_actuator01->SetMotion_axis(ChVector<>(0, 0, 1));

	//	/*link_actuator01->GetForce_D()->Set_active(true);
	//	auto modK = (ChFunction_Const*)(link_actuator01->GetForce_D()->Get_modul_K());
	//	link_actuator01->GetForce_D()->Set_K(1);
	//	modK->Set_yconst(100);*/

	//		
	//	link_actuator12->Initialize(arms[2], arms[1], false, ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx2*qy2), ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx2));
	//	link_actuator12->GetLimit_Rz()->Set_min(D2R*angLow);
	//	link_actuator12->GetLimit_Rz()->Set_max(D2R*angHigh);
	//	link_actuator12->GetLimit_Rz()->Set_active(true);
	//	link_actuator12->SetMotion_axis(ChVector<>(0, 0, 1));
	//	//link_actuator12->Set_mot_inertia(1);

	////if (!active)
	////{
	////	link_actuator01->GetLimit_Rz()->Set_min(D2R * GetInitialAngle(0));
	////	link_actuator01->GetLimit_Rz()->Set_max(D2R * GetInitialAngle(0));
	////	link_actuator12->GetLimit_Rz()->Set_min(D2R * GetInitialAngle(1));
	////	link_actuator12->GetLimit_Rz()->Set_max(D2R * GetInitialAngle(1));
	////}

	//link_actuator01->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);
	//link_actuator12->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);
	//
	//m_system->AddLink(link_actuator01);
	//m_system->AddLink(link_actuator12);
	//
	////////////
	////actuator 1
	//link_actuators[0] = link_actuator01;
	//link_actuators[1] = link_actuator12;

}

void Smarticle::Create() {
	jointClearance = 0;
	double l_mod;
	//double l_mod = l + 2 * r2 - jointClearance;

	ChQuaternion<> quat0 = Q_from_AngAxis(angles[0], VECT_Y);
	ChQuaternion<> quat1 = Q_from_AngAxis(angles[0], VECT_Z);
	ChQuaternion<> quat2 = Q_from_AngAxis(-angles[1], VECT_Y);

	//works with box
	//ChQuaternion<> quat0 = Angle_to_Quat(ANGLE, ChVector<>(0, -angles[0], 0));
	//ChQuaternion<> quat2 = Angle_to_Quat(ANGLE, ChVector<>(0, angles[1], 0));
	quat0.Normalize();
	quat1.Normalize();
	quat2.Normalize();

	ChVector<> armDimL = ChVector<>((w*0.7) / 2.0, r, r2); //size of arm for reference smarticle size
	ChVector<> armDimW = ChVector<>(w / 2.0, r, r2); //size of spine for reference smarticle size
	if (stapleSize)
	{
		
		if (!sphereVersion)
		{
			offPlaneoffset = 0;
			l_mod = l + 2 * r2 - jointClearance;
			CreateArms(0, l_mod, ChVector<>(-w / 2.0 - (l / 2.0)*cos(-angles[0]), 0, -(l_mod / 2.0 - r2)*sin(-angles[0])), quat0);
			CreateArms(1, w, ChVector<>(0, 0, 0));
			CreateArms(2, l_mod, ChVector<>(w / 2.0 + (l / 2.0)*cos(-angles[1]), 0, -(l_mod / 2.0 - r2)*sin(-angles[1])), quat2);
			armDimL.Set((0.7*w+2*r2-jointClearance)/2.0, armDimL.y(), armDimL.z());
		}
		//l_mod = l + 2*r - jointClearance;
		//CreateArmsPill(0, l_mod, ChVector<>(-w / 2.0 - (l / 2.0)*cos(-angles[0]), 0, -(l_mod / 2.0-r)*sin(-angles[0])), quat0);
		//CreateArmsPill(1, w, ChVector<>(0, 0, 0));
		//CreateArmsPill(2, l_mod, ChVector<>(w / 2.0 + (l / 2.0)*cos(-angles[1]), 0, -(l_mod / 2.0 -r)*sin(-angles[1])), quat2);

		if (sphereVersion)
		{
			l_mod = l;
			CreateArmsSphere(0, l_mod, ChVector<>(-w / 2.0 - (l / 2.0)*cos(-angles[0]), 0, -(l_mod / 2.0)*sin(-angles[0])), quat0);
			CreateArmsSphere(1, w, ChVector<>(0, 0, 0));
			CreateArmsSphere(2, l_mod, ChVector<>(w / 2.0 + (l / 2.0)*cos(-angles[1]), 0, -(l_mod / 2.0 - r)*sin(-angles[1])), quat2);
			
		}
		//l_smarticle = w_smarticle*mult_l; //r2=t2/2.0

	}
	else
	{
		offPlaneoffset = .00825 - r2;//smarticle arms arent centered in plane, arm is offset 8.25mm from front or offPlaneoffset-t2
		jointClearance = .0065;//6.5 mm in x dir jc in y dir is r2
		double armt = r*.96;
		double armt2 = .0032 / 2 * sizeScale; //8.06 mm with solar 3.2 without

		CreateArm3(0, l, armt, armt2, ChVector<>((-w / 2.0 + (jointClearance)-cos(-angles[0])*l / 2), 0, -(l / 2.0)*sin(-angles[0]) - offPlaneoffset), quat0);
		CreateArm3(1, w, r, r2, ChVector<>(0, 0, 0));
		CreateArm3(2, l, armt, armt2, ChVector<>((w / 2.0 - (jointClearance)+cos(-angles[1])*l / 2), 0, -(l / 2.0)*sin(-angles[1]) - offPlaneoffset), quat2);
	}

	refSmartMass = (utils::CalcBoxVolume(armDimL)*arm_density) * 2 + (utils::CalcBoxVolume(armDimW)*mid_density); //doesn't work for cylindrical arms
	//need this here
	if (read_from_file <= 0)//not reading
	{
		if (active == true)
		{
			if (genRand() < pctActive)
				active = true;
			else
				active = false;
		}
	}
	CreateActuators();

	// mass property
	mass = 0;
	for (size_t i = 0; i < numSegs; i++)
	{mass = mass + arms[i]->GetMass();}
	
	armsController = (new Controller(m_system, this));
	armsController->smartStr = this->smartStr;
	armsController->outputLimit = torqueLimit;
	armsController->omegaLimit = omegaLim;
	if (!active)
	{
		arms[0]->SetName("D_smarticle_arm");
		arms[1]->SetName("D_smarticle_cent");
		arms[2]->SetName("D_smarticle_arm");

		//armsController = NULL;
		setTextures();

		//link_actuators[0]->SetBroken(true);
		//link_actuators[0]->SetDisabled(1);
		//link_actuators[1]->SetDisabled(1);
	}
	//if (active)
	//{
	//	armsController = (new Controller(m_system,this));
	//	armsController->outputLimit= torqueLimit;
	//	armsController->omegaLimit = omegaLim;
	//}
	//else
	//{

	//	arm0->SetName("D_smarticle_arm");
	//	arm1->SetName("D_smarticle_cent");
	//	arm2->SetName("D_smarticle_arm");

	//	armsController = NULL;
	//	arm0_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_white_bordersBlack.png"));
	//	arm1_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_white_bordersBlack.png"));
	//	arm2_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_white_bordersBlack.png"));
	//}
}

void Smarticle::ChangeActive(bool m_active)
{
	if (m_active != active)
	{
		if (m_active)
		{
			arms[0]->SetName("smarticle_arm");
			arms[1]->SetName("smarticle_cent");
			arms[2]->SetName("smarticle_arm");

			for (int i = 0; i < numEngs; i++)
			{
				getLinkActuator(i)->GetLimit_Rz()->Set_min(D2R * angLow);
				getLinkActuator(i)->GetLimit_Rz()->Set_max(D2R * angHigh);
			}

			
		}
		else
		{
			arms[0]->SetName("D_smarticle_arm");
			arms[1]->SetName("D_smarticle_cent");
			arms[2]->SetName("D_smarticle_arm");

			//link_actuator12->Initialize(arm2, arm1, false, ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx2*qy2), ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx2));
			for (int i = 0; i < numEngs; i++)
			{
				getLinkActuator(i)->GetLimit_Rz()->Set_min(getLinkActuator(i)->GetRelAngle());
				getLinkActuator(i)->GetLimit_Rz()->Set_max(getLinkActuator(i)->GetRelAngle());
			}
			this->armsController->resetCumError = true;
		}

	}

	this->active = m_active;
	setTextures();
}
std::shared_ptr<ChFunction> Smarticle::GetActuatorFunction(int id) {
	if (id <= numEngs - 1)
	{
		return EngFunctions[id];
	}
	//else
	//{
	//	std::cout << "Error! "<< __FUNCTION__ <<" smarticle can only have actuators with ids from {0, 1}" << std::endl;
	//}
	return std::shared_ptr<ChFunction>(NULL);
}
void Smarticle::SetActuatorFunction(int id, std::shared_ptr<ChFunction> actuatorFunction) {

	GetActuatorFunction(id) = actuatorFunction;
	getLinkActuator(id)->Set_rot_funct(GetActuatorFunction(id));
}

std::shared_ptr<ChLinkEngine> Smarticle::getLinkActuator(int id)
{
	if (id <= numEngs - 1)
	{
		return link_actuators[id];
	}
	else
	{
		std::cout << "Error!" << __FUNCTION__ << " smarticle can only have actuators with ids from {0, 1}" << std::endl;
		exit(-1);
	}
}
double Smarticle::GetVolume() {
	//	return r * r2 * (w + 2 * (l + jointClearance));
	return (2 * r) * (2 * r2)* (w + 2 * (l + 2 * r2));
}
double Smarticle::GetMass() {
	//	return r * r2 * (w + 2 * (l + jointClearance));
	return mass;
}
double Smarticle::GetDensity(int id)
{
	if (id == 1)
		return mid_density;
	else
		return arm_density;
}
ChVector<> Smarticle::Get_cm() {
	ChVector<> cm = ChVector<>(0);
	for (size_t i = 0; i < numSegs; i++)
	{	cm = cm + arms[i]->GetMass()*arms[i]->GetPos();}
	return cm / mass;
	//return (arm0->GetMass() * arm0->GetPos() + arm1->GetMass() * arm1->GetPos() + arm2->GetMass() * arm2->GetPos()) / mass;

}
ChVector<> Smarticle::Get_COG() {
	ChVector<> cg = ChVector<>(0);
	for (size_t i = 0; i < numSegs; i++)
	{cg = cg + arms[i]->GetMass()*arms[i]->GetPos();}
	return cg;
	//return (arm0->GetMass() * arm0->GetPos() + arm1->GetMass() * arm1->GetPos() + arm2->GetMass() * arm2->GetPos());

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
	for (int i = 0; i < numEngs; i++)
	{
		initialAngs[i] = this->GetAngle(i);
	}

}
double Smarticle::GetInitialAngle(int id)
{
	return initialAngs[id];
}
void Smarticle::SetAngles(double mangle0, double mangle1, bool degrees)
{
	if (degrees)
	{
		angles[0] = mangle0*D2R;
		angles[1] = mangle1*D2R;
		return;
	}
	else
	{
		angles[0] = mangle0;
		angles[1] = mangle1;
	}
}
void Smarticle::SetAngle(int id, double mangle, bool degrees)
{

	if (degrees) { angles[id] = mangle*D2R; }
	else { angles[id] = mangle; }
}
void Smarticle::SetAngle(double mangle, bool degrees)
{
	for (int i = 0; i < numEngs; i++)
	{
		if (degrees) { angles[i] = mangle*D2R; }
		else { angles[i] = mangle; }
	}
}

double Smarticle::GetAngle(int id, bool degrees)
{
	if (degrees)
		return angles[id] * R2D;
	else
		return angles[id];
}

void Smarticle::addInterpolatedPathToVector(double a0i, double a1i, double a0f, double a1f)
{
	double dist0 = omegas[0] * dT;
	double dist1 = omegas[1] * dT;
	int n = int(std::max(abs((a0f - a0i) / dist0), abs((a1f - a1i) / dist1)));
	std::vector<double> a0 = linspace(a0i, a0f, n);
	std::vector<double> a1 = linspace(a1i, a1f, n);

	for (int i = 0; i < n; i++)
	{
		mv->emplace_back(a0.at(i), a1.at(i));
	}

}
std::vector<double> Smarticle::linspace(double a, double b, int n) {
	std::vector<double> vec;
	double step = (b - a) / (n - 1);

	for (int i = 0; i < n; ++i)
	{
		vec.push_back(a);
		a += step;           // could recode to better handle rounding errors
	}
	return vec;
}
bool Smarticle::NotAtDesiredPos(int id, double ang, double exp)//bad method name
{
	//GetLog() << "expAng" << id << ":" << GetExpAngle(id) << "\n    ";
	//returns true if anything else but 0 is returned from here
	return ChooseOmegaAmount(GetOmega(id), ang, exp);
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
	return this->getLinkActuator(id)->Get_mot_rot();
}
double Smarticle::GetNextAngle(int id)
{
	if (id == 0)
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
	double mdt, momega, mtorqueLimit, mangLow, mangHigh, mOmegaLim;
	smarticleMoves >>
		mdt >>
		momega >>
		mtorqueLimit >>
		mangLow >>
		mangHigh >>
		mOmegaLim;
	//printf("dt %f omega %f torqueThresh2 %f angLow %f angHigh %f", mdt, momega, mtorqueThresh2, mangLow, mangHigh);
	SetDefaultOmega(momega);
	SetOmega(momega);
	omegaLim = mOmegaLim;

	char ddCh;
	angHigh = mangHigh;
	angLow = mangLow;
	distThresh = mdt*momega;
	torqueLimit = mtorqueLimit;
	OTThresh = OTThresh*torqueLimit;
	MTThresh = MTThresh*torqueLimit;
	LTThresh = LTThresh*torqueLimit;
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
	SetAngles(0, 0, false);
	//for some reason it only works with vectors?
	ChVector<> angVals;
	smarticleMoves >> angVals.x() >> ddCh >> angVals.y() >> ddCh;
	angPair.first = angVals.x();
	angPair.second = angVals.y();
	firstAngPair.first = angVals.x();
	firstAngPair.second = angVals.y();
	//TODO need to rewrite the below way of reading file, very ugly!

	ot.emplace_back(GetAngle(0), GetAngle(1));
	ot.emplace_back(GetAngle(0), GetAngle(1));
	if (global.size() < 1)
	{
		ang1 = angPair.first;
		ang2 = angPair.second;
		global.push_back(angPair);

		//SetAngle(ang1, ang2);
		//SetAngles(ang1, ang2);
		//Global


		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x() >> ddCh >> angVals.y() >> ddCh;
			angPair.first = angVals.x();
			angPair.second = angVals.y();

			global.push_back(angPair);
			//GetLog() << angVals.x() << " " << angVals.y() << " ddch:" << ddCh << "\n";
			if (ddCh == '#')
				break;
		}
	}
	if (gui1.size() < 1)
	{
		//GUI1
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x() >> ddCh >> angVals.y() >> ddCh;
			angPair.first = angVals.x();
			angPair.second = angVals.y();

			gui1.push_back(angPair);
			//GetLog() << angVals.x() << " " << angVals.y() << " ddch:" << ddCh << "\n";
			if (ddCh == '#')
				break;
		}
	}

	if (gui2.size() < 1)
	{
		//GUI2
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x() >> ddCh >> angVals.y() >> ddCh;
			angPair.first = angVals.x();
			angPair.second = angVals.y();

			gui2.push_back(angPair);
			//GetLog() << angVals.x() << " " << angVals.y() << " ddch:" << ddCh << "\n";
			//exit(-1);
			if (ddCh == '#')
				break;
		}
	}

	if (gui3.size() < 1)
	{
		//GUI3
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x() >> ddCh >> angVals.y() >> ddCh;
			angPair.first = angVals.x();
			angPair.second = angVals.y();

			gui3.push_back(angPair);
			//GetLog() << angVals.x() << " " << angVals.y() << " ddch:" << ddCh << "\n";
			//exit(-1);
			if (ddCh == '#')
				break;
		}
	}
	if (extra1.size() < 1)
	{
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x() >> ddCh >> angVals.y() >> ddCh;
			angPair.first = angVals.x();
			angPair.second = angVals.y();

			extra1.push_back(angPair);
			//GetLog() << angVals.x() << " " << angVals.y() << " ddch:" << ddCh << "\n";
			//exit(-1);
			if (ddCh == '#')
				break;
		}
	}

	if (extra2.size() < 1)
	{
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x() >> ddCh >> angVals.y() >> ddCh;
			angPair.first = angVals.x();
			angPair.second = angVals.y();

			extra2.push_back(angPair);
			//GetLog() << angVals.x() << " " << angVals.y() << " ddch:" << ddCh << "\n";
			//exit(-1);
			if (ddCh == '#')
				break;
		}
	}
	if (midTorque.size() < 1)
	{
		//midtorque
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x() >> ddCh >> angVals.y() >> ddCh;
			angPair.first = angVals.x();
			angPair.second = angVals.y();

			midTorque.push_back(angPair);
			//GetLog() << angVals.x() << " " << angVals.y() << " ddch:" << ddCh << "\n";
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
	smarticleMoves.close();

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
	if (abs(deltaAng) > 2 * distThresh)
	{		//if destAng is larger, move with positive momega
		return sgn(deltaAng)*momega;
	}
	//if <= distThresh dont move, let omega = 0;
	return 0;
}
bool Smarticle::MoveLowStress(double tor0, double tor1, double timeSinceLastChange)
{
	double t0 = abs(tor0);
	double t1 = abs(tor1);


	if (LTRunning)
	{
		specialState = LTVal.at(LTValIdx);
		return true;
	}

	specialState = -1;
	return false;
}
bool Smarticle::MoveMidStress(double tor0, double tor1)
{
	double t0 = abs(tor0);
	double t1 = abs(tor1);

	if (t0 + t1 > MTThresh)
	{// MT<t0<OT //TODO perhaps make this function if(t0+t1>2*MT) since servo can only sense stress from both
		if (genRand() < activateStress)
		{
			specialState = MIDT;
			return true;
		}
	}
	return false;
}
bool Smarticle::MoveOverStress(double tor0, double tor1)
{

	double t0 = abs(tor0);
	double t1 = abs(tor1);


	if (OTRunning)
	{
		//specialState = OT;
		specialState = OTVal.at(OTValIdx);
		return true;
	}

	specialState = -1;
	return false;
}

//returns true if either arm is OT
bool Smarticle::ChangeArmColor(double torque01, double torque12, bool LA, bool MA, bool OA)
{
	bool LTMTcolor = false;
	//for vibration upon OT, change degreesToVibrate to amount you wish to vibrate
	double degreesToVibrate = 0;
	double moveAmt = degreesToVibrate*D2R;
	if (abs(torque01) > OTThresh)
	{
		//this->setCurrentMoveType(OT);
		//mv = &ot;
		if (!arm0OT && OA)//if not previously OT
		{
			textureAssets[0]->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));

			this->ot.clear();
			//this->ot.emplace_back(angles[0] + sign(torque01)*moveAmt, angles[1] + sign(torque12)*moveAmt);
			//this->ot.emplace_back(angles[0] + moveAmt, angles[1] + moveAmt);
			this->ot.emplace_back(angles[0], angles[1]);
			//this->ot.emplace_back(angles[0] - moveAmt, angles[1] - moveAmt);
			this->armsController->resetCumError = true;
		}
		//nothing needs to be done if prev OT
		arm0OT = true;
		//setstate OT
		//specialState = OT;
	}
	else
	{
		if (arm0OT) //it prev OT but currently not
		{
			arm0OT = false;
			if (OA)
				textureAssets[0]->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));
		}
		// nothing needs to be done if not prev OT
	}

	/////////////////////ARM2///////////////////////
	if (abs(torque12) > OTThresh)
	{
		//this->setCurrentMoveType(OT);
		//mv = &ot;
		if (!arm2OT && OA)//if not previously OT
		{
			textureAssets[2]->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));
			this->ot.clear();
			//this->ot.emplace_back(angles[0] + moveAmt, angles[1] + moveAmt);
			this->ot.emplace_back(angles[0], angles[1]);
			//this->ot.emplace_back(angles[0] - moveAmt, angles[1] - moveAmt);
			this->armsController->resetCumError = true;
		}
		arm2OT = true;
		//specialState = OT;
	}
	else
	{
		if (arm2OT) //it prev OT but currently not
		{
			arm2OT = false;
			if (OA)
				textureAssets[2]->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));
		}
		// nothing needs to be done if not prev OT
	}

	if (MA)
	{
		if (abs(torque01) + abs(torque12) > MTThresh && !arm2OT && !arm0OT)
		{
			textureAssets[0]->SetTextureFilename(GetChronoDataFile("cubetexture_orange_borderOrange.png"));
			textureAssets[2]->SetTextureFilename(GetChronoDataFile("cubetexture_orange_borderOrange.png"));
		}
	}
	if (LA)
	{
		if (abs(torque01) + abs(torque12) < LTThresh)
		{
			textureAssets[0]->SetTextureFilename(GetChronoDataFile("cubetexture_green_borderGreen.png"));
			textureAssets[2]->SetTextureFilename(GetChronoDataFile("cubetexture_green_borderGreen.png"));
		}
	}
	return false;

}
///////////////////////////////////////////////////////////////////////////
bool Smarticle::ChangeToOT(bool setOT, int index)
{
	bool LTMTcolor = false;
	//for vibration upon OT, change degreesToVibrate to amount you wish to vibrate
	double degreesToVibrate = 3;
	double moveAmt = degreesToVibrate*D2R;

	//determine if any are OT
	bool anyOT = false;
	for (int i = 0; i < numEngs; i++)
	{
		if (armOT[i])
			anyOT = true;
	}


	if (setOT && !anyOT)
	{
		this->setCurrentMoveType(OT);
		mv = &ot;

		this->ot.clear();
		this->ot.emplace_back(angles[0] + moveAmt, angles[1] + moveAmt);
		//this->ot.emplace_back(angles[0] + moveAmt, angles[1] + moveAmt);
		//this->ot.emplace_back(angles[0], angles[1]);
		//this->ot.emplace_back(angles[0] - moveAmt, angles[1] - moveAmt);
		this->armsController->resetCumError = true;
		
	}
	//do nothing if setOT and OT already exists
	if (!setOT)
	{
		anyOT = false;
		armOT[index] = setOT;
		for (int i = 0; i < numEngs; i++)
		{
			if (armOT[i])
				anyOT = true;
		}
		if (!anyOT)
		{
			this->setCurrentMoveType((MoveType)global_GUI_value);
		}
	}
	return false;


}
///////////////////////////////////////////////////////////////////////////
bool Smarticle::GetArmOT(int id)
{
	return armOT[id];
}
void Smarticle::SetArmOT(int id,bool isOT)
{
	armOT[id]=isOT;
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
		this->setCurrentMoveType(EXTRA1);
		mv = &extra1;
		break;
	case 6:
		this->setCurrentMoveType(EXTRA2);
		mv = &extra2;
		break;
	case 7:
		this->setCurrentMoveType(MIDT);
		mv = &midTorque;
		break;
	case 8:
		this->setCurrentMoveType(OT);
		mv = &ot;
		break;
	default:
		this->setCurrentMoveType(GLOBAL);
		mv = &global;
		break;
	}

}

void Smarticle::CheckLTTimer(double t1, double t2)
{
	//time between switchs = 2*LTMaxTime
	static double switchTime = 2 * LTMaxTime;
	if ((abs(t1) + abs(t2)) < LTThresh)
	{
		if (genRand() < activateStress)
		{
			LTRunning = true;
		}

	}
	if ((LTTimer > LTMaxTime)) //|| previousSuccessful
	{
		LTRunning = false;
		if (LTTimer > switchTime)
		{
			LTTimer = 0;
			LTRunning = false;
			this->LTValIdx = (LTValIdx + 1) % LTVal.size();//change OTval to next value in vector
		}

	}

	//if (LTRunning)
	//{
	LTTimer += dT;
	//}

}
void Smarticle::CheckOTTimer()
{
	if (GetArm0OT() || GetArm2OT())
	{
		if (genRand() < activateStress)
		{

			OTRunning = true;
		}

	}
	if (OTTimer > OTMaxTime) //|| previousSuccessful
	{
		OTTimer = 0;
		OTRunning = false;
		this->OTValIdx = (OTValIdx + 1) % OTVal.size();//change OTval to next value in vector
	}

	if (OTRunning)
	{
		OTTimer += dT;
	}
}

void Smarticle::ControllerMove(int guiState, double torque01, double torque12)
{

	/////////////////////////set which stress actions are actived/////////////////////////
	static const bool LowStressActive = false;
	static const bool MidStressActive = false;
	static const bool OverStressActive = false;
	//////////////////////////////////////////////////////////////////////////////////////

	if (active == false)
	{
		successfulMotion = false;
		return;
	}
	static bool moveDecided = false;
	bool sameMoveType = false;
	prevSuccessful = successfulMotion;
	successfulMotion = false;
	this->prevMoveType = this->moveType;


	///double t1 = abs(torque01);
	///double t2 = abs(torque12);

	//assigns OT and specialState to OT
	ChangeArmColor(torque01, torque12, LowStressActive, MidStressActive, OverStressActive);
	if (OverStressActive)
	{
		CheckOTTimer();
		moveDecided = MoveOverStress(torque01, torque12);
	}
	else//this is here so that if all actions are false, smarticles don't get stuck
	{
		if (!arm0OT || !arm2OT || !OverStressActive)
			specialState = -1;
	}


	if (MidStressActive && !moveDecided)
	{
		moveDecided = MoveMidStress(torque01, torque12);
	}


	if (LowStressActive && !moveDecided)
	{
		CheckLTTimer(torque01, torque12);
		moveDecided = MoveLowStress(torque01, torque12, LTTimer);
	}


	if (specialState != -1)
	{
		AssignState(specialState);
		this->moveType = (MoveType)specialState;
	}
	else
	{
		AssignState(guiState);
		this->moveType = (MoveType)guiState;
	}
	//!(moveType^prevMoveType)
	sameMoveType = (moveType == prevMoveType); // !(xor) gives true if values are equal, false if not
	if (!sameMoveType)
	{
		//GetLog() << "\n&&&&&&&&&&&RESET ERRORRR!&&&&&&&&&&&\n";
		this->armsController->resetCumError = true;
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

std::shared_ptr<ChBody> Smarticle::GetSmarticleBodyPointer()
{
	return arms[1];
}
int	Smarticle::GetID()
{
	return arms[0]->GetPhysicsItem()->GetIdentifier();
}
void Smarticle::setCurrentMoveType(MoveType newMoveType)
{
	this->moveType = newMoveType;
	//global_GUI_value = newMoveType;
}
double Smarticle::GetArmTorque(int index) ///gets torque exerted by controller on arm
{
	if (active)
	{
		double torque = armsController->GetEngine(index)->Get_tor_funct()->Get_y(m_system->GetChTime());
		//return (armsController->ycurr[index] - armsController->yold[index])*torque;
		return torque;
	}
	else
		return 0;
}
double Smarticle::GetTotalTorque()///sums up, for both arms, torque exerted by controller on arm
{
	double totalTorque = 0;
	for (int i = 0; i < numEngs; i++)
	{
		totalTorque += abs(GetArmTorque(i));
	}
	return totalTorque;
}
ChVector<> Smarticle::GetReactTorqueVector(int id)
{
	return getLinkActuator(id)->Get_react_torque();
}
double Smarticle::GetMotTorque(int id)
{
	return getLinkActuator(id)->Get_mot_torque(); //gets same value as one passed into motor from controller
}

void Smarticle::SetBodyFixed(bool mev) {
	for (size_t i = 0; i < numSegs; i++)
	{
		arms[i]->SetBodyFixed(mev);
	}
}
