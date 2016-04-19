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
	m_system->RemoveLink(link_actuator01);
	m_system->RemoveLink(link_actuator12);
	//m_system->RemoveLink(link_revolute01);
	//m_system->RemoveLink(link_revolute12);
	m_system->RemoveBody(arm0);
	m_system->RemoveBody(arm1);
	m_system->RemoveBody(arm2);
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
		double other_angle2){

	smarticleID = sID;
	density = other_density;
	mat_g = surfaceMaterial;
	l = other_l;
	w = other_w;
	r = other_r;
	r2 = other_r2;
	collisionEnvelope = other_envelope;
	initPos = pos;
	rotation = rot;
	angle1 = other_angle;
	angle2 = other_angle2;
	volume = GetVolume();



	mtextureOT->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));
	mtextureArm->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));
	mtextureMid->SetTextureFilename(GetChronoDataFile("cubetexture_blue_bordersBlueOriented.png"));
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
		double other_angle2){

	Properties(sID, other_density, surfaceMaterial, other_envelope, other_l, other_w, other_r, other_r2, pos, rot, other_angle, other_angle2);
	SetDefaultOmega(other_omega);
	//SetOmega(other_omega, true);
}
/////Will added Smarticle properties///////
void Smarticle::Properties(
	int sID,
	int mdumID,
	double other_density,
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
	double other_angHigh){

	Properties(sID, other_density, surfaceMaterial, other_envelope, other_l, other_w, other_r, other_r2, pos, rot, other_angle, other_angle2);
	SetDefaultOmega(other_omega);
	//SetOmega1(other_omega);
	//SetOmega2(other_omega);
	SetOmega(other_omega);
	moveType = (MoveType)global_GUI_value;
	prevMoveType = (MoveType) global_GUI_value;
	//moveTypeIdxs.resize(MoveType::OT, 0);
	moveTypeIdxs.resize(MoveType::OT+1, 0);
	arm0OT = false;
	arm2OT = false;
	torqueThresh2 = other_torThresh2;
	angLow = 0;
	angHigh = 120;
	dumID = mdumID;
	armBroken = false;

	std::tuple<double, double,double,double> a (0.0, 0.0,0.0,0.0);
	//torques = { a, a, a, a, a, a, a }; //probably a better way to do this....
	torques = { a, a }; //probably a better way to do this....
	nextOmega = { 0, 0 };
	nextAngle = { 0,0};
	currTorque = { 0, 0 };
	torqueAvg = std::make_tuple(0, 0, 0, 0);
	initialAng0 = other_angle;
	initialAng1 = other_angle2;
	percentToChangeStressState = .30; //.05
	
	OTThresh = .7;//for single arm 
	MTThresh = OTThresh/3.0;//for both arms?
	LTThresh = .05;//for both arms?

	//initialize OT timer params
	this->OTTimer = 0;
	this->OTMaxTime = .25;
	this->OTRunning = false;
	this->OTVal.emplace_back(GUI1);
	this->OTVal.emplace_back(EXTRA1);
	this->OTVal.emplace_back(EXTRA2);
	this->OTValIdx = genRandInt(0, OTVal.size() - 1);
	GetLog() << "smartRandidx:" << this->OTValIdx << "\n";

	//initialize LT timer params
	this->LTTimer = 0;
	this->LTMaxTime = .25;
	this->LTRunning = false;
	this->LTVal.emplace_back(GUI1);
	this->LTVal.emplace_back(GUI3);
	this->LTValIdx = genRandInt(0, LTVal.size() - 1);

}


void Smarticle::updateTorqueDeque()
{
	std::tuple < double, double, double, double > oldT = torques.back();
	torques.pop_back();
	torques.emplace_front(getLinkActuator(0)->Get_mot_torque(), getLinkActuator(1)->Get_mot_torque(), getLinkActuator(0)->Get_mot_rot_dt(), getLinkActuator(1)->Get_mot_rot_dt());
	updateTorqueAvg(oldT);

}
void Smarticle::updateTorqueAvg(std::tuple <double,double,double,double > oldT)
{ //already assuming it is an avg:
	size_t len =torques.size();
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

	vol = utils::CalcBoxVolume(ChVector<>(len/2.0, r, r2));
	gyr = utils::CalcBoxGyration(ChVector<>(len/2.0, r, r2)).Get_Diag();
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
		switch (armID) {
		case 0:
			arm0_textureAsset = std::make_shared<ChTexture>();
			arm0_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));
			arm->AddAsset(arm0_textureAsset);
			break;
		case 1:
			arm1_textureAsset = std::make_shared<ChTexture>();
			arm1_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_blue_bordersBlueOriented.png"));
			arm->AddAsset(arm1_textureAsset);
			break;
		case 2:
			arm2_textureAsset = std::make_shared<ChTexture>();
			arm2_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));
			arm->AddAsset(arm2_textureAsset);
			break;
		default:
			std::cerr << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
			break;
		}
	}

	arm->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(arm.get(), ChVector<>(len / 2.0, r, r2), ChVector<>(0, 0, 0),QUNIT,visualize);

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
	vol = utils::CalcBoxVolume(ChVector<>(len / 2.0, mr, mr2));
	gyr = utils::CalcBoxGyration(ChVector<>(len / 2.0, mr, mr2)).Get_Diag();
	// create body, set position and rotation, add surface property, and clear/make collision model
	auto arm = std::make_shared<ChBody>();

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
	//mat_g->SetFriction(.05);
	arm->SetMaterialSurface(mat_g);

	double mass = density * vol;
	//double mass = .005;//.043/3.0; //robot weight 43 grams
	arm->GetCollisionModel()->ClearModel();

	if (visualize)
	{
		switch (armID) {
		case 0:
			arm0_textureAsset = std::make_shared<ChTexture>();
			arm0_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));
			arm->AddAsset(arm0_textureAsset);
			break;
		case 1:
		{
			arm1_textureAsset = std::make_shared<ChTexture>();
			arm1_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_SmarticlePicture.png"));
			arm->AddAsset(arm1_textureAsset);

			//auto centerMesh = std::make_shared<ChObjShapeFile>();
			//centerMesh->SetFilename("C:/Users/root/Desktop/Case.obj");
			//centerMesh->SetColor(ChColor(.4, .2, .5, 0));
			//ChMatrix33<> a = centerMesh->Rot;
			////ChQuaternion<> b = a.ClipQuaternion(3,3);
			////b = b*Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(PI/2, PI/2, 0));
			////a.PasteQuaternion(b, 3, 3);
			////GetLog() << a;
			//a.Set_A_Rxyz(ChVector<>(PI / 2, PI / 2, PI / 2));
			//centerMesh->Rot = a;
			//centerMesh->Pos = ChVector<>(1, 1, 1);
			//GetLog() << "\n ROTATED\n";
			//GetLog() << a;
			//GetLog() << centerMesh->Rot;
			//centerMesh->SetFading(0.1);
			//arm->AddAsset(centerMesh);
			//ChFrame<> fr = arm->GetAssetsFrame();
			//fr.SetRot(a);
			break;
		}
		case 2:
			arm2_textureAsset = std::make_shared<ChTexture>();
			arm2_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));
			arm->AddAsset(arm2_textureAsset);
			break;
		default:
			std::cerr << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
			break;
		}
	}
	arm->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(arm.get(), ChVector<>(len / 2.0, mr, mr2), ChVector<>(0, 0, 0), QUNIT, visualize);

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

std::shared_ptr<ChBody> Smarticle::GetArm(int armID) {
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
	return std::make_shared<ChBody>();
}

std::shared_ptr<ChLinkLockRevolute> Smarticle::GetRevoluteJoint(int jointID) {

	switch (jointID) {
	case 0:
		return link_revolute01;
	case 1:
		return link_revolute12;
	default:
		std::cerr << "Error! smarticle can only have joints with ids from {0, 1}" << std::endl;
		break;
	}
	return std::shared_ptr<ChLinkLockRevolute>();
}
void Smarticle::SetSpeed(ChVector<> newSpeed)
{
	arm0->SetPos_dt(newSpeed);
	arm1->SetPos_dt(newSpeed);
	arm2->SetPos_dt(newSpeed);
}
void Smarticle::RotateSmarticle(ChQuaternion<> newRotation)
{
	arm0->SetRot(newRotation);
	arm1->SetRot(newRotation);
	arm2->SetRot(newRotation);
}
//returns (minx,maxx,miny,maxy)
std::vector<double> Smarticle::VertsMinMax()
{
	//GetLog() << this->GetArm(1)->GetPos();
	//armVerts[arm][verts];
	std::vector<double> xyminmax;
	double xmin = armVerts[0][0].x;
	double xmax = armVerts[0][0].x;
	double ymin = armVerts[0][0].y;
	double ymax = armVerts[0][0].y;
	for (int arm = 0; arm < 3; arm++)
	{
		for (int vertex = 0; vertex < 4; vertex++)
		{
			xmin = std::min(xmin, armVerts[arm][vertex].x);
			xmax = std::max(xmax, armVerts[arm][vertex].x);

			ymin = std::min(ymin, armVerts[arm][vertex].y);
			ymax = std::max(ymax, armVerts[arm][vertex].y);
		}
	}
	//.0065;
	//(2, l, armt, armt2, ChVector<>((w / 2.0 - (jointClearance)+cos(-angle2)*l / 2), 0, -(l / 2.0)*sin(-angle2) - offPlaneoffset), quat2);
	//GetLog() << "\n" << xmin << " " << xmax << " " << ymin << " " << ymax;

	xyminmax.emplace_back(xmin);
	xyminmax.emplace_back(xmax);
	xyminmax.emplace_back(ymin);
	xyminmax.emplace_back(ymax);
	//xyminmax.emplace_back((std::min(std::min(std::min(armVerts[2][3].x, armVerts[2][1].x), armVerts[0][3].x), armVerts[0][1].x)));
	//xyminmax.emplace_back((std::max(std::max(std::max(armVerts[2][3].x, armVerts[2][1].x), armVerts[0][3].x), armVerts[0][1].x)));

	//xyminmax.emplace_back((std::min(std::min(std::min(armVerts[2][3].y, armVerts[2][1].y), armVerts[0][3].y), armVerts[0][1].y)));
	//xyminmax.emplace_back((std::max(std::max(std::max(armVerts[2][3].y, armVerts[2][1].y), armVerts[0][3].y), armVerts[0][1].y)));

	return xyminmax;
}
double *Smarticle::Project(double minmax[]){
	
	return minmax;
}
void Smarticle::UpdateState()
{
	GetArm(0)->Update();
	GetArm(1)->Update();
	GetArm(2)->Update();
}
void Smarticle::SetEdges()
{
	UpdateState();
	//ChVector<>pos(mySmarticlesVec[0]->GetArm(1)->TransformPointLocalToParent(ChVector<>(w_smarticle / 2.0, 0, 0)));
	//ChVector<>pos2(mySmarticlesVec[0]->GetArm(2)->TransformPointLocalToParent(ChVector<>(l_smarticle / 2, 0, 0)));

	//vector3df armpos(pos.x, pos.y, pos.z);
	//vector3df armpos2(pos2.x, pos2.y, pos2.z);
	//application.GetVideoDriver()->setTransform(irr::video::ETS_WORLD, core::IdentityMatrix);
	//application.GetVideoDriver()->draw3DLine(armpos - vector3df(0, 10, 0),
	//	armpos + vector3df(0, 10, 0), irr::video::SColor(70, 255, 0, 0));

	//CreateArm2(0, l, armt, armt2, ChVector<>((-w / 2.0 + (jointClearance)-cos(-angle1)*l / 2), 0, -(l / 2.0)*sin(-angle1) - offPlaneoffset), quat0);
	//CreateArm2(1, w, r, r2, ChVector<>(0, 0, 0));
	//CreateArm2(2, l, armt, armt2, ChVector<>((w / 2.0 - (jointClearance)+cos(-angle2)*l / 2), 0, -(l / 2.0)*sin(-angle2) - offPlaneoffset), quat2);
	

	double armt = r;
	double armt2 = .00806 / 2 * sizeScale; //8.06 mm with solar 3.2 without
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

		arm0Front = ChVector<>(-arm0Front.y, arm0Front.x, arm0Front.z);
		arm0OuterEdge = ChVector<>(-arm0OuterEdge.y, arm0OuterEdge.x, arm0OuterEdge.z);
		arm0Back = ChVector<>(-arm0Back.y, arm0Back.x, arm0Back.z);
		arm0Edge = ChVector<>(-arm0Edge.y, arm0Edge.x, arm0Edge.z);

		//arm1
		arm1Front =(armVerts[1][0]) - (armVerts[1][3]);
		arm10Shared = (armVerts[1][1]) - (armVerts[1][0]);
		arm1Back = (armVerts[1][2]) - (armVerts[1][1]);
		arm12Shared = (armVerts[1][3]) - (armVerts[1][2]);

		arm1Front = ChVector<>(-arm1Front.y, arm1Front.x, arm1Front.z);
		arm10Shared = ChVector<>(-arm10Shared.y, arm10Shared.x, arm10Shared.z);
		arm1Back = ChVector<>(-arm1Back.y, arm1Back.x, arm1Back.z);
		arm12Shared = ChVector<>(-arm12Shared.y, arm12Shared.x, arm12Shared.z);



		//arm2
		arm2Front =			(armVerts[2][0]) - (armVerts[2][3]);
		arm2Edge =			(armVerts[2][1]) - (armVerts[2][0]);
		arm2Back=				(armVerts[2][2]) - (armVerts[2][1]);
		arm2OuterEdge = (armVerts[2][3]) - (armVerts[2][2]);

		arm2Front = ChVector<>(-arm2Front.y, arm2Front.x, arm2Front.z);
		arm2Edge = ChVector<>(-arm2Edge.y, arm2Edge.x, arm2Edge.z);
		arm2Back = ChVector<>(-arm2Back.y, arm2Back.x, arm2Back.z);
		arm2OuterEdge = ChVector<>(-arm2OuterEdge.y, arm2OuterEdge.x, arm2OuterEdge.z);


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
void Smarticle::RotateSmarticleBy(ChQuaternion<> newRotation)
{


	GetLog() << "arm1:" << arm1->GetPos();
	GetLog() << "arm2:" << arm2->GetPos();
	double l_mod = l + 2 * r2 - jointClearance;
	GetLog() << "\nl-mod" << l_mod << " w:" << w << " l:"<<l<<"\n";
	GetLog() << "local to parent" << arm1->TransformPointLocalToParent(ChVector<>(w / 2.0 - (jointClearance)+cos(-angle2)*l / 2, 0, -(l / 2.0)*sin(-angle2) - offPlaneoffset));
	GetLog() << "local to parent" << arm1->TransformPointLocalToParent(ChVector<>(w / 2.0+l, 0,0));
	GetLog() << "local to parent" << arm2->TransformPointLocalToParent(ChVector<>(l, 0, 0));
}
void Smarticle::TransportSmarticle(ChVector<> newPosition)
{
	arm0->SetPos(arm0->GetPos() - arm1->GetPos() + newPosition);
	arm2->SetPos(arm2->GetPos() - arm1->GetPos() + newPosition);
	arm1->SetPos(newPosition);

}
void Smarticle::CreateJoints() {
	// link 1
	link_revolute01 = std::make_shared<ChLinkLockRevolute>();
	link_revolute12 = std::make_shared<ChLinkLockRevolute>();
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

	link_actuator01 = std::make_shared<ChLinkEngine>();
	link_actuator12 = std::make_shared<ChLinkEngine>();
	
	ChVector<> pR01 (-w / 2.0, 0, 0);
	ChVector<> pR12 (w / 2.0, 0, 0);
	if (!stapleSize)
	{
		pR01 = ChVector<>(-(w / 2 - jointClearance), 0, -offPlaneoffset);
		pR12 = ChVector<>((w / 2 - jointClearance), 0, -offPlaneoffset);
	}

	ChQuaternion<> qx1 = Q_from_AngAxis(-PI_2, VECT_X);
	ChQuaternion<> qx2 = Q_from_AngAxis(PI_2, VECT_X);
	ChQuaternion<> qy1 = Q_from_AngAxis(GetAngle1(), VECT_Z);
	ChQuaternion<> qy2 = Q_from_AngAxis(GetAngle2(), VECT_Z);

	qx1.Normalize();
	qx2.Normalize();
	qy1.Normalize();
	qy2.Normalize();


	link_actuator01->Initialize(arm0, arm1, false, ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx1*qy1), ChCoordsys<>(rotation.Rotate(pR01) + initPos, rotation*qx1));
	m_system->AddLink(link_actuator01);
	link_actuator01->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);

	link_actuator12->Initialize(arm2, arm1, false, ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx2*qy2), ChCoordsys<>(rotation.Rotate(pR12) + initPos, rotation*qx2));
	m_system->AddLink(link_actuator12);
	link_actuator12->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);

	//auto mfun0 = std::dynamic_pointer_cast<ChFunction_Const>(link_actuator01->Get_tor_funct());
	//mfun0->Set_yconst(0);
	//auto mfun1= std::dynamic_pointer_cast<ChFunction_Const>(link_actuator12->Get_tor_funct());
	//mfun1->Set_yconst(0);
}

void Smarticle::Create() {
	jointClearance = 0;
	double l_mod;
	//double l_mod = l + 2 * r2 - jointClearance;

	//ChQuaternion<> quat0 = Q_from_AngAxis(angle1, VECT_Y);
	//ChQuaternion<> quat2 = Q_from_AngAxis(-angle2, VECT_Y);
	ChQuaternion<> quat0 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, -angle1, 0));
	ChQuaternion<> quat2 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, angle2, 0));
	quat0.Normalize();
	quat2.Normalize();


	if (stapleSize)
	{
		offPlaneoffset = 0;
		l_mod = l + 2 * r2 - jointClearance;
		CreateArm(0, l_mod, ChVector<>(-w / 2.0 - (l / 2.0)*cos(angle1), 0, -(l_mod / 2.0 - r2)*sin(angle1)), quat0);
		CreateArm(1, w, ChVector<>(0, 0, 0));
		CreateArm(2, l_mod, ChVector<>(w / 2.0 +  (l / 2.0)*cos(angle2), 0, -(l_mod / 2.0 - r2)*sin(angle2)), quat2);
	}
	else
	{
		offPlaneoffset = .00825-r2;//smarticle arms arent centered in plane, arm is offset 8.25mm from front or offPlaneoffset-t2
		jointClearance = .0065;//6.5 mm in x dir jc in y dir is r2
		double armt = r;
		double armt2 = .00806 / 2 * sizeScale; //8.06 mm with solar 3.2 without

		CreateArm2(0, l, armt, armt2, ChVector<>((-w / 2.0 + (jointClearance)-cos(-angle1)*l / 2), 0, -(l / 2.0)*sin(-angle1) - offPlaneoffset), quat0);
		CreateArm2(1, w, r, r2, ChVector<>(0, 0, 0));
		CreateArm2(2, l, armt, armt2, ChVector<>((w / 2.0 - (jointClearance)+cos(-angle2)*l / 2), 0, -(l / 2.0)*sin(-angle2) - offPlaneoffset), quat2);
	}

	CreateActuators();

	// mass property
	mass = arm0->GetMass() + arm1->GetMass() + arm2->GetMass();


	if (genRand() < pctActive)
	{
		active = true;
		armsController = (new Controller(m_system,this));
		armsController->outputLimit= torqueThresh2;
		armsController->omegaLimit = omegaLim;
	}
	else
	{
		active = false;
		//link_actuator01->SetDisabled(true);
		//link_actuator12->SetDisabled(true);
	}
}

std::shared_ptr<ChFunction> Smarticle::GetActuatorFunction(int actuatorID) {
	if (actuatorID == 0) {
		return function01;
	} else if (actuatorID == 1) {
		return function12;
	} else {
		std::cout << "Error! smarticle can only have actuators with ids from {0, 1}" << std::endl;
	}
	return std::shared_ptr<ChFunction>(NULL);
}

void Smarticle::SetActuatorFunction(int actuatorID, std::shared_ptr<ChFunction> actuatorFunction) {
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
	std::shared_ptr<ChLinkEngine> mlink_actuator;
	if (actuatorID == 0) {
		mlink_actuator = link_actuator01;
	} else {
		mlink_actuator = link_actuator12;
	}

//	auto mfun1 = std::dynamic_pointer_cast<ChFunction_Const>(mlink_actuator->Get_rot_funct());
//	mfun1->Set_yconst(diffTheta + mfun1->Get_yconst());
	//auto mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(mlink_actuator->Get_spe_funct());
	auto mfun1 = std::dynamic_pointer_cast<ChFunction_Const>(mlink_actuator->Get_tor_funct());

	
	//mfun2->Set_yconst(omega);
	mfun1->Set_yconst(omega);
}

void Smarticle::SetActuatorFunction(int actuatorID, double omega) {
	std::shared_ptr<ChLinkEngine> mlink_actuator;
	if (actuatorID == 0) {
		mlink_actuator = link_actuator01;
	} else {
		mlink_actuator = link_actuator12;
	}

	//auto mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(mlink_actuator->Get_spe_funct());
	auto mfun1 = std::dynamic_pointer_cast<ChFunction_Const>(mlink_actuator->Get_tor_funct());
	//mfun2->Set_yconst(omega);
	mfun1->Set_yconst(omega);
}

std::shared_ptr<ChLinkEngine> Smarticle::getLinkActuator(int id)
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
double Smarticle::GetInitialAngle(int id)
{
	if (id == 0)
		return initialAng0;
	else
		return initialAng1;
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
	OTThresh = OTThresh*torqueThresh2;
	MTThresh = MTThresh*torqueThresh2;
	LTThresh = LTThresh*torqueThresh2;
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
	smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
	angPair.first = angVals.x;
	angPair.second = angVals.y;
	firstAngPair.first = angVals.x;
	firstAngPair.second = angVals.y;
	//TODO need to rewrite the below way of reading file, very ugly!

	ot.emplace_back(GetAngle1(),GetAngle2());
	ot.emplace_back(GetAngle1(), GetAngle2());
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
	if (extra1.size() < 1)
	{
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
			angPair.first = angVals.x;
			angPair.second = angVals.y;

			extra1.push_back(angPair);
			//GetLog() << angVals.x << " " << angVals.y << " ddch:" << ddCh << "\n";
			//exit(-1);
			if (ddCh == '#')
				break;
		}
	}

	if (extra2.size() < 1)
	{
		while (smarticleMoves.good()) {
			smarticleMoves >> angVals.x >> ddCh >> angVals.y >> ddCh;
			angPair.first = angVals.x;
			angPair.second = angVals.y;

			extra2.push_back(angPair);
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

	if (t0+t1 > MTThresh)
	{// MT<t0<OT //TODO perhaps make this function if(t0+t1>2*MT) since servo can only sense stress from both
		if (genRand() < percentToChangeStressState)
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
		specialState = OTVal.at(OTValIdx);
		return true;
	}
	
	specialState = -1;
	return false;
}

//returns true if either arm is OT
bool Smarticle::ChangeArmColor(double torque01, double torque12)
{
	bool LTMTcolor = false;
	//for vibration upon OT, change degreesToVibrate to amount you wish to vibrate
	double degreesToVibrate = 0;
	double moveAmt = degreesToVibrate*D2R;
	

	if (abs(torque01) > OTThresh)
	{
		//this->setCurrentMoveType(OT);
		//mv = &ot;
		if (!arm0OT)//if not previously OT
		{
			arm0OT = true;
			arm0_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));

			this->ot.clear();
			//this->ot.emplace_back(GetAngle1() + sign(torque01)*moveAmt, GetAngle2() + sign(torque12)*moveAmt);
			this->ot.emplace_back(GetAngle1() + moveAmt, GetAngle2() + moveAmt);
			this->ot.emplace_back(GetAngle1(), GetAngle2());
			this->ot.emplace_back(GetAngle1() - moveAmt, GetAngle2() - moveAmt);
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
			arm0_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));

		}
		// nothing needs to be done if not prev OT
	}

	/////////////////////ARM2///////////////////////
	if (abs(torque12) > OTThresh)
	{
		//this->setCurrentMoveType(OT);
		//mv = &ot;
		if (!arm2OT)//if not previously OT
		{
			arm2OT = true;
			arm2_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));

			this->ot.clear();
			this->ot.emplace_back(GetAngle1() + moveAmt, GetAngle2() + moveAmt);
			this->ot.emplace_back(GetAngle1(), GetAngle2());
			this->ot.emplace_back(GetAngle1() - moveAmt, GetAngle2() - moveAmt);
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
			arm2_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_Smart_bordersBlack.png"));

		}
		// nothing needs to be done if not prev OT
	}
	
	if (LTMTcolor)
	{
		if (abs(torque01) + abs(torque12) > MTThresh && !arm2OT && !arm0OT)
		{
			arm0_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_orange_borderOrange.png"));
			arm2_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_orange_borderOrange.png"));
		}
		if (abs(torque01) + abs(torque12) < LTThresh)
		{
			arm0_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_green_borderGreen.png"));
			arm2_textureAsset->SetTextureFilename(GetChronoDataFile("cubetexture_green_borderGreen.png"));
		}

	}
		return false;

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
	if ((abs(t1)+abs(t2))<LTThresh)
	{
		if (genRand() < percentToChangeStressState)
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
		if (genRand() < percentToChangeStressState)
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
	static const bool OverStressActive = true;
	//////////////////////////////////////////////////////////////////////////////////////

	if (active == false)
	{
		successfulMotion = false;
		return;
	}
	bool moveDecided = false;
	bool sameMoveType = false;
	prevSuccessful = successfulMotion;
	successfulMotion = false;
	this->prevMoveType = this->moveType;


	///double t1 = abs(torque01);
	///double t2 = abs(torque12);

	//assigns OT and specialState to OT
	ChangeArmColor(torque01, torque12);
	
	if (OverStressActive && !moveDecided)
	{
		CheckOTTimer();
		moveDecided = MoveOverStress(torque01, torque12);
	}
	else//this is here so that if all actions are false, smarticles don't get stuck
	{
		if (!arm0OT || !arm2OT || !OverStressActive)
			specialState = -1;
	}
	

	if (MidStressActive  && !moveDecided)
	{
		moveDecided = MoveMidStress(torque01, torque12);
	}


	if (LowStressActive  && !moveDecided)
	{
		CheckLTTimer(torque01,torque12);
		moveDecided = MoveLowStress(torque01, torque12, LTTimer);
	}


	if (specialState != -1)
	{
		AssignState(specialState);
		this->moveType = (MoveType) specialState;
	}
	else
	{
		AssignState(guiState);
		this->moveType = (MoveType)guiState;
	}
	//!(moveType^prevMoveType)
	sameMoveType = (moveType==prevMoveType); // !(xor) gives true if values are equal, false if not
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
double Smarticle::GetReactTorqueLen(int index)
{
	if (index == 0)
		return GetReactTorqueLen01();
	else
		return GetReactTorqueLen12();
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
