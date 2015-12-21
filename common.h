#ifndef _common_h
#define _common_h

#include <stdlib.h>  // system, rand, srand, RAND_MAX
#include "utils/ChUtilsGeometry.h"
#include "utils/ChUtilsCreators.h"
#include "core/ChVector.h"
#include <vector>
#include <chrono>

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#define USE_PARALLEL false
#define irrlichtVisualization true
#define stapleSize false
	//extern double sizeScale;
	//extern double dT;
	//extern bool bucket_exist;
	//extern double vibAmp;//vibrate by some amount of degrees back and forth

	////////////

	extern double sizeScale;
	extern double dT;
	extern bool bucket_exist;



	enum SmarticleType { SMART_ARMS, SMART_U };
	enum BucketType { KNOBCYLINDER, HOOKRAISE, STRESSSTICK, CYLINDER, BOX, HULL, RAMP, HOPPER, DRUM };
	extern SmarticleType smarticleType;
	extern BucketType bucketType;


	////////////


	extern double bucket_rad;
	extern double vibAmp;
	extern double rampAngle;
	extern double rampInc;
	extern double drum_omega;
	extern double drum_freq;
	extern double inc;
	extern double p_gain;
	extern double i_gain;
	extern double d_gain;
	extern double fric;
	extern chrono::ChVector<> bucket_interior_halfDim;
	
	
	//extern SmarticleType smarticleType;
	//extern BucketType bucketType;
	extern chrono::ChSharedPtr<chrono::ChBody> bucket;
	extern chrono::ChSharedPtr<chrono::ChBody> bucket_bott;
	extern std::vector<chrono::ChSharedPtr<chrono::ChBody>> bucket_bod_vec;


#endif



////////////deprecated code which may still be useful in the future////////////////////////////
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//bool Smarticle::MoveToAngle2(std::vector<std::pair<double, double>> *v, double momega1, double momega2, MoveType mtype)
//{
//	if (!active)
//		return false;
//	if (arm0OT || arm2OT) //turn off motion at limit
//	{
//		this->SetActuatorFunction(0, 0);
//		this->SetActuatorFunction(1, 0);
//		return false;
//	}
//
//	//real ang01 and ang12
//	double ang01 = link_actuator01->Get_mot_rot();
//	double ang12 = link_actuator12->Get_mot_rot();
//
//	SetAngles(ang01, ang12);
//
//	//next ang01 and ang12
//	double nextAng01 = v->at((moveTypeIdxs.at(moveType) + 1) % v->size()).first;
//	double nextAng12 = v->at((moveTypeIdxs.at(moveType) + 1) % v->size()).second;
//
//
//	double expAng01 = v->at(moveTypeIdxs.at(mtype)).first;
//	double expAng12 = v->at(moveTypeIdxs.at(mtype)).second;
//
//
//	//GetLog() << "(ang1,ang2)=(" << ang01 * 180 / CH_C_PI << "," << ang12 * 180 / CH_C_PI << ") (nextang1,nextang2)=(" << nextAng01 << ","<< nextAng12<<")\n";
//	//exit(-1);
//	//different real - expected
//	//double ang01Diff = ang01-expAng01;
//	//double ang12Diff = ang12-expAng12;
//
//	double omega01 = ChooseOmegaAmount(momega1, ang01, expAng01);
//	double omega12 = ChooseOmegaAmount(momega2, ang12, expAng12);
//
//
//	//if within some threshold distance to where curr angle is supposed to be:
//	if (omega01 == 0)
//	{
//		if (omega12 == 0)
//		{
//			//arm01=right, arm12=right both angles should try to move to the next angle
//			this->SetActuatorFunction(0, ChooseOmegaAmount(momega1, ang01, nextAng01));
//			this->SetActuatorFunction(1, ChooseOmegaAmount(momega2, ang12, nextAng12));
//			return true;//was able to successfully move to next index
//		}
//		else
//		{
//			//arm01=right, arm12=wrong, arm12 must catch up
//			this->SetActuatorFunction(0, 0);
//			this->SetActuatorFunction(1, omega12);
//			return false;
//		}
//	}
//	// from this point we know arm01 is wrong
//	if (omega12 == 0)
//	{
//		//arm01=wrong, arm12=right
//		this->SetActuatorFunction(0, omega01);
//		this->SetActuatorFunction(1, 0);
//		return false;
//	}
//	else
//	{
//		//arm01=wrong, arm12=wrong
//		this->SetActuatorFunction(0, omega01);
//		this->SetActuatorFunction(1, omega12);
//		return false;
//	}
//}
//
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//void Smarticle::MoveLoop2(int guiState, double torque01, double torque12)
//{
//	//initialize boolean describing same moveType as last step
//	bool sameMoveType = false;
//	//static bool successfulMotion = false;
//	bool prevSucessful = successfulMotion;
//	successfulMotion = false;
//	//this pointer will point to the correct moveType vector
//	std::vector<std::pair<double, double>> *v;
//	//set prevMoveType to previous value
//	this->prevMoveType = this->moveType;
//	//get previous values from last timestep
//	double ang01 = link_actuator01->Get_mot_rot();
//	double ang12 = link_actuator12->Get_mot_rot();

//	double omega1Prev = link_actuator01->Get_mot_rot_dt();
//	double omega2Prev = link_actuator12->Get_mot_rot_dt();



//	//GetLog() << "\n*********" << torque01 << " " << torque12 << " thresh: " << torqueThresh2;


//	//determine moveType
//	switch (guiState)
//	{
//	case 0:
//		this->setCurrentMoveType(GLOBAL);
//		v = &global;
//		break;
//	case 1:
//		this->setCurrentMoveType(GUI1);
//		v = &gui1;
//		break;
//	case 2:
//		this->setCurrentMoveType(GUI2);
//		v = &gui2;
//		break;
//	case 3:
//		this->setCurrentMoveType(GUI3);
//		v = &gui3;
//		break;
//	case 4:
//		this->setCurrentMoveType(VIB);
//		v = &vib;
//		break;
//	case 5:
//		this->setCurrentMoveType(VIB);
//		v = &vib;
//		break;
//	default:
//		this->setCurrentMoveType(GLOBAL);
//		v = &global;
//		break;
//	}




//	static ChVector<> rel01 = link_actuator01->GetRelAxis();
//	static ChVector<> rel12 = link_actuator12->GetRelAxis();
//	static double relr01 = link_actuator01->GetDist();
//	static double relr12 = link_actuator12->GetDist();
//	//GetLog() << "\n" << "rel1:" << rel01 << "rel2:" << rel12 << "\n";


//	//if (fabs(rel01.y) > .05 || fabs(rel12.y) > .05)//if angle in x or y is > .1 radians, it is definitely broken
//	//{
//	//	GetLog() << "angle bad break! \n";
//	//	armBroken = true;
//	//}
//	//if (fabs(relr01) > .075*r2 || fabs(relr12) > .075*r2)//if distance between markers is .025% of thickness, break!
//	//{
//	//	GetLog() << "distance wrong break! \n";
//	//	armBroken = true;
//	//}

//	//ChVector<> rel01 = link_actuator01->GetRelRotaxis();
//	//ChVector<> rel12 = link_actuator12->GetRelRotaxis();
//	////GetLog() <<"\n"<< rel01 << "\n";
//	//if (fabs(rel01.x) > .1 || fabs(rel01.y) > .1
//	//	|| fabs(rel12.x) > .1 || fabs(rel12.x) > .1)//if angle in x or y is > .1 radians, it is definitely broken
//	//{
//	//	armBroken = true;
//	//}

//	/////////////torque color change was here/////////////
//	if (torque01 > torqueThresh2 || torque12 > torqueThresh2)//one arm is OT 
//	{
//		this->setCurrentMoveType(OT);
//		v = &ot;
//	}
//	ChangeArmColor(torque01, torque12);

//	//////////////////////////////////////////////////////


//	sameMoveType = !(moveType^prevMoveType); // !(xor) gives true if values are equal, false if not

//	switch (this->moveType) //have this in case I want to add different action based on move type
//	{
//	case GLOBAL://TODO implement different case if sameMoveType was wrong

//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		break;
//	case OT:
//		//if (sameMoveType){}
//		successfulMotion = MoveToAngle2(v, 0, 0, moveType);
//		break;
//	case GUI1:

//		//if (sameMoveType)
//		//{
//		//	if (omega1Prev == 0 && omega2Prev == 0)
//		//	{

//		//		successfulMotion = true;
//		//		break;
//		//	}
//		//}
//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		break;
//	case GUI2:

//		//if (sameMoveType){}
//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		break;
//	case GUI3:

//		//if (sameMoveType){}
//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		break;
//	case VIB:
//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		//GetLog() << "(0,1,2):" << v->at(0).first << "," << v->at(1).first << "," << v->at(2).first;
//		//exit(-1);
//		break;
//	}
//	//add 1 to size if move was successful (i.e. can move on to next move index if reached previous one)
//	if (successfulMotion&&active && (!arm0OT&&!arm2OT))
//	{
//		moveTypeIdxs.at(moveType) = ((moveTypeIdxs.at(moveType) + 1) % v->size());
//	}

//	return;
//}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//void Smarticle::MoveLoop2(int guiState = 0)
//{
//	//initialize boolean describing same moveType as last step
//	bool sameMoveType = false;
//	//static bool successfulMotion = false;
//	bool prevSucessful = successfulMotion;
//	successfulMotion = false;
//	//this pointer will point to the correct moveType vector
//	std::vector<std::pair<double, double>> *v;
//	//set prevMoveType to previous value
//	this->prevMoveType = this->moveType;
//	//get previous values from last timestep
//	double ang01 = link_actuator01->Get_mot_rot();
//	double ang12 = link_actuator12->Get_mot_rot();

//	double omega1Prev = link_actuator01->Get_mot_rot_dt();
//	double omega2Prev = link_actuator12->Get_mot_rot_dt();

//	double torque01 = fabs(link_actuator01->Get_react_torque().Length2()); //use length2 to avoid squareroot calculation be aware of blowing up because too high torque overflows double
//	double torque12 = fabs(link_actuator12->Get_react_torque().Length2());
//	//GetLog() << "\n*********" << torque01 << " " << torque12 << " thresh: " << torqueThresh2;


//	//determine moveType
//	switch (guiState)
//	{
//	case 0:
//		this->setCurrentMoveType(GLOBAL);
//		v = &global;
//		break;
//	case 1:
//		this->setCurrentMoveType(GUI1);
//		v = &gui1;
//		break;
//	case 2:
//		this->setCurrentMoveType(GUI2);
//		v = &gui2;
//		break;
//	case 3:
//		this->setCurrentMoveType(GUI3);
//		v = &gui3;
//		break;
//	case 4:
//		this->setCurrentMoveType(VIB);
//		v = &vib;
//		break;
//	case 5:
//		this->setCurrentMoveType(VIB);
//		v = &vib;
//		break;
//	default:
//		this->setCurrentMoveType(GLOBAL);
//		v = &global;
//		break;
//	}




//	static ChVector<> rel01 = link_actuator01->GetRelAxis();
//	static ChVector<> rel12 = link_actuator12->GetRelAxis();
//	static double relr01 = link_actuator01->GetDist();
//	static double relr12 = link_actuator12->GetDist();
//	//GetLog() << "\n" << "rel1:" << rel01 << "rel2:" << rel12 << "\n";


//	//if (fabs(rel01.y) > .05 || fabs(rel12.y) > .05)//if angle in x or y is > .1 radians, it is definitely broken
//	//{
//	//	GetLog() << "angle bad break! \n";
//	//	armBroken = true;
//	//}
//	//if (fabs(relr01) > .075*r2 || fabs(relr12) > .075*r2)//if distance between markers is .025% of thickness, break!
//	//{
//	//	GetLog() << "distance wrong break! \n";
//	//	armBroken = true;
//	//}

//	//ChVector<> rel01 = link_actuator01->GetRelRotaxis();
//	//ChVector<> rel12 = link_actuator12->GetRelRotaxis();
//	////GetLog() <<"\n"<< rel01 << "\n";
//	//if (fabs(rel01.x) > .1 || fabs(rel01.y) > .1
//	//	|| fabs(rel12.x) > .1 || fabs(rel12.x) > .1)//if angle in x or y is > .1 radians, it is definitely broken
//	//{
//	//	armBroken = true;
//	//}

//	if (torque01 > torqueThresh2 || torque12 > torqueThresh2)//one arm is OT 
//	{
//		this->setCurrentMoveType(OT);
//		v = &ot;

//		if (torque01 > torqueThresh2) // arm 0 is overtorqued
//		{
//			if (!arm0OT)//if arm was previously not OT, add color asset
//			{
//				arm0OT = true;
//				arm0->AddAsset(mtextureOT);
//				this->ot.clear();
//				this->ot.emplace_back(GetAngle1(), GetAngle2());
//			}
//		}
//		if (torque12 > torqueThresh2)// arm 2 is overtorqued
//		{
//			if (!arm2OT)
//			{
//				arm2OT = true;
//				arm2->AddAsset(mtextureOT);
//				this->ot.clear();
//				this->ot.emplace_back(GetAngle1(), GetAngle2());
//			}
//		}

//	}
//	else //arms are not OT 
//	{
//		if (arm0OT)//if arm was previously OT, pop off previous armOT color asset
//		{
//			arm0OT = false;
//			arm0->GetAssets().pop_back();
//		}
//		if (arm2OT)
//		{
//			arm2OT = false;
//			arm2->GetAssets().pop_back();
//		}

//	}

//	sameMoveType = !(moveType^prevMoveType); // !(xor) gives true if values are equal, false if not

//	switch (this->moveType) //have this in case I want to add different action based on move type
//	{
//	case GLOBAL://TODO implement different case if sameMoveType was wrong

//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		break;
//	case OT:
//		//if (sameMoveType){}
//		successfulMotion = MoveToAngle2(v, 0, 0, moveType);
//		break;
//	case GUI1:

//		if (sameMoveType)
//		{
//			if (omega1Prev == 0 && omega2Prev == 0)
//			{

//				successfulMotion = true;
//				break;
//			}
//		}
//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		break;
//	case GUI2:

//		//if (sameMoveType){}
//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		break;
//	case GUI3:

//		//if (sameMoveType){}
//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		break;
//	case VIB:
//		successfulMotion = MoveToAngle2(v, omega1, omega2, moveType);
//		//GetLog() << "(0,1,2):" << v->at(0).first << "," << v->at(1).first << "," << v->at(2).first;
//		//exit(-1);
//		break;
//	}
//	//add 1 to size if move was successful (i.e. can move on to next move index if reached previous one)
//	if (successfulMotion&&active && (!arm0OT&&!arm2OT))
//	{
//		moveTypeIdxs.at(moveType) = ((moveTypeIdxs.at(moveType) + 1) % v->size());
//	}

//	return;
//}


///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//this code was so that I could run through smarticle vector only once
//true if removed, false if not fixed
//bool FixSmarticles(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mySmarticlesVec, Smarticle &sPtr, int tstep, int idx) //TODO remake method
//{
//	//recycle smarticles if bucket or hopper
//	if (bucketType == HOPPER && bucket_exist == false) //if hopper, put smarticles back inside after reaching below hopper if bucket_bott still exists delete
//	{
//		recycleSmarticles(mphysicalSystem, sPtr);
//		return false;
//	}
//	if (sPtr.armBroken)
//	{
//		sPtr.~Smarticle();
//		mySmarticlesVec.erase(mySmarticlesVec.begin() + idx);
//		GetLog() << "\narm broken removing smarticle\n";
//		return true;
//	}
//	if (bucketType != HOPPER)
//	{
//		if (sPtr.GetArm(1)->GetPos().z < -3.0*bucket_interior_halfDim.z) //if far below bucket
//		{
//			sPtr.~Smarticle();
//			mySmarticlesVec.erase(mySmarticlesVec.begin() + idx);
//			GetLog() << "\nbelow bucket\n";
//			return true;
//		}
//		if ((bucketType == CYLINDER || bucketType == STRESSSTICK || bucketType == HOOKRAISE || bucketType == KNOBCYLINDER) && !IsInRadial(sPtr.Get_cm(), bucket_bott->GetPos() + ChVector<>(0, 0, bucket_interior_halfDim.z), ChVector<>(bucket_rad, bucket_bott->GetPos().z, bucket_bott->GetPos().z + 3 * bucket_interior_halfDim.z))) //if outside radius
//		{
//			sPtr.~Smarticle();
//			//mySmarticlesVec.erase(mySmarticlesVec.begin() + idx);
//			GetLog() << "\noutside radius removing bucketType==CYLINDER || STRESSSTICK!\n";
//			return true;
//		}
//	}
//	else //when does this happen?
//	{
//
//		if (!IsInRadial(sPtr.Get_cm(), bucket_bott->GetPos(), ChVector<>(2 * bucket_rad, -4.0*bucket_interior_halfDim.z, 4.0*bucket_interior_halfDim.z)))
//		{
//			sPtr.~Smarticle();
//			mySmarticlesVec.erase(mySmarticlesVec.begin() + idx);
//			GetLog() << "\noutside radius else\n";
//			return true;
//		}
//	}
//	return false;
//
//}

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//double Find_Max_Z(CH_SYSTEM& mphysicalSystem) {
//	std::string smarticleTypeName;
//	if (smarticleType == SMART_ARMS) {
//		smarticleTypeName = "smarticle_arm";
//	}
//	else if (smarticleType == SMART_U) {
//		smarticleTypeName = "smarticle_u";
//	}
//	else {
//		std::cout << "Error! Smarticle type is not set correctly" << std::endl;
//	}
//	double zMax = -999999999;
//
//
//
//	//std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
//	std::vector<ChSharedPtr<ChBody> >::iterator ibody = mphysicalSystem.Get_bodylist()->begin();
//	for (size_t i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
//		//ChBody* bodyPtr = *(myIter + i);
//		ChSharedPtr<ChBody> bodyPtr = *(ibody + i);
//		if (strcmp(bodyPtr->GetName(), smarticleTypeName.c_str()) == 0) {
//			if (zMax < bodyPtr->GetPos().z) {
//				//zMax = bodyPtr->GetPos().z;
//				zMax = bodyPtr->GetPos().z - bucket_bott->GetPos().z;
//			}
//		}
//	}
//	return zMax;
//}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//void recycleSmarticles(CH_SYSTEM& mphysicalSystem, Smarticle& sPtr)
//{
//	double pos = -.75*bucket_interior_halfDim.z;//z position below which smarticles are regenerated above pile inside container
//	double ang = 2 * CH_C_PI / 5;
//	double rp = MyRand()*ang / 4; //add slight offset to angInc to allow particles not always fall in nearly same position
//	static int recycledSmarticles = 0;
//	static int inc = 0;
//	//ChVector<> myPos = bucket_ctr + ChVector<>(sin(ang * i + phase) *(bucket_rad / 2 + w*MyRand()), //TODO for hopper no -w/2.0
//	//	cos(ang*i + phase)*(bucket_rad / 2 + w*MyRand() - w / 2.0),
//	//	zpos);
//	if (sPtr.GetArm(1)->GetPos().z < pos)
//	{
//		if (bucketType == HOPPER)
//		{
//			sPtr.TransportSmarticle(bucket_ctr + ChVector<>(
//				sin(ang*inc + rp)*(bucket_rad / 2 + 4 * w_smarticle*(MyRand() - 1 / 2.0)),
//				cos(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(MyRand() - 1 / 2.0)),
//				bucket_interior_halfDim.z * 2
//				));
//
//			//sPtr->SetSpeed(sPtr->GetArm(1)->GetPos_dt() / 4);
//			sPtr.SetSpeed(ChVector<>(0, 0, -9.8*.01 / 2.0 - w_smarticle / .01));
//		}
//		else
//		{
//
//			sPtr.TransportSmarticle(bucket_ctr + ChVector<>(
//				sin(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(MyRand() - 1 / 2.0)),
//				cos(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(MyRand() - 1 / 2.0)),
//				bucket_interior_halfDim.z*1.75
//				));
//			//sPtr->TransportSmarticle(ChVector<>
//			//	(ChVector<>(sPtr->GetArm(1)->GetPos().x,
//			//	sPtr->GetArm(1)->GetPos().y,
//			//	bucket_interior_halfDim.z*1.75)));
//			//sPtr->SetSpeed(sPtr->GetArm(1)->GetPos_dt()/2);
//		}
//
//
//		recycledSmarticles++;
//		inc = (inc + 1) % 5;
//	}
//	printFlowRate(mphysicalSystem.GetChTime(), recycledSmarticles);
//}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//creates an approximate cylinder from a n-sided regular polygon
//num_boxes = number of boxes to use
//bucket_rad = radius of cylinder, center point to midpoint of side a side
//ChSharedPtr<ChBody> create_cylinder_from_blocks(int num_boxes, int id, bool overlap, CH_SYSTEM* mphysicalSystem, ChSharedPtr<ChMaterialSurfaceBase> wallMat)
//{
//	ChSharedPtr<ChBody> cyl_container;
//	if (USE_PARALLEL) {
//		cyl_container = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
//	}
//	else {
//		cyl_container = ChSharedPtr<ChBody>(new ChBody);
//	}
//	cyl_container->SetIdentifier(id);
//	//cyl_container->SetMass(mass);
//	cyl_container->SetPos(bucket_ctr);
//	cyl_container->SetRot(QUNIT);
//	cyl_container->SetBodyFixed(false);
//	cyl_container->SetCollide(true);
//	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
//	double wallt = t / 5; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
//	double half_height = bucket_interior_halfDim.z;
//	double box_side = bucket_rad * 2.0 * tan(CH_C_PI / num_boxes);//side length of cyl
//	double o_lap = 0;
//	if (overlap){ o_lap = t * 2; }
//	double ang = 2.0 * CH_C_PI / num_boxes;
//	ChVector<> box_size = (0, 0, 0); //size of plates
//	ChVector<> pPos = (0, 0, 0);  //position of each plate
//	ChQuaternion<> quat = QUNIT; //rotation of each plate
//	cyl_container->GetCollisionModel()->ClearModel();
//	cyl_container->SetMaterialSurface(wallMat);
//	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_pinkwhite.png"));
//	for (int i = 0; i < num_boxes; i++)
//	{
//
//		box_size = ChVector<>((box_side + wallt) / 2.0,
//			wallt,
//			half_height + o_lap);
//
//		pPos = bucket_ctr + ChVector<>(sin(ang * i) * (wallt + bucket_rad),
//			cos(ang*i)*(wallt + bucket_rad),
//			half_height);
//
//		quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));
//
//		//this is here to make half the cylinder invisible.
//		bool m_visualization = false;
//		if (ang*i < 3 * CH_C_PI / 4 || ang*i > 5 * CH_C_PI / 4)
//		{
//			m_visualization = true;
//			cyl_container->AddAsset(bucketTexture);
//		}
//		cyl_container->GetCollisionModel()->SetEnvelope(collisionEnvelope);
//		utils::AddBoxGeometry(cyl_container.get_ptr(), box_size, pPos, quat, m_visualization);
//
//	}
//	//Add ground piece
//	//
//	//utils::AddBoxGeometry(cyl_container.get_ptr(), Vector(bucket_rad, bucket_rad + t, t), Vector(0, 0, -t), QUNIT, true);
//
//	//checks top,bottom, and middle location
//	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos() + Vector(0,0,2 * bucket_interior_halfDim.z + 2.0 * bucket_half_thick), Q_from_AngAxis(CH_C_PI/2.0, VECT_Y),true);
//	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos(), Q_from_AngAxis(CH_C_PI / 2.0, VECT_Y));
//	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos() + Vector(0, 0, bucket_interior_halfDim.z), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
//
//	//ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, bucket_interior_halfDim.z);
//
//
//	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad + 2 * t, t, ChVector<>(0, 0, -t), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
//	//add up volume of bucket and multiply by rho to get mass;
//	double cyl_volume = CH_C_PI*(2 * box_size.z - 2 * t)*(2 * box_size.z - 2 * t)*((2 * bucket_rad + 2 * t)*(2 * bucket_rad + 2 * t) - bucket_rad*bucket_rad) + (CH_C_PI)*(bucket_rad + 2 * t)*(bucket_rad + 2 * t) * 2 * t;
//	cyl_container->SetMass(rho_cylinder*cyl_volume);
//
//	//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
//	cyl_container->GetCollisionModel()->BuildModel();
//
//	mphysicalSystem->AddBody(cyl_container);
//	return cyl_container;
//}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//ChSharedPtr<ChBody> Create_hopper(CH_SYSTEM* mphysicalSystem, ChSharedPtr<ChMaterialSurfaceBase> wallMat, double w1, double w2, double w3, double h1, double h2, bool overlap)
//{
//	ChSharedPtr<ChBody> hopper;
//	if (USE_PARALLEL) {
//		hopper = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
//	}
//	else {
//		hopper = ChSharedPtr<ChBody>(new ChBody);
//	}
//
//
//	double hw1 = w1;
//	double hw2 = w2;
//	double hw3 = w3;
//	double hh1 = h1*.5;
//	double hh2 = h2*.5;
//	double ht = bucket_half_thick;
//
//	hopper->SetPos(bucket_ctr);
//	hopper->SetRot(QUNIT);
//	hopper->SetBodyFixed(true);
//	hopper->SetCollide(true);
//
//
//	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
//	double o_lap = 0;
//	if (overlap){ o_lap = 2 * t; }
//
//	hopper->GetCollisionModel()->ClearModel();
//	hopper->SetMaterialSurface(wallMat);
//	double mtheta = atan((hw1 - hw3) / h1);
//
//
//	bucketTexture->SetTextureFilename(GetChronoDataFile("greenwhite.png"));
//	hopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);
//	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(hw1 + ht, 0, h1 + hh2), QUNIT, true); // upper part, max_x plate
//
//	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(-hw1 - ht, 0, h1 + hh2), QUNIT, true); // upper part, min_x plate
//	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, hw2 + ht, h1 + hh2), QUNIT, true); // upper part, min_x plate
//	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, -hw2 - ht, h1 + hh2), QUNIT, false); // upper part, min_x plate
//
//	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, -hw2 - ht, hh1), QUNIT, false); // upper part, min_x plate
//	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, hw2 + ht, hh1), QUNIT, true); // upper part, min_x plate
//
//	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(hw3 + hh1 * tan(mtheta) + ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(mtheta, VECT_Y), true); // upper part, min_x plate
//	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(-hw3 - hh1 * tan(mtheta) - ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(-mtheta, VECT_Y), true); // upper part, min_x plate
//	hopper->AddAsset(bucketTexture);
//
//	double estimated_volume = 8 * (w1 * t * h1); // Arman : fix this
//	hopper->SetMass(rho_cylinder*estimated_volume);
//	hopper->GetCollisionModel()->BuildModel();
//	mphysicalSystem->AddBody(hopper);
//	return hopper;
//}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//void InitializeMbdPhysicalSystem_Parallel(ChSystemParallelDVI& mphysicalSystem, int argc, char* argv[]) {
//	// initializd random seeder
//	MySeed(964);
//	// Desired number of OpenMP threads (will be clamped to maximum available)
//	int threads = 1;
//	// Perform dynamic tuning of number of threads?
//	bool thread_tuning = true;
//
//	//	uint max_iteration = 20;//10000;
//	int max_iteration_normal = 50;
//	int max_iteration_sliding = 50;
//	int max_iteration_spinning = 0;
//	int max_iteration_bilateral = 50;
//
//	// ----------------------
//	// Set params from input
//	// ----------------------
//
//	SetArgumentsForMbdFromInput(argc, argv, threads, max_iteration_sliding, max_iteration_bilateral, dT, numLayers, armAngle, read_from_file, pctActive, angle1, angle2);
//
//	// ----------------------
//	// Set number of threads.
//	// ----------------------
//
//	//
//	int max_threads = omp_get_num_procs();
//	if (threads > max_threads)
//		threads = max_threads;
//	mphysicalSystem.SetParallelThreadNumber(threads);
//	omp_set_num_threads(threads);
//
//	mphysicalSystem.GetSettings()->perform_thread_tuning = thread_tuning;
//	mphysicalSystem.GetSettings()->min_threads = std::max(1, threads / 2);
//	mphysicalSystem.GetSettings()->max_threads = std::min(max_threads, int(3.0 * threads / 2));
//	const std::string simulationParams = out_dir + "/simulation_specific_parameters.txt";
//	simParams.open(simulationParams.c_str(), std::ios::app);
//	// ---------------------
//	// Print the rest of parameters
//	// ---------------------
//	simParams << std::endl <<
//		" number of threads: " << threads << std::endl <<
//		" max_iteration_normal: " << max_iteration_normal << std::endl <<
//		" max_iteration_sliding: " << max_iteration_sliding << std::endl <<
//		" max_iteration_spinning: " << max_iteration_spinning << std::endl <<
//		" max_iteration_bilateral: " << max_iteration_bilateral << std::endl <<
//		" l_smarticle: " << l_smarticle << std::endl <<
//		" l_smarticle mult for w (w = mult x l): " << l_smarticle / w_smarticle << std::endl <<
//		" dT: " << dT << std::endl <<
//		" tFinal: " << tFinal << std::endl <<
//		" vibrate start: " << vibrateStart << std::endl <<
//		" read from file: " << read_from_file << std::endl <<
//		" arm angle: " << angle1 << " " << angle2 << std::endl << std::endl;
//
//
//	// ---------------------
//	// Edit mphysicalSystem settings.
//	// ---------------------
//
//	double tolerance = 0.001;  // 1e-3;  // Arman, move it to paramsH
//	mphysicalSystem.Set_G_acc(ChVector<>(0, 0, gravity));
//
//	mphysicalSystem.GetSettings()->solver.solver_mode = SLIDING;                              // NORMAL, SPINNING
//	mphysicalSystem.GetSettings()->solver.max_iteration_normal = max_iteration_normal;        // max_iteration / 3
//	mphysicalSystem.GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;      // max_iteration / 3
//	mphysicalSystem.GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;    // 0
//	mphysicalSystem.GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;  // max_iteration / 3
//	mphysicalSystem.GetSettings()->solver.tolerance = tolerance;
//	mphysicalSystem.GetSettings()->solver.alpha = 0;  // Arman, find out what is this
//	mphysicalSystem.GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
//	mphysicalSystem.ChangeSolverType(APGD);  // Arman check this APGD APGDBLAZE
//	//  mphysicalSystem.GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;
//
//	mphysicalSystem.GetSettings()->collision.collision_envelope = collisionEnvelope;
//	mphysicalSystem.GetSettings()->collision.bins_per_axis = _make_int3(40, 40, 40);  // Arman check
//}

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//void UpdateSmarticles(
//	CH_SYSTEM& mphysicalSystem,
//	Smarticle& sPtr) {
//
//	double current_time = mphysicalSystem.GetChTime();
//	sPtr.MoveLoop2(Smarticle::global_GUI_value);
//}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//void drawOTArms(Smarticle &sPtr)
//{
//	if (sPtr.GetArm0OT())
//	{
//		sPtr.GetArm(0)->AddAsset(Smarticle::mtextureOT);
//	}
//	else
//	{
//		sPtr.GetArm(0)->AddAsset(Smarticle::mtextureArm);
//	}
//	if (sPtr.GetArm2OT())
//	{
//		sPtr.GetArm(2)->AddAsset(Smarticle::mtextureOT);
//	}
//	else
//	{
//		sPtr.GetArm(2)->AddAsset(Smarticle::mtextureArm);
//	}

//	app->AssetBind(sPtr.GetArm(0));
//	app->AssetBind(sPtr.GetArm(2));
//	app->AssetUpdate(sPtr.GetArm(0));
//	app->AssetUpdate(sPtr.GetArm(2));
//}
//void drawOTArms()
//{
//	for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
//	{
//		Smarticle* sPtr = sv->at(i);

//		if (sPtr->GetArm0OT())
//		{
//			sPtr->GetArm(0)->AddAsset(Smarticle::mtextureOT);
//		}
//		else
//		{
//			sPtr->GetArm(0)->AddAsset(Smarticle::mtextureArm);
//		}
//		if (sPtr->GetArm2OT())
//		{
//			sPtr->GetArm(2)->AddAsset(Smarticle::mtextureOT);
//		}
//		else
//		{
//			sPtr->GetArm(2)->AddAsset(Smarticle::mtextureArm);
//		}

//		app->AssetBind(sPtr->GetArm(0));
//		app->AssetBind(sPtr->GetArm(2));
//		app->AssetUpdate(sPtr->GetArm(0));
//		app->AssetUpdate(sPtr->GetArm(2));
//	}
//}


// double Find_Z_Region_Heights(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mSmartVec)
// {
// 	double sqSize = w_smarticle/2; // size of squares in grid
// 	int rowSize = ceil(bucket_rad*2/sqSize);
// 	static std::vector<double> zHeights(rowSize*rowSize);
// 	double zmax = 0;
// 	for (size_t i = 0; i < mySmarticlesVec.size(); i++)
// 	{
// 		Smarticle* sPtr = mySmarticlesVec[i];
// 		zCom += sPtr->Get_cm().z-bucketMin.z;
//
// 	//isinradial rad parameter is Vector(bucketrad,zmin,zmax)
// 		ChVector<> pos = sPtr->Get_cm() - ChVector<>(0,0,bucket->GetPos());
// 		int xpos = int(pos.x/rowSize);
// 		int ypos = int(pos.y/rowSize);
// 		int vecPos = rowSize*xpos+y;
// 		zHeights[vecPos]=pos.z;
// 	}
// 	for (size_t i = 0; i < zHeights.size(); i++)
// 	{
// 		zmax= zmax+zHeights.at(i);
// 	}
// }
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//void SetEnvelopeForSystemObjects(ChSystem& mphysicalSystem) {
//	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
//	for (int i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
//		(*myIter)->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
//		myIter++;
//	}
//
//}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&