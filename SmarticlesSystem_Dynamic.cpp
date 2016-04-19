//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2010-2012 Alessandro Tasora
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

///////////////////////////////////////////////////
//
//   Demo code about
//
//     - Smarticle Simulation
//
//       (This is just a possible method of integration
//       of Chrono::Engine + Irrlicht: many others
//       are possible.)
//
//	 CHRONO
//   ------
//   Multibody dinamics engine
//
// ------------------------------------------------
//             www.deltaknowledge.com
// ------------------------------------------------
////*************** chrono parallel
#include <omp.h>
//#include "chrono_parallel/physics/ChSystemParallel.h"
//#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"
#include <iostream>
//#include <IStream>

#include "utils/ChUtilsCreators.h"     //Arman: why is this
#include "utils/ChUtilsInputOutput.h"  //Arman: Why is this
#include "utils/ChUtilsGenerators.h"
#include "common.h"
#include <ctime>
#include "core/ChFileutils.h" // for MakeDirectory
#include "IrrGui.h"
#include "Smarticle.h"
#include "SmarticleU.h"
#include "CheckPointSmarticles.h"
//#include <vld.h> //TODO used to find memory leaks
#include <memory>


#if irrlichtVisualization

#ifdef CHRONO_OPENGL
#undef CHRONO_OPENGL
#endif

//#include "unit_IRRLICHT/ChIrrApp.h"
#include "chrono_irrlicht/ChBodySceneNode.h"  //changed path from unit to chrono to reflect changes in updated chrono
#include "chrono_irrlicht/ChBodySceneNodeTools.h"
//#include "unit_IRRLICHT/ChIrrTools.h"
#include "chrono_irrlicht/ChIrrWizard.h"
#include "core/ChRealtimeStep.h"
//#include <irrlicht.h>
#include "assets/ChTexture.h"

using namespace chrono;
using namespace chrono::irrlicht;
using namespace irr;

using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;
using namespace irr::io;
using namespace irr::gui;
#endif


#ifdef CHRONO_OPENGL

#include "chrono_opengl/ChOpenGLWindow.h"
#endif



//#define GLEW_STATIC
//#include <GL/glew.h>
//#pragma comment(lib, "opengl32.lib")
//#pragma comment(lib, "glu32.lib")
#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#if USE_PARALLEL
#define CH_SYSTEM ChSystemParallelDVI
#else
#define CH_SYSTEM ChSystem
#endif

//***********************************
// Use the namespace of Chrono
//enum SmarticleType { SMART_ARMS, SMART_U };
//enum BucketType { KNOBCYLINDER, HOOKRAISE, STRESSSTICK, CYLINDER, BOX, HULL, RAMP, HOPPER, DRUM };
SmarticleType smarticleType = SMART_ARMS;//SMART_U;
BucketType bucketType = BOX;
std::vector<std::shared_ptr<ChBody>> sphereStick;
std::shared_ptr<ChBody> bucket;
std::shared_ptr<ChBody> bucket_bott;

double Find_Max_Z(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mSmartVec);
//double Find_Max_Z(CH_SYSTEM& mphysicalSystem);
std::ofstream simParams;
double sizeScale = 1;
int appWidth = 1280;
int appHeight = 720;
//double gravity = -9.81 * sizeScale;
double gravity = -9.81;
double vibration_freq = 30;
double fric =.3814; //keyboard box friction = .3814
double omega_bucket = 2 * PI * vibration_freq;  // 30 Hz vibration similar to Gravish 2012, PRL
double mGamma = 2.0 * gravity;
double vibration_amp = mGamma / (omega_bucket*omega_bucket);
unsigned int largeID = 10000000;


//double dT = std::min(0.001, 1.0 / vibration_freq / 200);;//std::min(0.0005, 1.0 / vibration_freq / 200);
double dT = 0.0005;//std::min(0.0005, 1.0 / vibration_freq / 200);
double contact_recovery_speed = .5* sizeScale;
double tFinal = 100;
double vibrateStart= 100;


double rho_cylinder = 1180.0;
//double rho_smarticle = 7850.0 / (sizeScale * sizeScale * sizeScale);
//double rho_cylinder = 1180.0 / (sizeScale * sizeScale * sizeScale);
std::shared_ptr<ChMaterialSurface> mat_g;
int numLayers = 100;
double armAngle = 90;
double sOmega = 5;  // smarticle omega

std::shared_ptr<ChLinkEngine> bucket_actuator;

double gaitChangeLengthTime = .5;


////////////////rescaled robot geometry (3.93) based on w_smarticle scaling
////////////////robot dim is l/w =1, w=.046 t=.031 t2=.021
#if stapleSize
	double bucket_rad = sizeScale*0.02;
	double w_smarticle = sizeScale * 0.0117; // sizeScale * 0.0117
	double l_smarticle = 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	double t_smarticle = sizeScale * .00127;
	double t2_smarticle = sizeScale * .0005;
	ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, 2*bucket_rad/sizeScale);
	double rho_smarticle = 7850.0;
#else
	
	double w_smarticle = sizeScale * 0.05316/ 1;
	double l_smarticle = 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	//real value
	double t_smarticle = sizeScale * .029982 / 1; //height of solar panels
	double t2_smarticle = sizeScale * .02122 / 1;
	double bucket_rad = sizeScale*w_smarticle*2; //3
	ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, 2 * bucket_rad / sizeScale);
	double rho_smarticle = 443.0;
	ChVector<>boxdim(.28/1.5, .55245, 2 * bucket_rad / 8);
#endif

	double p_gain = .2;   //.32           //.2
	double i_gain = .225;	 //.4 //.225       //.225
	double d_gain = 0.01; //.1  //.01      //.01

	// double t_smarticle 	= sizeScale * .00254;
	// double t2_smarticle	= sizeScale * .001;
double vol = (t2_smarticle) * (t_smarticle)* (w_smarticle + 2 * (l_smarticle));
////robot smarticle geometry
//double w_smarticle = 0.046; //4.6cm
//double l_smarticle = 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
//double t_smarticle = .03;
//double t2_smarticle = .021;

double collisionEnvelope = .1 * t2_smarticle;

bool bucket_exist = true;

int read_from_file = 0; //-1 write 0 do nothing 1 read
bool povray_output = false;
int out_fps = 30;
const std::string out_dir = "PostProcess";
const std::string pov_dir_mbd = out_dir + "/povFilesSmarticles";
int numPerLayer =5;
bool placeInMiddle = false;	/// if I want make a single smarticle on bottom surface
ChVector<> bucket_ctr = ChVector<>(0,0,0);
//ChVector<> Cbucket_interior_halfDim = sizeScale * ChVector<>(.05, .05, .025);
//double bucket_rad = sizeScale*0.034;
//double bucket_rad = sizeScale*0.02;
//double bucket_rad = sizeScale*0.022;
//	double bucket_rad = sizeScale*0.04;

std::vector<std::shared_ptr<ChBody>> bucket_bod_vec;

//ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(.1, .1, .05);
double bucket_half_thick = sizeScale * .005;
double percentToMoveToGlobal = .25;
double max_z = 0;
double rampInc = 1;
double drum_freq = 1;
double box_ang =-40 *D2R;
double drum_omega = drum_freq*2*PI;
double pctActive = 1.0;
double inc = 0.00001;
double angle1 = 90;
double angle2 = 90;
double vibAmp = 5 * D2R; //vibrate by some amount of degrees back and forth
double videoFrameInterval = 134;
auto bucketTexture = std::make_shared<ChTexture>();
auto sphereTexture = std::make_shared<ChTexture>();
auto groundTexture = std::make_shared<ChTexture>();
auto floorTexture = std::make_shared<ChTexture>();

// =====================================================================================================
class MyBroadPhaseCallback : public collision::ChBroadPhaseCallback {
public:
	/// Callback used to report 'near enough' pairs of models.
	/// This must be implemented by a child class of ChBroadPhaseCallback.
	/// Return false to skip narrow-phase contact generation for this pair of bodies.
	virtual bool BroadCallback(collision::ChCollisionModel* mmodelA,  ///< pass 1st model
		collision::ChCollisionModel* mmodelB)   ///< pass 2nd model
	{
		return (!(abs(mmodelA->GetPhysicsItem()->GetIdentifier() - mmodelB->GetPhysicsItem()->GetIdentifier()) < 2));
	}
};
class ChFunctionCustom : public ChFunction{
public:
	ChFunctionCustom(){ y = 0; y_dx = 0; y_dxdx = 0; }
	virtual ~ChFunctionCustom(){};
	void Copy(ChFunction* source) {
		Set_y(source->Get_y(0));
		Set_y_dx(source->Get_y_dx(0));
		Set_y_dxdx(source->Get_y_dxdx(0));
	}
	virtual ChFunction* new_Duplicate() {
		ChFunctionCustom* m_func;
		m_func = new ChFunctionCustom;
		m_func->Copy(this);
		return (m_func);
	}
	virtual int Get_Type() { return 1; }
	void Set_y(double x){ y = x; }
	void Set_y_dx(double x){ y_dx = x; }
	void Set_y_dxdx(double x){ y_dxdx = x; }
	virtual double Get_y(double x) {return y;}
	virtual double Get_y_dx(double x) {return y_dx;}
	virtual double Get_y_dxdx(double x) { return y_dxdx; }

private:
	double y;
	double y_dx;
	double y_dxdx;
};

class ext_force :public ChReportContactCallback{

public:
	double n_contact_force = 0;
	ChVector<> t_contact_force = (0, 0, 0);
	double m_contact_force = 0;
	double maxHeight = 0;
	virtual bool ReportContactCallback(
		const ChVector<>& pA,             ///< get contact pA
		const ChVector<>& pB,             ///< get contact pB
		const ChMatrix33<>& plane_coord,  ///< get contact plane coordsystem (A column 'X' is contact normal)
		const double& distance,           ///< get contact distance
		const ChVector<>& react_forces,   ///< get react.forces (if already computed). In coordsystem 'plane_coord'
		const ChVector<>& react_torques,  ///< get react.torques, if rolling friction (if already computed).
		ChContactable* contactobjA,  ///< get model A (note: some containers may not support it and could be zero!)
		ChContactable* contactobjB   ///< get model B (note: some containers may not support it and could be zero!)
		)
	{
		unsigned int ia = contactobjA->GetPhysicsItem()->GetIdentifier();// reports force BY ib ON ia.
		unsigned int ib = contactobjB->GetPhysicsItem()->GetIdentifier();
		//GetLog() << "Report Callback\n";

		//normal force
		if (ia >= largeID)
		{

			//GetLog() << "running method";
			//double a = react_forces.Length();
			this->m_contact_force += react_forces.Length();
			//TODO get normal forces only!

			//n_contact_force += react_forces.y;
			//GetLog() << "Normal Force: " << m_contact_force << "\n";
			//t_contact_force += Vector(react_forces.y, react_forces.x, react_forces.z); ///x(output)=y(system) y(output)=x(system)  z(output) = z(sys)
		}

		return true;
	}

};
// =============================================================================\

double showForce(CH_SYSTEM *msys)
{
		
		ext_force ef;
		msys->GetContactContainer()->ReportAllContacts(&ef);
		return ef.m_contact_force; //TODO return max height too
}
// =============================================================================
void MySeed(double s = time(NULL)) { srand(s); }
double MyRand() { return float(rand()) / RAND_MAX; }
// =============================================================================
void SetArgumentsForMbdFromInput(int argc, char* argv[], int& threads, int& max_iteration_sliding, int& max_iteration_bilateral, double& dt, int& num_layers, double& mangle,int& readFile,double& mpctActive, double& mangle1, double& mangle2) {
  if (argc > 1) {
	const char* text = argv[1];
	double mult_l = atof(text);
	l_smarticle = mult_l * w_smarticle;
  }
	if (argc > 2){
		const char* text = argv[2];
		dt = atof(text);
	}
	if (argc > 3){
		const char* text = argv[3];
		numLayers = atoi(text);
	}

	if (argc > 4){
		const char* text = argv[4];
		readFile = atoi(text);
	}
	if (argc > 5){
		const char* text = argv[5];
		mpctActive = atof(text);
	}
	if (argc > 6){
		const char* text = argv[6];
		angle1 = atof(text);
		//angle1 = angle1*D2R;
	}
	if (argc > 7){
		const char* text = argv[7];
		angle2 = atof(text);
		//angle2 = angle2*D2R;
	}
	if (argc > 8){
		const char* text = argv[8];
		box_ang = atof(text)*D2R;
	}
	if (argc > 9){
		const char* text = argv[9];
		numPerLayer = atoi(text);
	}
	/// if parallel, get solver setting
  //if (USE_PARALLEL) {
	 // if (argc > 8) {
		//const char* text = argv[8];
		//threads = atoi(text);
	 // }
	 // if (argc > 9) {
		//const char* text = argv[9];
		//max_iteration_sliding = atoi(text);
	 // }
	 // if (argc > 10) {
		//const char* text = argv[10];
		//max_iteration_bilateral = atoi(text);
	 // }
  //}
}
// =============================================================================
void InitializeMbdPhysicalSystem_NonParallel(ChSystem& mphysicalSystem, int argc, char* argv[]) {
	// initializd random seeder
	MySeed();
	ChSetRandomSeed(time(NULL));

  // ---------------------
  // Print the rest of parameters
  // ---------------------
	const std::string simulationParams = out_dir + "/simulation_specific_parameters.txt";
	simParams.open(simulationParams.c_str(), std::ios::app);
	int dummyNumber0;
	int dummyNumber1;
	int dummyNumber2;
	int max_threads = omp_get_num_procs();
	int threads = 1;
	if (threads > max_threads)
		threads = max_threads;
	mphysicalSystem.SetParallelThreadNumber(threads);
	omp_set_num_threads(threads);


  SetArgumentsForMbdFromInput(argc, argv, dummyNumber0, dummyNumber1, dummyNumber2, dT,numLayers, armAngle, read_from_file,pctActive,angle1,angle2);
	vol = (t2_smarticle) * (t_smarticle)* (w_smarticle + 2 * (l_smarticle));
	simParams << std::endl <<
		"l_smarticle: " << l_smarticle << std::endl <<
		"l_smarticle mult for w (w = mult x l): " << l_smarticle / w_smarticle << std::endl <<
		"read from file: " << read_from_file << std::endl <<
		"dT: " << dT << std::endl << std::endl <<
		"tFinal: " << tFinal << std::endl <<
		"vibrate start: " << vibrateStart << std::endl <<
		"Active Percent: " << pctActive << std::endl <<
		"Start Angles: " << angle1 << " " << angle2<< std::endl;

	simParams << "Smarticle volume: " << vol << std::endl;
	simParams << "Smarticle mass: " << vol*rho_smarticle << std::endl;

	//copy smarticle checkpoint if used to PostProcess folder
	if (read_from_file==1)
	{
		const std::string copyCheckpoint = std::string("cp smarticles.csv " + out_dir);
		std::system(copyCheckpoint.c_str());
	}
	//copy smarticleMoves to PostProcess Folder
	const std::string copyMoves = std::string("cp smarticleMoves.csv " + out_dir);
	std::system(copyMoves.c_str());

  // ---------------------
  // Edit mphysicalSystem settings.
  // ---------------------

  // Modify some setting of the physical system for the simulation, if you want
	mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR);
	//mphysicalSystem.SetIntegrationType(ChSystem::INT_EULER_IMPLICIT_PROJECTED);
	mphysicalSystem.SetIterLCPmaxItersSpeed(80+1.1*numLayers);
  mphysicalSystem.SetIterLCPmaxItersStab(0);   // unuseful for Anitescu, only Tasora uses this
  mphysicalSystem.SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  mphysicalSystem.SetIterLCPwarmStarting(true);
  mphysicalSystem.SetUseSleeping(false);
  mphysicalSystem.Set_G_acc(ChVector<>(0, 0, gravity));
	vol = (t2_smarticle)* (t_smarticle)* (w_smarticle + 2 * (l_smarticle));
	//mphysicalSystem.SetTolForce(.0005);
	//mphysicalSystem.SetTol(.0001);
	//mphysicalSystem.SetMinBounceSpeed(.3);
	simParams.close();
}
// =============================================================================

#if irrlichtVisualization
void AddParticlesLayer1(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec, ChIrrApp& application,double timeForDisp) {
#else
void AddParticlesLayer1(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec,double timeForDisp) {
#endif
	
	ChVector<> dropSpeed = VNULL;
	ChQuaternion<> myRot = QUNIT;
	double z;
	double zpos;
	size_t smarticleCount = mySmarticlesVec.size();
	double ang = 2*PI / numPerLayer;
	double w = w_smarticle;
	if (smarticleCount < numPerLayer){ z = w_smarticle / 1; }
	//else{ z = max_z; }
	else{ z = Find_Max_Z(mphysicalSystem,mySmarticlesVec); }
	double phase = genRand(PI_2);
	ChVector<> myPos;
	for (int i = 0; i < numPerLayer; i++)
	{
		phase = genRand(PI_2);
		zpos = std::min(3 * bucket_interior_halfDim.z, z) + w_smarticle;

		switch (bucketType){
		case DRUM:
			zpos = SaturateValue(z, bucket_rad);
			myPos = bucket_ctr + ChVector<>(genRand(bucket_interior_halfDim.z / 2.5),
				genRand(bucket_interior_halfDim.z / 2.5),
				zpos);
			break;
		case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
			if (!placeInMiddle)
			{
				myPos = bucket_ctr + ChVector<>(sin(ang * i + phase) *(bucket_rad / 2.2),
				cos(ang*i + phase)*(bucket_rad / 2.2),
				std::max(bucket_interior_halfDim.z*2,zpos));
				dropSpeed = ChVector<>(0, 0, gravity*timeForDisp / 2.0 - 2 * w_smarticle / timeForDisp);
				myRot = ChQuaternion<>(genRand(-1, 1), genRand(-1, 1), genRand(-1, 1), genRand(-1, 1));
			}
			else////////////place in center of bucket on bucket bottom
			{
				myPos = bucket_ctr + ChVector<>(0,-t_smarticle*1.45,bucket_bott->GetPos().z + t_smarticle );
				dropSpeed = VNULL;
				myRot = Q_from_AngAxis(-PI_2, VECT_X);
			}
			////////////////////////////////////



			break;
		case HOPPER:
			myPos = bucket_ctr + ChVector<>(sin(ang * i + phase) *(bucket_rad / 2 + genRand(w)),
				cos(ang*i + phase)*(bucket_rad / 2 + genRand(w)),
				zpos);
			break;
		case BOX:
		{
			//myPos = bucket_ctr + ChVector<>((2*MyRand()-1)*.9*bucket_interior_halfDim.x,
			//	(2*MyRand()-1)*.9*bucket_interior_halfDim.y ,
			//	bucket_interior_halfDim.z+w_smarticle);
			//myRot = ChQuaternion<>(2 * MyRand() - 1, 2 * MyRand() - 1, 2 * MyRand() - 1, 2 * MyRand() - 1);
			////dropSpeed = ChVector<>(0, 0, gravity*timeForDisp / 2.0 - 2 * w_smarticle / timeForDisp);
			//dropSpeed = VNULL;
			///////////////////////////place in specific location///////////////////////////
			dropSpeed = VNULL;
			myRot = Q_from_AngAxis(PI_2, VECT_X);
			ChQuaternion<> buckRot = bucket->GetRot();
			double buckRotAngx = Quat_to_Angle(ANGLESET_RXYZ,buckRot).x;
			//myRot = buckRot*Angle_to_Quat(ANGLESET_RXYZ,ChVector<>(PI/2,PI,0));
			//myRot = buckRot*Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(genRand(0, 2*PI), genRand(0, 2*PI), genRand(0, 2*PI)));
			
			myRot = buckRot*Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(PI/2,0, 0));
			double xPos = genRand(-3,3)*t2_smarticle/1.25;
			double yPos = (i-4.2) * 2 * t2_smarticle;
			myPos = bucket_ctr + ChVector<>(xPos, yPos, ( - yPos - 2*bucket_half_thick)*tan(buckRotAngx)+t_smarticle/1.99);
			//myPos = bucket_ctr + ChVector<>(xPos, yPos, (-yPos - 2 * bucket_half_thick)*tan(buckRotAngx) + 3*t_smarticle) ;
			///////////////////////////place in specific location///////////////////////////
			break;
		}

		default:
			myPos = bucket_ctr + ChVector<>(sin(ang * i + phase) *(bucket_rad / 2 + genRand(-w/2,w/2)),
				cos(ang*i + phase)*(bucket_rad / 2 + genRand(-w / 2, w / 2)),
				zpos);
			myRot = ChQuaternion<>(genRand(-1, 1), genRand(-1, 1), genRand(-1, 1), genRand(-1, 1));
			break;
		}
			
			myRot.Normalize();
			/////////////////flat rot/////////////////
			Smarticle * smarticle0 = new Smarticle(&mphysicalSystem);
			smarticle0->Properties(mySmarticlesVec.size(), mySmarticlesVec.size() * 4,
				rho_smarticle, mat_g,
				collisionEnvelope,
				//l_smarticle+t2_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
				l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
				sOmega,
				true,
				myPos,
				myRot,
				angle1*D2R,angle2*D2R);

			if (genRand()<1)//to reduce amount visualized amount
				smarticle0->visualize = true;
			smarticle0->populateMoveVector();
			smarticle0->SetAngles(angle1, angle2,true);
			smarticle0->SetInitialAngles();
			smarticle0->Create();
			smarticle0->setCurrentMoveType((MoveType) Smarticle::global_GUI_value);
			smarticle0->vib.emplace_back(angle1*D2R, angle2*D2R);

			//must put this line because of how linspace is written;
			smarticle0->AssignState(VIB);
			smarticle0->GenerateVib(angle1*D2R, angle2*D2R);
			smarticle0->AssignState(Smarticle::global_GUI_value);
			//smarticle0->ss.emplace_back(angle1, angle2);
			//smarticle0->midTorque.emplace_back(angle1*D2R + vibAmp, angle2*D2R + vibAmp);
			//smarticle0->midTorque.emplace_back(angle1*D2R + vibAmp, angle2*D2R + vibAmp);

			mySmarticlesVec.emplace_back((Smarticle*)smarticle0);
			smarticle0->SetSpeed(dropSpeed);
			
#if irrlichtVisualization
			application.AssetBindAll();
			application.AssetUpdateAll();
			//application.AssetUpdate(smarticle0->GetArm(0));
			//application.AssetUpdate(smarticle0->GetArm(1));
			//application.AssetUpdate(smarticle0->GetArm(2));
			
#endif
	}


}

std::shared_ptr<ChBody> create_drum(int num_boxes, int id, bool overlap, CH_SYSTEM* mphysicalSystem, std::shared_ptr<ChMaterialSurfaceBase> wallMat,int ridges = 5)
{//essentially the same as create cyl container except made it bigger and added ridges

	auto drum = std::make_shared<ChBody>();

	double radMult =1.5;
	drum->SetIdentifier(id);
	drum->SetPos(bucket_ctr);
	drum->SetRot(QUNIT);
	drum->SetBodyFixed(false);
	drum->SetCollide(true);
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double wallt = t / 5.0; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
	double half_height = bucket_interior_halfDim.z/(radMult*2);
	double box_side = bucket_rad * radMult*2 * tan(PI / num_boxes);//side length of cyl
	double o_lap = 0;
	if (overlap){ o_lap = t * 2; }
	double ang = 2.0 * PI / num_boxes;
	ChVector<> box_size = (0, 0, 0); //size of plates
	ChVector<> ridge_size = (0, 0, 0); //size of plates
	ChVector<> pPos = (0, 0, 0);  //position of each plate
	ChQuaternion<> quat = QUNIT; //rotation of each plate
	drum->GetCollisionModel()->ClearModel();
	drum->SetMaterialSurface(wallMat);
	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
	int ridgeNum = num_boxes/ridges;
	for (int i = 0; i < num_boxes; i++)
	{

		box_size = ChVector<>((box_side + wallt) / 2.0,
			wallt,
			half_height + o_lap);
		pPos = bucket_ctr + ChVector<>(sin(ang * i) * (wallt + bucket_rad*radMult),
			cos(ang*i)*(wallt + bucket_rad*radMult),
			0);

		quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));

		drum->AddAsset(bucketTexture);

		drum->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		utils::AddBoxGeometry(drum.get(), box_size, pPos, quat);



		if (i%ridgeNum == 0)
		{
			ridge_size = ChVector<>((box_side) / 4.0,
				w_smarticle/8.0,
				half_height + o_lap);
			pPos = bucket_ctr + ChVector<>(sin(ang * i) * (-wallt + bucket_rad*radMult),
				cos(ang*i)*(-wallt + bucket_rad*radMult),
				0);

			drum->GetCollisionModel()->SetEnvelope(collisionEnvelope);
			//utils::AddBoxGeometry(drum.get(), ridge_size, pPos, quat);
		}
		drum->SetRot(Q_from_AngAxis(PI_2, VECT_X));


	}
//TODO add bucketVolume as global variable and set it in each function to calculate for each shape volumefraction seamlessly
	//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);

	drum->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	//front wall made invisible so we can see inside
	utils::AddBoxGeometry(drum.get(), ChVector<>(wallt + bucket_rad*radMult, wallt + bucket_rad*radMult, wallt), bucket_ctr + VECT_Z*(half_height + 2 * o_lap - 2 * t), QUNIT,false);
	drum->GetCollisionModel()->BuildModel();
	drum->GetCollisionModel()->SetFamily(1);
	drum->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	mphysicalSystem->AddBody(drum);


	bucket_bott->SetBodyFixed(true);
	bucket_bott->SetCollide(true);
	bucket_bott->GetCollisionModel()->ClearModel();
	bucket_bott->SetPos(bucket_ctr);
	bucket_bott->SetMaterialSurface(mat_g);
	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));//custom file
	bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	bucket_bott->AddAsset(bucketTexture);
	utils::AddBoxGeometry(bucket_bott.get(), ChVector<>(wallt + bucket_rad*radMult, wallt + bucket_rad*radMult, wallt), bucket_ctr - VECT_Z*(half_height + 2 * o_lap - 2 * t), QUNIT);


	bucket_bott->GetCollisionModel()->BuildModel();
	mphysicalSystem->AddBody(bucket_bott);
	bucket_bott->SetRot(Q_from_AngAxis(PI_2, VECT_X));
	bucket_bott->GetCollisionModel()->SetFamily(1);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	return drum;
}
std::shared_ptr<ChBody> create_complex_convex_hull(CH_SYSTEM* mphysicalSystem, std::shared_ptr<ChMaterialSurfaceBase> wallMat, double numBoxes)
{

	auto convexShape = std::make_shared<ChBody>();
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code

	//cyl_container->SetMass(mass);
	convexShape->SetPos(bucket_ctr);
	convexShape->SetRot(QUNIT);
	convexShape->SetBodyFixed(false);
	convexShape->SetCollide(true);

	std::vector<ChVector<>> points;

	convexShape->GetCollisionModel()->ClearModel();

	double ang = 2 * PI / numBoxes;
	for (size_t i = 0; i < numBoxes; i++)
	{
		convexShape->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	}

	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));

	ChVector<> rampSize(w_smarticle * 5, w_smarticle * 5, t);


	//utils::AddBoxGeometry(convexShape.get(), rampSize, rampPos, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(box_ang, 0, 0)), true);



	convexShape->GetCollisionModel()->BuildModel();
	convexShape->AddAsset(floorTexture);
	mphysicalSystem->AddBody(convexShape);
	return convexShape;
}
std::shared_ptr<ChBody> create_ramp(int id, CH_SYSTEM* mphysicalSystem, std::shared_ptr<ChMaterialSurfaceBase> wallMat)
{
	std::shared_ptr<ChBody> ramp;
	ramp = std::make_shared<ChBody>();
	double t = bucket_half_thick/2; //bucket thickness redefined here for easier to read code
	double w = w_smarticle *3;
	double h = w/2;

	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
	//cyl_container->SetMass(mass);
	ramp->SetPos(bucket_ctr);
	ramp->GetCollisionModel()->SetFamily(1);
	ramp->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	ramp->SetBodyFixed(false);
	ramp->SetCollide(true);

	ChVector<> rampPos(0, 0, -sin(box_ang)*w - t);

	ramp->GetCollisionModel()->ClearModel();
	ramp->SetMaterialSurface(wallMat);



	ramp->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	ramp->GetCollisionModel()->SetEnvelope(collisionEnvelope);

	utils::AddBoxGeometry(ramp.get(), ChVector<>(w, w, t), rampPos, QUNIT, true);//bottom
	utils::AddBoxGeometry(ramp.get(), ChVector<>(t, w, h), rampPos - VECT_X*(w - t) + VECT_Z*(h - t), QUNIT, true);// -x
	utils::AddBoxGeometry(ramp.get(), ChVector<>(t, w, h), rampPos + VECT_X*(w - t) + VECT_Z*(h - t), QUNIT, true);// +x
	utils::AddBoxGeometry(ramp.get(), ChVector<>(w, t, h), rampPos - VECT_Y*(w - t) + VECT_Z*(h - t), QUNIT, true);//high side -y
	//utils::AddBoxGeometry(ramp.get(), ChVector<>(w, t, h), rampPos + VECT_Y*(w - t) + VECT_Z*(h - t), QUNIT, true); //down side +y


	///bucket_bott in ramp is the +y side!
	bucket_bott->SetBodyFixed(true);
	bucket_bott->SetCollide(true);
	bucket_bott->GetCollisionModel()->ClearModel();
	bucket_bott->SetPos(bucket_ctr);
	bucket_bott->SetMaterialSurface(mat_g);
	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));//custom file
	bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(bucket_bott.get(), ChVector<>(w, t, h), rampPos + VECT_Y*(w - t) + VECT_Z*(h - t), QUNIT, true); //down side +y
	bucket_bott->AddAsset(floorTexture);
	bucket_bott->GetCollisionModel()->SetFamily(1);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	bucket_bott->GetCollisionModel()->BuildModel();
	bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(box_ang, 0, 0)));
	mphysicalSystem->AddBody(bucket_bott);

	ramp->GetCollisionModel()->BuildModel();
	ramp->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(box_ang, 0, 0)));
	ramp->AddAsset(floorTexture);
	mphysicalSystem->AddBody(ramp);
	return ramp;
}
// =============================================================================


std::shared_ptr<ChBody> create_cylinder_from_blocks2(int num_boxes, int id, bool overlap, CH_SYSTEM* mphysicalSystem, std::shared_ptr<ChMaterialSurfaceBase> wallMat)
{
	auto cyl_container = std::make_shared<ChBody>();
	cyl_container->SetIdentifier(id);
	//cyl_container->SetMass(mass);
	cyl_container->SetPos(bucket_ctr);
	cyl_container->SetRot(QUNIT);
	cyl_container->SetBodyFixed(false);
	cyl_container->SetCollide(true);
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double wallt = t / 5; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
	double half_height = bucket_interior_halfDim.z;
	double box_side = bucket_rad * 2.0 * tan(PI / num_boxes);//side length of cyl
	double o_lap = 0;
	if (overlap){ o_lap = t * 2; }
	double ang = 2.0 * PI / num_boxes;
	ChVector<> box_size = (0, 0, 0); //size of plates
	ChVector<> pPos = (0, 0, 0);  //position of each plate
	ChQuaternion<> quat = QUNIT; //rotation of each plate
	cyl_container->GetCollisionModel()->ClearModel();
	cyl_container->SetMaterialSurface(wallMat);
	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_pinkwhite.png"));
	for (int i = 0; i < num_boxes; i++)
	{

		box_size = ChVector<>((box_side + wallt) / 2.0,
			wallt,
			half_height + o_lap);

		pPos = bucket_ctr + ChVector<>(sin(ang * i) * (wallt + bucket_rad),
			cos(ang*i)*(wallt + bucket_rad),
			half_height);

		quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));

		//this is here to make half the cylinder invisible.
		bool m_visualization = false;
		if (ang*i < 3 * PI / 4 || ang*i > 5 * PI / 4)
		{
			m_visualization = true;
			cyl_container->AddAsset(bucketTexture);
		}
		cyl_container->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		utils::AddBoxGeometry(cyl_container.get(), box_size, pPos, quat, m_visualization);

	}

	double cyl_volume = PI*(2 * box_size.z - 2 * t)*(2 * box_size.z - 2 * t)*((2 * bucket_rad + 2 * t)*(2 * bucket_rad + 2 * t) - bucket_rad*bucket_rad) + (PI)*(bucket_rad + 2 * t)*(bucket_rad + 2 * t) * 2 * t;
	cyl_container->SetMass(rho_cylinder*cyl_volume);

	//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	cyl_container->GetCollisionModel()->BuildModel();


	for (int i = 0; i < num_boxes; i++)
	{
		auto wallPiece = std::make_shared<ChBody>();
		box_size = ChVector<>((box_side + wallt) / 2.0,
			wallt,
			half_height + o_lap);

		pPos = bucket_ctr + ChVector<>(sin(ang * i) * (wallt + bucket_rad),
			cos(ang*i)*(wallt + bucket_rad),
			half_height);

		quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));

		//this is here to make half the cylinder invisible.
		bool m_visualization = false;
		if (ang*i < 3 * PI / 4 || ang*i > 5 * PI / 4)
		{
			m_visualization = true;
			wallPiece->AddAsset(bucketTexture);
		}
		wallPiece->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		wallPiece->SetPos(pPos);
		wallPiece->GetCollisionModel()->ClearModel();
		utils::AddBoxGeometry(wallPiece.get(), box_size, ChVector<>(0,0,0), quat, m_visualization);
		wallPiece->SetMass(rho_cylinder*cyl_volume);
		//wallPiece->SetPos(ChVector<>(0,0,0));


		//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
		wallPiece->GetCollisionModel()->BuildModel();
		wallPiece->GetCollisionModel()->SetFamily(1);
		wallPiece->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
		mphysicalSystem->AddBody(wallPiece);

		wallPiece->SetRot(QUNIT);
		wallPiece->SetBodyFixed(true);
		wallPiece->SetCollide(true);
		bucket_bod_vec.emplace_back(wallPiece);

	}

	return cyl_container;
}

// =============================================================================
std::shared_ptr<ChBody> Create_hopper2(CH_SYSTEM* mphysicalSystem, std::shared_ptr<ChMaterialSurfaceBase> wallMat, double theta, double holeSize, bool overlap)
{
	auto hopper = std::make_shared<ChBody>();

	holeSize = holeSize / 2;
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double r = bucket_rad;
	double ang = theta*D2R;
	double h = bucket_interior_halfDim.z;
	double sH = (h-t) / sin(ang); //side height
	//double w = holeSize / 2 + cos(ang)*h;
	double w =cos(ang)*h-t/2;
	hopper->SetPos(bucket_ctr);
	hopper->SetRot(QUNIT);
	hopper->SetBodyFixed(true);
	hopper->SetCollide(true);



	double o_lap = 0;
	if (overlap){ o_lap = 2 * t; }

	hopper->GetCollisionModel()->ClearModel();
	hopper->SetMaterialSurface(wallMat);

	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_black_bordersBlack.png"));
	hopper->AddAsset(bucketTexture);
	hopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);


	ChVector<> hop_pos = bucket_ctr;
	hopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	//bucket_rad + 2 * bucket_half_thick, bucket_rad + 2 * bucket_half_thick, bucket_half_thick
	//utils::AddBoxGeometry(ramp.get(), ChVector<>(w, w, t), rampPos, QUNIT, true);//bucket_bottom


	utils::AddBoxGeometry(hopper.get(), ChVector<>(t, r + t, h), hop_pos - VECT_X*(sin(ang)*h + holeSize+cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, ang, 0)), true);// -x negAngle
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(t, r + t, h), hop_pos - VECT_X*(sin(ang)*(h)+holeSize+cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, ang, 0)), true);// -x negAngle +theta
	utils::AddBoxGeometry(hopper.get(), ChVector<>(t, r + t, h), hop_pos + VECT_X*(sin(ang)*h + holeSize+cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, -ang, 0)), true);// +x  -theta = /


	utils::AddBoxGeometry(hopper.get(), ChVector<>(2 * sin(ang)*h + holeSize+2*cos(ang)*t, t, h), hop_pos - VECT_Y*(r + o_lap) + VECT_Z*(h - o_lap), QUNIT, false);//camera pos -y
	utils::AddBoxGeometry(hopper.get(), ChVector<>(2 * sin(ang)*h + holeSize+2*cos(ang)*t, t, h), hop_pos + VECT_Y*(r + o_lap) + VECT_Z*(h - o_lap), QUNIT, true); //camera target +y

	hopper->GetCollisionModel()->BuildModel();
	mphysicalSystem->AddBody(hopper);
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(r, t, h + o_lap), ChVector<>(-r/2, 0, h+o_lap), QUNIT, true); // front plate, max_x plate

	//utils::AddBoxGeometry(hopper.get(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(-hw1 - ht, 0, h1 + hh2), QUNIT, true); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, hw2 + ht, h1 + hh2), QUNIT, true); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, -hw2 - ht, h1 + hh2), QUNIT, false); // upper part, min_x plate

	//utils::AddBoxGeometry(hopper.get(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, -hw2 - ht, hh1), QUNIT, false); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, hw2 + ht, hh1), QUNIT, true); // upper part, min_x plate

	//utils::AddBoxGeometry(hopper.get(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(hw3 + hh1 * tan(mtheta) + ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(mtheta, VECT_Y), true); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(-hw3 - hh1 * tan(mtheta) - ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(-mtheta, VECT_Y), true); // upper part, min_x plate
	//hopper->AddAsset(bucketTexture);
	return hopper;
}
void CreateBucket_bott(CH_SYSTEM& mphysicalSystem)
{
	bucket_bott->SetBodyFixed(true);
	bucket_bott->SetCollide(true);
	bucket_bott->GetCollisionModel()->ClearModel();
	bucket_bott->SetPos(bucket_ctr);
	bucket_bott->SetMaterialSurface(mat_g);
	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));//custom file
	bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(bucket_bott.get(), Vector(bucket_rad + 2 * bucket_half_thick, bucket_rad + 2 * bucket_half_thick, bucket_half_thick), Vector(0, 0, -bucket_half_thick), QUNIT, true);
	bucket_bott->AddAsset(floorTexture);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	bucket_bott->GetCollisionModel()->BuildModel();
	mphysicalSystem.AddBody(bucket_bott);
}

// =============================================================================
void CreateMbdPhysicalSystemObjects(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec) {
	/////////////////
	// Ground body
	////////////////
	mat_g->SetFriction(fric); //steel- plexiglass   (plexiglass was outer cylinder material)
	mat_g->SetKfriction(fric);
	// ground
	ChVector<> boxDim = sizeScale * ChVector<>(0.1, 0.1, .002);
	ChVector<> boxLoc = sizeScale * ChVector<>(0, 0, -5.0*bucket_interior_halfDim.z);
	auto ground = std::make_shared<ChBody>();
	bucket = std::make_shared<ChBody>();
	bucket_bott = std::make_shared<ChBody>();
	ground->SetMaterialSurface(mat_g);
	ground->SetPos(boxLoc);

	// ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(true);

	ground->GetCollisionModel()->ClearModel();
	ground->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddCylinderGeometry(ground.get(), boxDim.x, boxDim.z, ChVector<>(0, 0, 0), Q_from_AngAxis(PI_2, VECT_X));
	ground->GetCollisionModel()->SetFamily(1);
	ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	ground->GetCollisionModel()->BuildModel();
	mphysicalSystem.AddBody(ground);

	groundTexture->SetTextureFilename(GetChronoDataFile("greenwhite.png"));
	ground->AddAsset(groundTexture);

	// 1: create bucket
	ChVector<> dim = bucket_interior_halfDim;
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));
		switch (bucketType)		//http://www.engineeringtoolbox.com/friction-coefficients-d_778.html to get coefficients
		{
		case BOX:
		{
			//dim = (2 * dim.x, 2 * dim.y, dim.z / 8);
			//dim.x = dim.x;
			dim = (.28 / 2);//28 cm across
			dim.y = .55245;
			dim.z = dim.z / 8;
			bucket = utils::CreateBoxContainer(&mphysicalSystem, 1000000000, mat_g,
				dim, bucket_half_thick, bucket_ctr, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(box_ang, 0, 0)), true, false, true, false);
			//bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
			bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));
			bucket->AddAsset(bucketTexture);
			bucket->SetCollide(true);
			bucket->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
			bucket_bott->GetCollisionModel()->SetFamily(1);
			bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
		
			break;
		}
		case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
			CreateBucket_bott(mphysicalSystem);
			bucket = create_cylinder_from_blocks2(25, 1, true, &mphysicalSystem, mat_g);
			break;

		case RAMP:
			bucket = create_ramp(1,&mphysicalSystem,mat_g);
			break;

		case HOPPER:
			//bucket = Create_hopper(&mphysicalSystem, mat_g, bucket_interior_halfDim.x*2, bucket_interior_halfDim.y, 0.5 * bucket_interior_halfDim.x, bucket_interior_halfDim.z, 2 * bucket_interior_halfDim.z, true);
			bucket = Create_hopper2(&mphysicalSystem,mat_g,30,2*w_smarticle,true);
			CreateBucket_bott(mphysicalSystem);
			break;

		case HULL:
			break;

		case DRUM:
			bucket = create_drum(25, 1, true, &mphysicalSystem, mat_g);
		}
		//mat_g->SetFriction(fric); //steel - steel


	bucket->SetBodyFixed(true);
	bucket->SetCollide(true);
	bucket->GetCollisionModel()->SetFamily(1);
	bucket->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
}
// =============================================================================

void SavePovFilesMBD(CH_SYSTEM& mphysicalSystem,
                     int tStep) {
  int out_steps = std::ceil((1.0 / dT) / out_fps);
  //printf("tStep %d , outstep %d, num bodies %d chrono_time %f\n", tStep, out_steps, mphysicalSystem.Get_bodylist()->size(), mphysicalSystem.GetChTime());

  static int out_frame = 0;

  // If enabled, output data for PovRay postprocessing.
  if (povray_output && tStep % out_steps == 0) {
    char filename[100];
    sprintf(filename, "%s/data_%03d.dat", pov_dir_mbd.c_str(), out_frame + 1);
    utils::WriteShapesPovray(&mphysicalSystem, filename);

    out_frame++;
  }
}
// =============================================================================

double Find_Max_Z(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mySmarticlesVec) {
	std::string smarticleTypeName;
	if (smarticleType == SMART_ARMS) {
		smarticleTypeName = "smarticle_arm";
	}
	else if (smarticleType == SMART_U) {
		smarticleTypeName = "smarticle_u";
	}
	else {
		std::cout << "Error! Smarticle type is not set correctly" << std::endl;
	}
	double zMax = -999999999;
	//std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
	std::vector<std::shared_ptr<ChBody> >::iterator ibody = mphysicalSystem.Get_bodylist()->begin();

	for (size_t i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
		//ChBody* bodyPtr = *(myIter + i);
		std::shared_ptr<ChBody> bodyPtr = *(ibody + i);
		if (strcmp(bodyPtr->GetName(), smarticleTypeName.c_str()) == 0) {
			if (zMax < bodyPtr->GetPos().z) {
				//zMax = bodyPtr->GetPos().z;
				zMax = bodyPtr->GetPos().z - bucket_bott->GetPos().z;
			}
		}
	}
	return zMax;
}
// =============================================================================
bool IsIn(ChVector<> pt, ChVector<> min, ChVector<> max) {

	if ((pt < max) && (pt > min)) {
		return true;
	}

	return false;
}
// =============================================================================
//used for now cylinder
// isinradial rad parameter is Vector(bucketrad, zmin, zmax)
bool IsInRadial(ChVector<> pt, ChVector<> centralPt, ChVector<> rad)
{
	ChVector<> dist = pt - centralPt;
	double xydist = std::sqrt(dist.x * dist.x + dist.y * dist.y);
	//xydist = (std::sqrt((pt.x - centralPt.x)*(pt.x - centralPt.x) + (pt.y - centralPt.y)*(pt.y - centralPt.y)));
//	if ((xydist < rad.x) && (pt.z >= rad.y) && (pt.z < rad.z)) {
//		return true;
//	}
//	return false;
	if (xydist >= rad.x) { /*GetLog() << "\noutside radius\n";*/ return false; } // if outside radius
	if (pt.z < rad.y || pt.z >rad.z){ /*GetLog() <<  "outside z";*/ return false; }
	return true;
}
void printFlowRate(double time,int count) //SAVE smarticle gaitType out for reference!
{
	static bool started= false;
	const std::string flowRate = out_dir + "/flowrate.txt";

	std::ofstream flowRate_of;
	if (!started)
	{
		flowRate_of.open(flowRate.c_str());
		started = true;
	}
	else
	{
		flowRate_of.open(flowRate.c_str(), std::ios::app);
	}
	flowRate_of << time << ", " << count << ", " << Smarticle::global_GUI_value<<std::endl;

	flowRate_of.close();
}
// =============================================================================
void drawGlobalCoordinateFrame(CH_SYSTEM& mphysicalSystem)
{	
	double len = w_smarticle*2;
	double rad = t_smarticle/2;
	ChVector<> pos = bucket_ctr + ChVector<>(2.5*bucket_rad, 0, bucket_interior_halfDim.z);

	auto xaxis = std::make_shared<ChBody>();
	auto yaxis = std::make_shared<ChBody>();
	auto zaxis = std::make_shared<ChBody>();


	xaxis->SetPos(pos);						yaxis->SetPos(pos);							zaxis->SetPos(pos);
	xaxis->SetCollide(false);			yaxis->SetCollide(false);				zaxis->SetCollide(false);
	xaxis->SetBodyFixed(true);		yaxis->SetBodyFixed(true);			zaxis->SetBodyFixed(true);
	xaxis->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	yaxis->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	zaxis->GetCollisionModel()->SetEnvelope(collisionEnvelope);

	utils::AddCylinderGeometry(xaxis.get(), rad, len, ChVector<>(len - rad / 2, 0, 0) + pos, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, PI_2)), true);//bottom
	utils::AddCylinderGeometry(yaxis.get(), rad, len, ChVector<>(0, len-rad/2, 0) + pos, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, 0)), true);//bottom, true);//bottom
	utils::AddCylinderGeometry(zaxis.get(), rad, len, ChVector<>(0, 0, len - rad) + pos, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(PI_2, 0, 0)), true);//bottom

	xaxis->AddAsset(std::make_shared<ChColorAsset>(1.0f, 0, 0));
	yaxis->AddAsset(std::make_shared<ChColorAsset>(0, 1.0f, 0));
	zaxis->AddAsset(std::make_shared<ChColorAsset>(0, 0, 1.0f));
	mphysicalSystem.AddBody(xaxis); mphysicalSystem.AddBody(yaxis); mphysicalSystem.AddBody(zaxis);

}



void recycleSmarticles(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mySmarticlesVec)
{
	double pos = -.75*bucket_interior_halfDim.z;//z position below which smarticles are regenerated above pile inside container
	double ang = 2 * PI / numPerLayer;
	double rp = genRand(ang/4) ; //add slight offset to angInc to allow particles not always fall in nearly same position
	static int recycledSmarticles = 0;
	static int inc = 0;
	for (size_t i = 0; i < mySmarticlesVec.size(); i++)
	{
		Smarticle* sPtr = mySmarticlesVec[i];
		if (sPtr->GetArm(1)->GetPos().z < pos)
		{


			if (bucketType == HOPPER)
			{
				sPtr->TransportSmarticle(bucket_ctr + ChVector<>(
					sin(ang*inc + rp)*(bucket_rad / 2 + 4*w_smarticle*(genRand(-0.5,0.5))),
					cos(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(genRand(-0.5, 0.5))),
					bucket_interior_halfDim.z*2
					));

				//sPtr->SetSpeed(sPtr->GetArm(1)->GetPos_dt() / 4);
				sPtr->SetSpeed(ChVector<>(0, 0, -9.8*.01 / 2.0 - w_smarticle / .01));
			}
			else
			{

				sPtr->TransportSmarticle(bucket_ctr + ChVector<>(
					sin(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(genRand(-0.5, 0.5))),
					cos(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(genRand(-0.5, 0.5))),
					bucket_interior_halfDim.z*1.75
					));
			}


			recycledSmarticles++;
			inc = (inc+1)%numPerLayer;
		}
	}
	printFlowRate(mphysicalSystem.GetChTime(), recycledSmarticles);
}
// =============================================================================
void FixBodies(CH_SYSTEM& mphysicalSystem, int tStep) {
	//std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
	std::vector<std::shared_ptr<ChBody> >::iterator ibody = mphysicalSystem.Get_bodylist()->begin();

	for (size_t i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
		//ChBody* bodyPtr = *(myIter + i);
		auto bodyPtr = *(ibody + i);
		if (bodyPtr->GetPos().z < -5.0*bucket_interior_halfDim.z) {
			bodyPtr->SetBodyFixed(true);
		}
	}
}
// =============================================================================

void FixRotation(CH_SYSTEM& mphysicalSystem, Smarticle* sPtr) //reduces rotation speed by half 
{
	if (sPtr->GetArm(0)->GetRot_dt().GetVector().Length2() > 10000)//added this because small arms can start to spin uncontrollably
	{
		sPtr->GetArm(0)->SetRot_dt(sPtr->GetArm(0)->GetRot() / 2);
		sPtr->GetArm(1)->SetRot_dt(sPtr->GetArm(1)->GetRot() / 2);
		sPtr->GetArm(2)->SetRot_dt(sPtr->GetArm(2)->GetRot() / 2);
		GetLog() << "\n\n******WARNING******\n arm is rotating too fast and FixRotation method is running\n******WARNING******\n";
	}
	
}
void EraseSmarticle(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*>::iterator& myIter, Smarticle& sPtr, std::vector<Smarticle*> &mySmarticlesVec)
{
	sPtr.~Smarticle();
	myIter = mySmarticlesVec.erase(myIter); 

	//sPtr->~Smarticle();
	//myIter = mySmarticlesVec.erase(myIter);
}
void FixSmarticles(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mySmarticlesVec, double tstep) { ///remove all traces of smarticle from system //TODO REMAKE METHOD
	if (bucketType == HOPPER && bucket_exist == false) //if hopper, put smarticles back inside after reaching below hopper if bucket_bott still exists delete
	{
		recycleSmarticles(mphysicalSystem, mySmarticlesVec);
	}

	std::vector<Smarticle*>::iterator myIter;
	for (myIter = mySmarticlesVec.begin(); myIter != mySmarticlesVec.end();)
	{
		Smarticle* sPtr = *(myIter);
		//if (sPtr->armBroken)
		//{
		//	EraseSmarticle(mphysicalSystem, myIter, *sPtr, mySmarticlesVec);
		//	GetLog() << "\nArm broken removing smarticle \n";
		//	continue;
		//}
		FixRotation(mphysicalSystem, sPtr);

		//if smarticles are too low and not hopper
		if (bucketType != HOPPER)
		{
			if (sPtr->GetArm(1)->GetPos().z < -3.0*bucket_interior_halfDim.z)
			{
				EraseSmarticle(mphysicalSystem, myIter, *sPtr, mySmarticlesVec);
				GetLog() << "\nRemoving Smarticle far below container \n";
				continue;
			}
		}

		switch (bucketType)
		{
		case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
			if (!IsInRadial(sPtr->Get_cm(), bucket_bott->GetPos() + ChVector<>(0, 0, bucket_interior_halfDim.z), ChVector<>(bucket_rad*3, bucket_bott->GetPos().z, bucket_bott->GetPos().z + 4 * bucket_interior_halfDim.z)))
			{
				EraseSmarticle(mphysicalSystem, myIter, *sPtr, mySmarticlesVec);
				GetLog() << "\nRemoving smarticle outside system \n";
				continue;
			}
			++myIter;
			break;


		case HOPPER:
			if (!IsInRadial(sPtr->GetArm(1)->GetPos(), bucket->GetPos(), ChVector<>(2 * bucket_rad, -4.0*bucket_interior_halfDim.z, 4.0*bucket_interior_halfDim.z)))
			{
				EraseSmarticle(mphysicalSystem, myIter, *sPtr, mySmarticlesVec);
				GetLog() << "\nRemoving smarticle outside hopper \n";

			}
			else{ ++myIter; }
			break;

		default:
			++myIter;
			break;
		}
	}


}

void PrintStress(CH_SYSTEM* mphysicalSystem, int tstep, double zmax,double cylrad) //TODO include knobs in calculation
{

	const static std::string stress = out_dir + "/Stress.txt";
	std::ofstream stress_of;
	if (tstep == 0) {
		stress_of.open(stress.c_str());
	}
	else {
		stress_of.open(stress.c_str(), std::ios::app);	}
	ChVector<> temp = bucket_bod_vec.at(1)->GetPos();
	double currBuckRad = sqrt(temp.x*temp.x + temp.y*temp.y) - bucket_half_thick / 5.0;//bucket_half_thick/5 is how wall thickness is defined!
	//GetLog() << bucket_half_thick<< "thick\n";
	//showForce(mphysicalSystem)/(PI*2*cylrad*zmax)
	double force = showForce(mphysicalSystem);
	//GetLog() << "\nforce:" << force;
	stress_of << mphysicalSystem->GetChTime() << ", " << force <<","<< Smarticle::global_GUI_value <<", "<< currBuckRad<< std::endl;
	stress_of.close();
}
void PrintStress2(CH_SYSTEM* mphysicalSystem, int tstep, double zmax, double cylrad, std::vector<Smarticle*> mySmarticlesVec) //TODO include knobs in calculation
{

	bool printAllSmarticleInfo = true;
	static int frame = 0;

	const static std::string stress = out_dir + "/Stress.txt";
	std::ofstream stress_of;
	if (tstep == 0) {
		stress_of.open(stress.c_str());
		stress_of << dT << ", " <<out_fps<<", " <<videoFrameInterval<< ", " << bucket_rad << ", "<< bucketType <<std::endl;
	}
	else {
		stress_of.open(stress.c_str(), std::ios::app);
	}
	//GetLog() << bucket_half_thick<< "thick\n";
	//showForce(mphysicalSystem)/(PI*2*cylrad*zmax)
	ChVector<> temp;
	double currBuckRad;
	switch (bucketType)
	{
	case CYLINDER:
		temp = bucket_bod_vec.at(1)->GetPos();
		currBuckRad = sqrt(temp.x*temp.x + temp.y*temp.y) - bucket_half_thick / 5.0;//bucket_half_thick/5 is how wall thickness is defined!
		stress_of << mphysicalSystem->GetChTime() << ", " << 0 << ", " << Smarticle::global_GUI_value << ", " << currBuckRad << ", " << 0 << std::endl; //final 0 is a placeholder 
		break;
	case STRESSSTICK: case KNOBCYLINDER:
		temp = bucket_bod_vec.at(1)->GetPos();
		currBuckRad = sqrt(temp.x*temp.x + temp.y*temp.y) - bucket_half_thick / 5.0;//bucket_half_thick/5 is how wall thickness is defined!
		stress_of << mphysicalSystem->GetChTime() << ", " << showForce(mphysicalSystem) << ", " << Smarticle::global_GUI_value << ", " << currBuckRad << ", " << 0 << std::endl;
		break;
	case BOX:
		stress_of << mphysicalSystem->GetChTime() << ", " << 0 << ", " << Smarticle::global_GUI_value << ", " << box_ang << ", " << 0<< std::endl;
		break;
	}


	if (printAllSmarticleInfo)
	{
		for (size_t i = 0; i < mySmarticlesVec.size(); i++)
		{
			stress_of << mySmarticlesVec[i]->GetAngle1(true) << ", " << mySmarticlesVec[i]->GetAngle2(true) << ", " << mySmarticlesVec[i]->moveType << ", " << mySmarticlesVec[i]->Get_cm().z << std::endl;
		}
		stress_of << "#EF" <<frame<< std::endl;
	}
	stress_of.close();
	frame = frame + 1;
}
void PrintFractions(CH_SYSTEM& mphysicalSystem, int tStep, std::vector<Smarticle*> mySmarticlesVec) {

	const static std::string vol_frac = out_dir + "/volumeFraction.txt";
	static int stepSave = 10;
	if (tStep % stepSave != 0) return;
	double zComz = 0;
	double meanOT = 0;
	//static std::shared_ptr<ChBody> grid;  //uncomment to visualize vol frac boxes
	//static bool a = false;						//uncomment to visualize vol frac boxes
	std::ofstream vol_frac_of;
	if (tStep == 0) {
	  vol_frac_of.open(vol_frac.c_str());
	} else {
	  vol_frac_of.open(vol_frac.c_str(), std::ios::app);
	}

	//double sqSize = w_smarticle; // try increasing!
	//int rowSize = std::ceil(bucket_rad*2/sqSize);
	double sqSizex = (w_smarticle+2*t2_smarticle); // try increasing!
	double sqSizey = (l_smarticle+2*t2_smarticle); // try increasing!
	int colSize = std::ceil(bucket_rad*2/sqSizex);
	int rowSize = std::ceil(bucket_rad*2/sqSizey);
	std::pair<int,double> p(0,0.0);
	// std::vector<std::pair<int,double>> zHeights(rowSize*rowSize,p);
	 std::vector<std::pair<int,double>> zHeights(rowSize*colSize,p);

	double zmax = 0;
	double max2 = 0;
	int xpos = 0;
	int ypos = 0;
	int vecPos=0;
	ChVector<> com;
		ChVector<> pos;
	double zMax =0;
	int countInside =0;

	//double zMax = Find_Max_Z(mphysicalSystem,mySmarticlesVec);
	ChVector<> bucketMin = bucket_bott->GetPos();


	// *** remember, 2 * bucket_half_thick is needed since bucket is initialized inclusive. the half dims are extended 2*bucket_half_thick from each side
	ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, bucket_interior_halfDim.z);
	double totalVolume2 = 0;
	int countInside2 = 0;
	double volumeFraction = 0;
	

	switch (bucketType)
	{
	case BOX:
		zMax = Find_Max_Z(mphysicalSystem, mySmarticlesVec);
		zMax = std::min(zMax, bucketMin.z + 2 * bucket_interior_halfDim.z);
		for (size_t i = 0; i < mySmarticlesVec.size(); i++) {
			Smarticle* sPtr = mySmarticlesVec[i];
			if (IsIn(sPtr->Get_cm(), bucketCtr - bucket_interior_halfDim, bucketCtr + bucket_interior_halfDim + ChVector<>(0, 0, 2.0 * bucket_half_thick))) {
				countInside2++;
				totalVolume2 += sPtr->GetVolume();
			}
		}

		volumeFraction = totalVolume2 / (4.0 * bucket_interior_halfDim.x * bucket_interior_halfDim.y * (zMax - bucketMin.z));
		break;
	case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
		for (size_t i = 0; i < mySmarticlesVec.size(); i++) {
			Smarticle* sPtr = mySmarticlesVec[i];
			//isinradial rad parameter is Vector(bucketrad,zmin,zmax)
			if (IsInRadial(sPtr->Get_cm(), bucketCtr, ChVector<>(bucket_rad, bucketMin.z, bucketMin.z + 2.0*bucket_interior_halfDim.z))) {
				countInside2++;
				//com = sPtr->Get_cm()-ChVector<>(0,0,bucketMin.z);
				com = sPtr->Get_cm() - ChVector<>(0, 0, bucketMin.z);
				zComz += com.z;
				max2 = std::max(max2, com.z);
				if (max2>zMax)
				{
					double temp = zMax;
					zMax = max2;
					max2 = temp;
				}
				meanOT += sPtr->GetReactTorqueLen01() + sPtr->GetReactTorqueLen12();
				//zMax = std::max(zMax, sPtr->GetArm(1)->GetPos().z- bucketMin.z);

			}
		}
		volumeFraction = countInside2*vol / (max2*PI*bucket_rad*bucket_rad);
		//GetLog() << vol << " " << countInside2 << " " << bucket_rad << " " << zMax << " " << volumeFraction << "\n";
		//GetLog() << "phi=" << volumeFraction << "\n";
		zComz = zComz / countInside2;
		meanOT = meanOT / (countInside2 * 2.0); //multiply by 2 (2 arms for each smarticle)
		break;
	default:
		break;
	}
	vol_frac_of << mphysicalSystem.GetChTime() << ", " << countInside2 << ", " << volumeFraction << ", " << zMax << ", " << zComz << ", " << meanOT<< ", " << Smarticle::global_GUI_value << std::endl;
	vol_frac_of.close();
	return;
}

// move bucket
void vibrate_bucket(double t) {
		double phase = -omega_bucket*vibrateStart;
		double x_bucket = vibration_amp*sin(omega_bucket * t + phase);
		double xDot_bucket = vibration_amp*omega_bucket*cos(omega_bucket * t + phase);
		double xDDot_bucket = vibration_amp*omega_bucket*omega_bucket*-1 * sin(omega_bucket * t + phase);
		bucket->SetPos(ChVector<>(0, 0, x_bucket));
		bucket->SetPos_dt(ChVector<>(0, 0, xDot_bucket));
		bucket->SetPos_dtdt(ChVector<>(0, 0, xDDot_bucket));
		bucket->SetRot(QUNIT);
		if (bucketType == CYLINDER || bucketType == STRESSSTICK || bucketType == HOOKRAISE || bucketType==KNOBCYLINDER)
	{
		bucket_bott->SetPos(ChVector<>(0, 0, x_bucket));
		bucket_bott->SetPos_dt(ChVector<>(0, 0, xDot_bucket));
		bucket_bott->SetPos_dtdt(ChVector<>(0, 0, xDDot_bucket));
		bucket_bott->SetRot(QUNIT);
	}
}
void vibrate_bucket(double t, std::shared_ptr<ChBody> body) {
		body->SetBodyFixed(false);
		double phase = -omega_bucket*vibrateStart;
		double x_bucket = vibration_amp*sin(omega_bucket * t + phase);
		double xDot_bucket = vibration_amp*omega_bucket*cos(omega_bucket * t + phase);
		double xDDot_bucket = vibration_amp*omega_bucket*omega_bucket*-1 * sin(omega_bucket * t + phase);
		//body->SetPos(ChVector<>(0, 0, x_bucket));
		body->SetPos_dt(ChVector<>(0, 0, xDot_bucket));
		body->SetPos_dtdt(ChVector<>(0, 0, xDDot_bucket));
		body->SetRot(QUNIT);
		bucket->SetPos(ChVector<>(0, 0, x_bucket));
}
void rotate_bucket(double t)//method is called on each iteration to rotate drum at an angular velocity of drum_omega
{
	static std::shared_ptr<ChFunction_Const> mfun2;
	if (bucketType == DRUM)
	{
		bucket->SetBodyFixed(false);
		bucket_bott->SetBodyFixed(true);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_spe_funct());
		drum_omega = drum_freq*PI * 2;
		mfun2->Set_yconst(drum_omega);
	}
	else if (bucketType == BOX)
	{
		bucket->SetBodyFixed(false);
		bucket_bott->SetBodyFixed(true);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
		//set rotation amount equal to box_ang-currentangle, therefore no rotation if 2 values are equal
		mfun2->Set_yconst(box_ang-Quat_to_Angle(ANGLESET_RXYZ,bucket->GetRot()).x);
	}
}
//set up actuat
void setUpBucketActuator(CH_SYSTEM& mphysicalSystem)
{
	std::shared_ptr<ChFunction_Const> mfun2; //needs to be declared outside switch
	ChVector<> pR01(0, 0, 0);
	ChQuaternion<> qx;
	bucket_actuator = std::make_shared<ChLinkEngine>();
	switch (bucketType)
	{
	case DRUM:
		
		
		qx = Q_from_AngAxis(PI_2, VECT_Z);
		bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), bucket->GetRot()));
		bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
		mphysicalSystem.AddLink(bucket_actuator);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_spe_funct());
		drum_omega = drum_freq*PI * 2;
		mfun2->Set_yconst(drum_omega);
		break;

	case BOX:
		qx = Q_from_AngAxis(PI_2, VECT_Y);
		bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), qx));
		bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
		//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
		mphysicalSystem.AddLink(bucket_actuator);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
		mfun2->Set_yconst(0);
		break;
	}


}

// =============================================================================

	//CImage image;
	//image.Attach(hBitmap);
	//image.Save("c:\\pngPicture.png");

	//hres = CreateStreamOnHGlobal(0, TRUE, &pStream);
	//hr = myImage.Save(pStream, Gdiplus::ImageFormatPNG);
void UpdateSmarticles(
		CH_SYSTEM& mphysicalSystem,
		std::vector<Smarticle*> mySmarticlesVec) {

	//static double torquethresh = mySmarticlesVec[0]->torqueThresh2;

	double t = mphysicalSystem.GetChTime(); 
	
	for (size_t i = 0; i < mySmarticlesVec.size(); i++) {
		mySmarticlesVec[i]->updateTorqueDeque();
		double tor1 = std::get<0>(mySmarticlesVec[i]->torqueAvg);
		double tor2 = std::get<1>(mySmarticlesVec[i]->torqueAvg);
		

		
		int moveType=0;
		///////////////////random chance at current timestep for smarticle to not move to globalValue, models real life delay for smarticles to start motion to current state
		if (genRand() < percentToMoveToGlobal)
			moveType = Smarticle::global_GUI_value; 
		else
			moveType = mySmarticlesVec[i]->prevMoveType;

		mySmarticlesVec[i]->ControllerMove(moveType, tor1, tor2);

		//GetLog() << "\nangle(1,2):" << mySmarticlesVec[i]->GetAngle1(true) << "\t" << mySmarticlesVec[i]->GetAngle2(true);
		//GetLog() << "\n\ntor(1,2):" <<tor1 << "\t" << tor2;
		mySmarticlesVec[i]->steps = mySmarticlesVec[i]->steps + 1;
	}
}
// =============================================================================
bool SetGait(double time)
{

	if (time <= 3)
		Smarticle::global_GUI_value = 2;
	else if (time > 3 && time <= 6)
		Smarticle::global_GUI_value = 3;
	else if (time > 6 && time <= 9)
		Smarticle::global_GUI_value = 2;
	else if (time > 9 && time <= 12)
		Smarticle::global_GUI_value = 3;
	else if (time > 12 && time <= 15)
		Smarticle::global_GUI_value = 2;
	else if (time > 15 && time <= 18)
		Smarticle::global_GUI_value = 3;
	else if (time > 18 && time <= 21)
		Smarticle::global_GUI_value = 2;
	else if (time > 21 && time <= 24)
		Smarticle::global_GUI_value = 3;
	else if (time > 24 && time <= 27)
		Smarticle::global_GUI_value = 2;
	else if (time > 27 && time <= 30)
		Smarticle::global_GUI_value = 3;
	//else if (time > 30 && time <= 33)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 33 && time <= 36)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 36 && time <= 39)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 39 && time <= 42)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 42 && time <= 45)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 45 && time <= 48)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 48 && time <= 51)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 51 && time <= 54)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 54 && time <= 57)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 57 && time <= 60)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 60 && time <= 63)
	//	Smarticle::global_GUI_value = 3;
	//else	
	//	return true;

	//if (time <= 5)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 5 && time <= 10)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 10 && time <= 15)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 15 && time <= 20)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 20 && time <= 25)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 25 && time <= 30)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 30 && time <= 35)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 35 && time <= 40)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 40 && time <= 45)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > 45 && time <= 50)
	//	Smarticle::global_GUI_value = 3;
	//else if (time > 50 && time <= 55)
	//	Smarticle::global_GUI_value = 2;
	//else	
	//	return true;

	return false;


	// false;
}

int main(int argc, char* argv[]) {
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	//ChTimerParallel step_timer;
	Smarticle::global_GUI_value = 1;
	
	//set chrono dataPath to data folder placed in smarticle directory so we can share created files
#if defined(_WIN64)
	char* pPath = getenv("USERNAME");
	GetLog()<<pPath;
	std::string fp;
	if (strcmp(pPath, "root") == 0)
		fp = std::string("D:\\ChronoCode\\chronoPkgs\\Smarticles\\data\\");
	else
		fp = std::string("D:\\GT Coursework\\smarticles\\data\\");
	//fp = __FILE__+fp;
	
	SetChronoDataPath(fp);
#else
	SetChronoDataPath("/home/wsavoie/Documents/ChronoPkgs/Smarticles/data/");
#endif

	time_t rawtimeCurrent;
	struct tm* timeinfoDiff;
	// --------------------------
	// Create output directories.
	// --------------------------


	if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
		std::cout << "Error creating directory " << out_dir << std::endl;
		return 1;
	}

	if (povray_output) {
		if (ChFileutils::MakeDirectory(pov_dir_mbd.c_str()) < 0) {
			std::cout << "Error creating directory " << pov_dir_mbd << std::endl;
			return 1;
		}
	}

	const std::string rmCmd = std::string("rm -rf ") + pov_dir_mbd;
	std::system(rmCmd.c_str());

	if (povray_output) {
		if (ChFileutils::MakeDirectory(pov_dir_mbd.c_str()) < 0) {
			std::cout << "Error creating directory " << pov_dir_mbd << std::endl;
			return 1;
		}
	}
	sphereTexture->SetTextureFilename(GetChronoDataFile("sphereTexture.png"));
	const std::string simulationParams = out_dir + "/simulation_specific_parameters.txt";
	simParams.open(simulationParams.c_str());
	simParams << "Job was submitted at date/time: " << asctime(timeinfo) << std::endl;
	simParams.close();
	// define material property for everything
	mat_g = std::make_shared<ChMaterialSurface>();
	mat_g->SetFriction(fric); // .6 for wall to staple using tan (theta) tested on 7/20

	// Create a ChronoENGINE physical system
	CH_SYSTEM mphysicalSystem;


	InitializeMbdPhysicalSystem_NonParallel(mphysicalSystem, argc, argv);
	GetLog() << "\npctActive" << pctActive << "\n";
	Smarticle::pctActive = pctActive;
	MyBroadPhaseCallback mySmarticleBroadphaseCollisionCallback;
	mphysicalSystem.GetCollisionSystem()->SetBroadPhaseCallback(&mySmarticleBroadphaseCollisionCallback);


	std::vector<Smarticle*> mySmarticlesVec;
	CreateMbdPhysicalSystemObjects(mphysicalSystem, mySmarticlesVec);


#ifdef CHRONO_OPENGL
	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	//	ChVector<> CameraLocation = ChVector<>(0, -10, 4);
	//	ChVector<> CameraLookAt = ChVector<>(0, 0, -1);
	ChVector<> CameraLocation = sizeScale * ChVector<>(-.1, -.06, .1);
	ChVector<> CameraLookAt = sizeScale * ChVector<>(0, 0, -.01);
	gl_window.Initialize(appWidth, appHeight, "Dynamic Smarticles", &mphysicalSystem);
	gl_window.SetCamera(CameraLocation, CameraLookAt, ChVector<>(0, 0, 1)); //camera
	gl_window.viewer->render_camera.camera_scale = 2.0 / (1000.0)*sizeScale;
	gl_window.viewer->render_camera.near_clip = .001;
	gl_window.SetRenderMode(opengl::WIREFRAME);

	// Uncomment the following two lines for the OpenGL manager to automatically
	// return 0;
#endif

#if irrlichtVisualization
	std::cout << "@@@@@@@@@@@@@@@@  irrlicht stuff  @@@@@@@@@@@@@@@@" << std::endl;
	// Create the Irrlicht visualization (open the Irrlicht device,
	// bind a simple user interface, etc. etc.)
	ChIrrApp application(&mphysicalSystem, L"Dynamic Smarticles",
		core::dimension2d<u32>(appWidth, appHeight), false, true);
	////////////!@#$%^
	// Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
	ChIrrWizard::add_typical_Logo(application.GetDevice());
	ChIrrWizard::add_typical_Sky(application.GetDevice());
	ChIrrWizard::add_typical_Lights(application.GetDevice(),
		core::vector3df(0, 0, 4*bucket_rad)*sizeScale,
		core::vector3df(0, 0, 1*bucket_rad)*sizeScale);

	ChIrrWizard::add_typical_Lights(application.GetDevice(),
		core::vector3df(.0139, -.39, -.0281)*sizeScale,
		core::vector3df(0.0139, .298, -.195)*sizeScale);

	ChIrrWizard::add_typical_Lights(application.GetDevice());
	

	RTSCamera* camera = new RTSCamera(application.GetDevice(), application.GetDevice()->getSceneManager()->getRootSceneNode(),
		application.GetDevice()->getSceneManager(), -1, -50.0f, 0.5f, 0.0005f);
	camera->setUpVector(core::vector3df(0, 0, 1));
//camera->setPosition(core::vector3df(0, -.17, .07));
	//camera->setPosition(core::vector3df(0, bucket_interior_halfDim.y, 4 * bucket_rad / 2.0));
	//camera->setTarget(core::vector3df(0, 2*bucket_interior_halfDim.y, -bucket_rad / 2.1)); //	camera->setTarget(core::vector3df(0, 0, .01));

	camera->setPosition(core::vector3df(.0139, -.445, .7));
	camera->setTarget(core::vector3df(0.0139, .243, -.195)); //	camera->setTarget(core::vector3df(0, 0, .01));



	ChIrrWizard::add_typical_Lights(application.GetDevice(),
		camera->getPosition(),
		camera->getTarget());

	camera->setNearValue(0.0005f);
	camera->setMinZoom(0.1f);
	camera->setZoomSpeed(0.1f);
	
	drawGlobalCoordinateFrame(mphysicalSystem);

	//framerecord
	application.SetVideoframeSaveInterval(videoFrameInterval);//only save out frames to make it 30fps




	// Use this function for adding a ChIrrNodeAsset to all items
	// If you need a finer control on which item really needs a visualization
	// proxy in
	// Irrlicht, just use application.AssetBind(myitem); on a per-item basis.
	application.AssetBindAll();
	// Use this function for 'converting' into Irrlicht meshes the assets
	// into Irrlicht-visualizable meshes
	application.AssetUpdateAll();
	application.SetStepManage(true);
	application.SetTimestep(dT);  // Arman modify

	std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << std::endl;
	IrrGui receiver(&application, &mySmarticlesVec);
	// scan all contact
	// note how to add the custom event receiver to the default interface:
	application.SetUserEventReceiver(&receiver);



	receiver.drawAngle();//initialize draw angle
	application.SetShowInfos(true);

#endif


	int stepEnd = int(tFinal / dT);  // 1.0e6;//2.4e6;//600000;//2.4e6 * (.02 * paramsH.sizeScale) /

	auto truss = std::make_shared<ChBody>();
	auto link_engine = std::make_shared<ChLinkEngine>();
	auto sinefunc = std::make_shared<ChFunction_Sine>();
	auto func2 = std::make_shared<ChFunction_Const>();
	auto knobcylinderfunc = std::make_shared<ChFunction_Const>();
	auto func = std::make_shared<ChFunctionCustom>();
	auto knobstick = std::make_shared<ChBody>();
	auto stick = std::make_shared<ChBody>();

	//std::shared_ptr<ChLinkLinActuator> pris_engine(new ChLinkLinActuator);
	std::shared_ptr<ChLinkLinActuator> pris_engine;
	std::shared_ptr<ChLinkLockPrismatic> link_prismatic;
	double rad = 0;


	switch (bucketType)
	{
	case STRESSSTICK: case HOOKRAISE:// case KNOBCYLINDER:
	{

		stick->SetRot(QUNIT);
		stick->SetBodyFixed(true);
		stick->SetMaterialSurface(mat_g);
		stick->AddAsset(groundTexture);
		stick->GetCollisionModel()->ClearModel();
		stick->SetMass(4);


		double mult = 4.0;
		
		if (stapleSize)
		{
			rad = t_smarticle*mult/10.0;
		}
		else
		{
			rad = t_smarticle*mult / 4.0;
		}

		double stickLen = bucket_interior_halfDim.z*1.5;

		int sphereNum = stickLen / (t_smarticle / 2);
		if (bucketType==STRESSSTICK)
			sphereNum = stickLen / (t_smarticle / 2);;

		//double sphereStickHeight = t_smarticle*mult / 2.0 * (sphereNum + 1); //shouldnt need extra 2*rad offset because of how z is defined using i below
		for (size_t i = 0; i < sphereNum; i++)
		{
			
			stick->GetCollisionModel()->SetEnvelope(collisionEnvelope);
			//utils::AddSphereGeometry(stick.get(), t_smarticle / 2, bucket_ctr + ChVector<>(0, 0, t_smarticle*(i + 1 / 2.0)), QUNIT, true); // upper part, min_x plate
			//utils::AddSphereGeometry(stick.get(), t_smarticle / 5, bucket_ctr + ChVector<>(0, 0, t_smarticle*(i + 1 / 5.0)), QUNIT, true); // upper part, min_x plate
			//utils::AddSphereGeometry(stick.get(), t2_smarticle/2.0, bucket_ctr + ChVector<>(0, 0, t2_smarticle*(i + 1 /2.0)), QUNIT, true); // upper part, min_x plate
		
			//if you change z height between spheres, you must change sphereStickHeight above!
			utils::AddSphereGeometry(stick.get(), rad, bucket_ctr + ChVector<>(0, 0, stickLen / sphereNum * (i)), Angle_to_Quat(ANGLESET_RXYZ, ChVector<double>(0, 0, PI)), true);
			sphereStick.emplace_back(stick);
		
			
		}
		if (bucketType == HOOKRAISE)
		{
			int hookNum;
			hookNum = 8 - stapleSize * 2;

			for (size_t i = 0; i < hookNum; i++)
			{
				if (stapleSize)
				{

					//AddBoxGeometry
					utils::AddSphereGeometry(stick.get(), rad, bucket_ctr + ChVector<>(rad*(i + 1)*2, 0, stickLen / sphereNum), Angle_to_Quat(ANGLESET_RXYZ, ChVector<double>(PI, 0, 0)), true);
					
				}
				else
				{
					utils::AddSphereGeometry(stick.get(), rad, bucket_ctr + ChVector<>(rad*(i + 1)*2, 0, stickLen / sphereNum), Angle_to_Quat(ANGLESET_RXYZ, ChVector<double>(PI, 0, 0)), true);
				}
				sphereStick.emplace_back(stick);

				


				
			}
		}

		stick->GetCollisionModel()->BuildModel();
		stick->GetCollisionModel()->SetFamily(1);
		stick->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
		stick->GetPhysicsItem()->SetIdentifier(largeID + 1);
		stick->SetCollide(true);
		mphysicalSystem.AddBody(stick);


		truss->SetBodyFixed(true);
		truss->GetCollisionModel()->ClearModel();
		utils::AddCylinderGeometry(truss.get(), t2_smarticle / 2, bucket_interior_halfDim.z * 1, bucket_ctr + ChVector<>(0, 0, bucket_interior_halfDim.z), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(-PI_2, 0, 0)), false);
		truss->GetCollisionModel()->BuildModel();
		truss->AddAsset(sphereTexture);
		truss->SetCollide(false);
		mphysicalSystem.AddBody(truss);

		link_prismatic = std::make_shared<ChLinkLockPrismatic>();
		link_prismatic->Initialize(stick, truss, true, ChCoordsys<>(), ChCoordsys<>(ChVector<>(0, 0, 0), QUNIT));  // set prism as vertical (default would be aligned to z, horizontal
		mphysicalSystem.AddLink(link_prismatic);
		

		pris_engine = std::make_shared<ChLinkLinActuator>();
		pris_engine->Initialize(stick, truss, true, ChCoordsys<>(stick->GetPos() + ChVector<>(0, 0, -stickLen), QUNIT), ChCoordsys<>(stick->GetPos() + ChVector<>(0, 0, stickLen), QUNIT));
		
		


		GetLog() << "StickLen:" << stickLen;

		//func = std::dynamic_pointer_cast<ChFunctionCustom>(pris_engine->Get_dist_funct());
		func->Set_y(0);
		func->Set_y_dx(2.5-.5); //the value in this is always -2.5+(value specified), dont know where -2.5 comes from....
		pris_engine->Set_dist_funct(func);


		pris_engine->SetDisabled(true);
		mphysicalSystem.AddLink(pris_engine);

		break;
	}

	case DRUM: case BOX:
	{
		setUpBucketActuator(mphysicalSystem);
		break;
	}

	
	case HOPPER:
	{
		truss->SetBodyFixed(true);
		truss->GetCollisionModel()->ClearModel();
		utils::AddCylinderGeometry(truss.get(), t2_smarticle / 2, bucket_interior_halfDim.z * 1, bucket_ctr + ChVector<>(0, 0, bucket_interior_halfDim.z), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(-PI_2, 0, 0)), false);
		truss->GetCollisionModel()->BuildModel();
		truss->AddAsset(sphereTexture);
		truss->SetCollide(false);
		mphysicalSystem.AddBody(truss);


		link_prismatic = std::make_shared<ChLinkLockPrismatic>();
		link_prismatic->Initialize(bucket, truss, false, ChCoordsys<>(), ChCoordsys<>(ChVector<>(0, 0, 0), QUNIT));  // set prism as vertical (default would be aligned to z, horizontal
		mphysicalSystem.AddLink(link_prismatic);
		pris_engine = std::make_shared<ChLinkLinActuator>();
		pris_engine->Initialize(bucket, truss, false, ChCoordsys<>(VNULL, QUNIT), ChCoordsys<>(VNULL, QUNIT));
		sinefunc->Set_amp(vibAmp);
		sinefunc->Set_w(omega_bucket);
		sinefunc->Set_phase(-omega_bucket*vibrateStart);
		pris_engine->Set_dist_funct(sinefunc);

		pris_engine->SetDisabled(true);
		mphysicalSystem.AddLink(pris_engine);
		break;

	}
	case KNOBCYLINDER:
	{		
		truss->SetBodyFixed(true);
		truss->GetCollisionModel()->ClearModel();
		utils::AddCylinderGeometry(truss.get(), t2_smarticle/2, bucket_interior_halfDim.z*1, bucket_ctr+ChVector<>(0,0,bucket_interior_halfDim.z), Angle_to_Quat(ANGLESET_RXYZ,ChVector<>(PI_2,0,0)), true);
		truss->GetCollisionModel()->BuildModel();
		truss->AddAsset(sphereTexture);
		truss->SetCollide(false);
		mphysicalSystem.AddBody(truss);


		double mult = 2.0;
		double knobRad;
		double stickLen = bucket_interior_halfDim.z*1.5;
		int sphereNum = stickLen / (t_smarticle/2);

		//double sphereStickHeight = t_smarticle*mult / 2.0 * (sphereNum + 1); //shouldnt need extra 2*rad offset because of how z is defined using i below
		double sphereStickHeight = 2*bucket_interior_halfDim.z; //shouldnt need extra 2*rad offset because of how z is defined using i below
		
		knobstick->GetCollisionModel()->ClearModel();
		knobstick->SetRot(QUNIT);
		knobstick->SetBodyFixed(true);
		knobstick->SetMaterialSurface(mat_g);
		knobstick->AddAsset(groundTexture);
		knobstick->SetMass(1);

		for (size_t i = 0; i < sphereNum; i++)
		{
			knobstick->GetCollisionModel()->SetEnvelope(collisionEnvelope);
			
			if (stapleSize)
			{
				rad = t_smarticle*mult*1.5;
				knobRad = t_smarticle*2;
			}
			else
			{
				rad = t_smarticle*mult / 2;
				knobRad = t_smarticle;

			}
			//if you change z height between spheres, you must change sphereStickHeight above!
			utils::AddSphereGeometry(knobstick.get(), rad, bucket_ctr + ChVector<>(0, 0, sphereStickHeight / sphereNum * (i)),QUNIT, true);
			sphereStick.emplace_back(knobstick);
		}
			unsigned int kpr = 4;//knobs per row
			unsigned int rows = 15; //knob per z
			double ang = 2 * PI / kpr;
			double hp = (sphereStickHeight - 2 * rad) / rows;//height between rows
			double pOffset = PI/kpr; //phase offset
			for (size_t row = 0; row < rows; row++)
			{
				for (size_t col = 0; col < kpr; col++)
				{
					double theta = col*ang +row*pOffset;
					//utils::AddSphereGeometry(knobstick.get(), knobRad, bucket_ctr + ChVector<>(rad*cos(col*ang + row*pOffset), rad*sin(col*ang + row*pOffset), hp*(row + 1)), Angle_to_Quat(col*ang + row*pOffset, VECT_Y), true);
					utils::AddBoxGeometry(knobstick.get(), ChVector<>(knobRad*1.5, rad / 4, rad / 8), bucket_ctr + ChVector<>(rad*cos(theta), rad*sin(theta), hp*(row + 1)), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, theta + row % 2 * PI_2)), true);
					sphereStick.emplace_back(knobstick);
				}
			}

		knobstick->GetCollisionModel()->BuildModel();
		knobstick->SetCollide(true);
		knobstick->GetCollisionModel()->SetFamily(1);
		knobstick->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
		knobstick->GetPhysicsItem()->SetIdentifier(largeID + 1);
		mphysicalSystem.AddBody(knobstick);
	
		

		//auto link_engine = std::make_shared<ChLinkEngine>();

		double knobAmp = PI_2;
		double knobW = 0;//// rod rotating speed knobW = PI
		double knobPhase = -knobW*vibrateStart;
		//knobcylinderfunc->Set_amp(knobAmp);
		//knobcylinderfunc->Set_w(knobW);
		//knobcylinderfunc->Set_phase(knobPhase);



		link_engine->Initialize(knobstick, truss,
			ChCoordsys<>(ChVector<>(0, 0, 0), QUNIT));
		link_engine->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
		link_engine->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		//link_engine->Set_rot_funct(knobcylinderfunc);

		auto mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(link_engine->Get_spe_funct());
		mfun2->Set_yconst(knobW);
		link_engine->SetDisabled(true);
		mphysicalSystem.AddLink(link_engine);
		
		//auto link = std::make_shared<ChLinkLockRevolute>();
		//link->Initialize(knobstick, truss, ChCoordsys<>(VNULL,,ChVector<>(1,0,0)));
		//link->SetMotion_axis(ChVector<>(0, 0,1));
		//mphysicalSystem.AddLink(link);
		
		break;

	}
		default:
			break;
	}
	double timeForVerticalDisplacement = 0.015;
	if (bucketType == DRUM)
		timeForVerticalDisplacement = 0.095; // 1.5 for safety proximity .015
	int numGeneratedLayers = 0;


	if (read_from_file==1)
	{
		CheckPointSmarticlesDynamic_Read(mphysicalSystem, mySmarticlesVec,application);
		application.AssetBindAll();
		application.AssetUpdateAll();
		numGeneratedLayers = numLayers;
	}
//  for (int tStep = 0; tStep < 1; tStep++) {
	//START OF LOOP 
	application.DrawAll();
	
	for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
		double t = mphysicalSystem.GetChTime();
		if (read_from_file != 1)
		{
			if ((fmod(t, timeForVerticalDisplacement) < dT) &&
				(numGeneratedLayers < numLayers)){
#if irrlichtVisualization
				AddParticlesLayer1(mphysicalSystem, mySmarticlesVec, application, timeForVerticalDisplacement);
#else
				AddParticlesLayer1(mphysicalSystem, mySmarticlesVec);
#endif
				numGeneratedLayers++;
			}

		}

		///add method about system actuation
		if (bucketType == DRUM || bucketType == BOX)
		{
			rotate_bucket(t);
		}

		//vibration movement
		if (t > vibrateStart && t < vibrateStart + 3)
		{
			switch (bucketType)
			{
			case HOOKRAISE: case STRESSSTICK:
			{

				//if (pris_engine->IsDisabled())
				//{
				//	stick->SetBodyFixed(false);
				//	pris_engine->SetDisabled(false);
				//
				//}
				//pris_engine->GetDist_dt();
				////pris_engine->GetRelC_dt()
				////GetLog() << pris_engine->GetDist_dt() << "\n";
				break;
			}
			case KNOBCYLINDER:
			{

				if (link_engine->IsDisabled())
				{
					knobstick->SetBodyFixed(false);
					link_engine->SetDisabled(false);
				}

				break;
			}

			case CYLINDER:
			{
				//vibrate cylinder
				for (size_t i = 0; i < bucket_bod_vec.size(); i++)
				{
					vibrate_bucket(t, bucket_bod_vec.at(i));
				}
				vibrate_bucket(t, bucket_bott);
				break;
			}
			case HOPPER:
			{
				bucket->SetBodyFixed(false);
				pris_engine->SetDisabled(false);
				break;
			}
			default:
				break;
			}
		}


		if ( (fmod(t, timeForVerticalDisplacement) < dT)	&&		(mySmarticlesVec.size() < numPerLayer*numLayers) &&	(numGeneratedLayers == numLayers) )
			AddParticlesLayer1(mphysicalSystem, mySmarticlesVec, application, timeForVerticalDisplacement);
		//SavePovFilesMBD(mphysicalSystem, tStep);
		//step_timer.start("step time");

	#ifdef CHRONO_OPENGL
			if (gl_window.Active()) {
				gl_window.DoStepDynamics(dT);
				gl_window.Render();
			}
	#else
	#if irrlichtVisualization
			if (!(application.GetDevice()->run())) break;
			//if ((application.GetDevice()->isWindowMinimized())) break;
			application.GetVideoDriver()->beginScene(true, true,
				video::SColor(255, 140, 161, 192));			

			for (size_t i = 0; i < mySmarticlesVec.size(); ++i)
			{
				application.AssetUpdate(mySmarticlesVec[i]->GetArm(0));
				application.AssetUpdate(mySmarticlesVec[i]->GetArm(2));
			}
			
			application.DoStep();
			UpdateSmarticles(mphysicalSystem, mySmarticlesVec);

			receiver.drawSmarticleAmt(numGeneratedLayers);
			receiver.drawSuccessful();
			application.DrawAll();

			application.GetVideoDriver()->endScene();
			
	#else
			
			mphysicalSystem.DoStepDynamics(dT);
			UpdateSmarticles(mphysicalSystem, mySmarticlesVec);
	#endif
#endif

		if (SetGait(t) == true)
			break;


		if (bucketType == STRESSSTICK || bucketType == KNOBCYLINDER|| bucketType==CYLINDER|| bucketType==BOX)
		{
			double zmax = Find_Max_Z(mphysicalSystem, mySmarticlesVec);
			PrintStress2(&mphysicalSystem, tStep, zmax, rad, mySmarticlesVec);
		}
		
		FixSmarticles(mphysicalSystem, mySmarticlesVec, tStep);


		PrintFractions(mphysicalSystem, tStep, mySmarticlesVec);
		
	  time(&rawtimeCurrent);
	  double timeDiff = difftime(rawtimeCurrent, rawtime);
	  //step_timer.stop("step time");
		receiver.drawCamera();		
		receiver.dtPerFrame = videoFrameInterval;
		receiver.fps = out_fps;
		receiver.screenshot(receiver.dtPerFrame);


	  std::cout.flush();

		//CheckPointSmarticlesDynamic_Write(mySmarticlesVec,
	 // 		tStep,
	 // 		mat_g,
	 // 		l_smarticle,
	 // 		w_smarticle,
	 // 		t_smarticle,
	 // 		t2_smarticle,
	 // 		collisionEnvelope,
	 // 		rho_smarticle);

  }
	simParams.open(simulationParams.c_str(), std::ios::app);
	simParams << "Smarticle OT: " <<	 mySmarticlesVec.at(0)->torqueThresh2 << std::endl;
  for (int i = 0; i < mySmarticlesVec.size(); i++) {
	  delete mySmarticlesVec[i];

  }
	if (receiver.saveFrame)
	{
		receiver.SaveToMovie();
	}
  mySmarticlesVec.clear();

	simParams << "completed"<<std::endl;
  simParams.close();
  return 0;
}
