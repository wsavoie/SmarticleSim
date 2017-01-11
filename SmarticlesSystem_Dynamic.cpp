//https://www.lri.fr/~hansen/cmaes_inmatlab.html
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
//   Multibody dynamics engine
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
#include "SystemGeometry.h"
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
using namespace irrlicht;
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
//enum BucketType { KNOBCYLINDER, HOOKRAISE, STRESSSTICK, CYLINDER, BOX, HULL, RAMP, HOPPER, DRUM,FLATHOPPER};
SmarticleType smarticleType = SMART_ARMS;//SMART_U;
BucketType bucketType = BOX;
std::vector<std::shared_ptr<ChBody>> /*sphereStick*/;
//std::shared_ptr<ChBody> bucket;
//std::shared_ptr<ChBody> bucket_bott;
int overlaptest(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
double Find_Max_Z(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mSmartVec);
//double Find_Max_Z(CH_SYSTEM& mphysicalSystem);

std::ofstream flowRate_of;
std::ofstream simParams;
std::ofstream stress_of;
std::ofstream vol_frac_of;
std::ofstream ringPos;

double sizeScale = 1;
int appWidth = 1280;
int appHeight = 720;
bool saveFrame = false;
//double gravity = -9.81 * sizeScale;
double gravity = -9.81;

double vibrateStart=0.9;
double smart_fric = .2;//.3814; //keyboard box friction = .3814
double vibration_freq = 30;
double omega_bucket = 2 * PI * vibration_freq;  // 30 Hz vibration similar to Gravish 2012, PRL
double mGamma = 2.0 * gravity;
double vibration_amp = mGamma / (omega_bucket*omega_bucket);
unsigned int largeID = 10000000;
unsigned int smartIdCounter = 4; //start at non-zero value to not collide
//double dT = std::min(0.001, 1.0 / vibration_freq / 200);;//std::min(0.0005, 1.0 / vibration_freq / 200);
double dT = 0.0005;//std::min(0.0005, 1.0 / vibration_freq / 200);
double contact_recovery_speed = .5* sizeScale;
double tFinal = 100;
double ringRad = 0.192 / 2.0;



//double rho_smarticle = 7850.0 / (sizeScale * sizeScale * sizeScale);
//double rho_cylinder = 1180.0 / (sizeScale * sizeScale * sizeScale);
std::shared_ptr<ChMaterialSurface> mat_smarts;
//std::shared_ptr<ChMaterialSurface> mat_wall;
int numLayers = 100;
double armAngle = 90;
double sOmega = 5;  // smarticle omega



double gaitChangeLengthTime = .5;


////////////////rescaled robot geometry (3.93) based on w_smarticle scaling
////////////////robot dim is l/w =1, w=.046 t=.031 t2=.021
#if stapleSize
	
	double w_smarticle = sizeScale * 0.0117; // sizeScale * 0.0117
	double l_smarticle = 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	double t_smarticle = sizeScale * .00127;
	double t2_smarticle = sizeScale * .0005;
	double rho_smarticle = 7850.0;
	double rho_smarticleArm = 7850.0;
	double rho_smarticleMid = 7850.0;
	double p_gain = 0.025;   //.1//.2         //0.133
	double i_gain = 10;// 0.03;	 //.5//.225						//0.05
	double d_gain = 0.1; //.0025 //.01       //0.0033
	
#else
	
	double w_smarticle = sizeScale * 0.05316/ 1;
	double l_smarticle = 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	//real value
	double t_smarticle = sizeScale * .029982 / 1.0; //height of solar panels
	double t2_smarticle = sizeScale * .02122 / 1.0;
	double rho_smarticleArm = 925;
	double rho_smarticleMid = 739.18;

	double p_gain =1;   //.1//.2         //0.133
	double i_gain = 1;// 0.03;	 //.5//.225						//0.05
	double d_gain = 1; //.0025 //.01       //0.0033


#endif

	// double t_smarticle 	= sizeScale * .00254;
	// double t2_smarticle	= sizeScale * .001;
double vol = (t2_smarticle) * (t_smarticle)* (w_smarticle + 2 * (l_smarticle));
////robot smarticle geometry
//double w_smarticle = 0.046; //4.6cm
//double l_smarticle = 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
//double t_smarticle = .03;
//double t2_smarticle = .021;

double collisionEnvelope = .1 * t2_smarticle;

SystemGeometry *sys;
bool bucket_exist = true;

int read_from_file = 0; //-1 write 0 do nothing 1 read 2 read and write 
bool povray_output = false;
int out_fps = 30;
const std::string out_dir = "PostProcess";
const std::string pov_dir_mbd = out_dir + "/povFilesSmarticles";
int numPerLayer =5;
bool placeInMiddle = false;	/// if I want make a single smarticle on bottom surface

//ChVector<> Cbucket_interior_halfDim = sizeScale * ChVector<>(.05, .05, .025);
//double bucket_rad = sizeScale*0.034;
//double bucket_rad = sizeScale*0.02;
//double bucket_rad = sizeScale*0.022;
//	double bucket_rad = sizeScale*0.04;

std::vector<std::shared_ptr<ChBody>> bucket_bod_vec;

//ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(.1, .1, .05);

double percentToMoveToGlobal = 1.0/800.0;
double percentToChangeStressState = 0;
double max_z = 0;
double rampInc = 1;
double drum_freq = 1; //omega=2*PI*freq
double box_ang =-40 *D2R;
double drum_omega = drum_freq*2*PI;
double pctActive = 1.0;
double inc = 0.00001;
double angle1 = 90;
double angle2 = 90;
double vibAmp = 5 * D2R; //vibrate by some amount of degrees back and forth
int videoFrameInterval = 1/(out_fps*dT); //dt = [sec/step], fps=[frames/sec] --> 1/(dt*fps)=[(sec*steps)/(sec*frames)]=[steps/frame]

int smarticleHopperCount = 0;

// =====================================================================================================
class MyChCustomCollisionPointCallback : public ChSystem::ChCustomCollisionPointCallback
	{
	public:
		/// Callback used to report contact points being added to the container.
		/// This must be implemented by a child class of ChAddContactCallback

		virtual void ContactCallback(
			const collision::ChCollisionInfo& mcontactinfo,  ///< get info about contact (cannot change it)
			ChMaterialCouple& material                       ///< you can modify this!
			)
		{
			//int bucketId = bucket->GetIdentifier();
			//if (mcontactinfo.modelA->GetPhysicsItem()->GetIdentifier() == bucketId || mcontactinfo.modelB->GetPhysicsItem()->GetIdentifier() == bucketId)
			//{
				//GetLog() << material.sliding_friction << " :kfric sfric:" << material.static_friction<<"\n"; 
			//	if (bucket->GetMaterialSurface()->GetKfriction() == 0 || bucket->GetMaterialSurface()->GetSfriction() == 0)
			//	{
			//		material.static_friction = 0;
			//		material.sliding_friction = 0;

			//	}
			//}
 		}

};

class MyBroadPhaseCallback : public collision::ChBroadPhaseCallback {
public:
	//skips collision between smarticle arms and center body 
	virtual bool BroadCallback(collision::ChCollisionModel* mmodelA,  ///< pass 1st model
		collision::ChCollisionModel* mmodelB)   ///< pass 2nd model
	{
		return (!(abs(mmodelA->GetPhysicsItem()->GetIdentifier() - mmodelB->GetPhysicsItem()->GetIdentifier()) < 2));
	}
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

void placeSmarticles(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec, ChIrrApp& application, Smarticle * smarticle0)
{
	//be careful about hlaf size and fullsize dimensions
	//angle1=left angle2=right
	//bucketX = half size

	double r = sys->bucket_rad;
retry:


	double xposi, yposi, zposi, bcx, bcy, x0, y0, lft, rgt, top, bot,bucketX,bucketY =0;
	switch (bucketType)
	{
		case HOPPER:
		{
			//date 7/8/16 in lab notebook #3 for geometry
			//zpos->xpos ypos->ypos xpos->zpos
			double h = sys->bucket_interior_halfDim.z;
			double t = sys->bucket_half_thick; //bucket thickness redefined here for easier to read code
			double sH = (h - t) / sin(box_ang); //side height


			bucketX = sH;
			bucketY = r;
			bcx = sys->bucket->GetPos().z;
			bcy = sys->bucket->GetPos().y;
			x0 = 1.5*sH*genRand();
			y0 = genRand(0, r);
			lft = bcy + (bucketX);
			rgt = bcy - (bucketX);
			top = bcx + (bucketY);
			bot = bcx - (bucketY);
			/*myPos = sys->bucket_ctr + ChVector<>(sin(ang * i + phase) *(bucket_rad / 2 + genRand(w)),
				cos(ang*i + phase)*(bucket_rad / 2 + genRand(w)),
				zpos);*/


			xposi = x0 + bcx - bucketX + sys->bucket_half_thick;
			yposi = y0 + bcy - bucketY + sys->bucket_half_thick;
			zposi = (-yposi - 2 * sys->bucket_half_thick)*tan(Quat_to_Angle(ANGLESET_RXYZ, sys->bucket->GetRot()).x) + t_smarticle / 1.99; //tangent may need to be fixed see buckrotx above
			break;
		}
		case BOX: case FLATHOPPER:
		{
			bucketX = sys->boxdim.x;
			bucketY = sys->boxdim.y;
			bcx = sys->bucket->GetPos().x;
			bcy = sys->bucket->GetPos().y;
			x0 = (bucketX - 4 * sys->bucket_half_thick) * 2 * genRand();
			y0 = (bucketY - 4 * sys->bucket_half_thick) * 2 * genRand();
			lft = bcx + (-sys->bucket_half_thick + bucketX);
			rgt = bcx - (-sys->bucket_half_thick + bucketX);
			top = bcy + (bucketY + 2 * sys->bucket_half_thick)*cos(box_ang);
			bot = bcy - (bucketY + 2 * sys->bucket_half_thick)*cos(box_ang);


			xposi = x0 + bcx - bucketX + sys->bucket_half_thick;
			yposi = y0*cos(box_ang) + bcy - (bucketY + 2 * sys->bucket_half_thick)*cos(box_ang);
			zposi = (-yposi - 2 * sys->bucket_half_thick)*tan(Quat_to_Angle(ANGLESET_RXYZ, sys->bucket->GetRot()).x) + t_smarticle / 1.99; //tangent may need to be fixed see buckrotx above
			break;
		}
		default:
		{
			bucketX = sys->boxdim.x;
			bucketY = sys->boxdim.y;
			bcx = sys->bucket->GetPos().x;
			bcy = sys->bucket->GetPos().y;
			x0 = (bucketX - 4 * sys->bucket_half_thick) * 2 * genRand();
			y0 = (bucketY - 4 * sys->bucket_half_thick) * 2 * genRand();
			lft = bcx + (-sys->bucket_half_thick + bucketX);
			rgt = bcx - (-sys->bucket_half_thick + bucketX);
			top = bcy + (bucketY + 2 * sys->bucket_half_thick)*cos(box_ang);
			bot = bcy - (bucketY + 2 * sys->bucket_half_thick)*cos(box_ang);

			xposi = x0 + bcx - bucketX + sys->bucket_half_thick;
			yposi = y0*cos(box_ang) + bcy - (bucketY + 2 * sys->bucket_half_thick)*cos(box_ang);
			zposi = (-yposi - 2 * sys->bucket_half_thick)*tan(Quat_to_Angle(ANGLESET_RXYZ, sys->bucket->GetRot()).x) + t_smarticle / 1.99; //tangent may need to be fixed see buckrotx above
			break;
		}
	}
	
	int overlap = 0;
	int m1 = 0; //vertex of n1;
	int m2 = 0; //vertex of n2;

	//points matlab [floor((0:11)./4)',mod(0:11,4)',floor((0:11)./4)',mod((0:11)+1,4)']
	//box sides

	auto & cs = smarticle0; //current smarticle 
	if (mySmarticlesVec.size() > 0)// more than 1 smarticle exists in system 
	{
		for (int k = 0; k < mySmarticlesVec.size(); k++)//I think I want to go loop through all but current smarticle 
		{

			cs->TransportSmarticle(ChVector<>(xposi, yposi, zposi));

			//cs->RotateSmarticleBy(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x, /*genRand(0, 2*PI)*/ 0, PI/2)));
			cs->SetEdges();

			application.DrawAll();

			application.GetVideoDriver()->endScene();
			application.GetVideoDriver()->beginScene(true, true,
				video::SColor(255, 140, 161, 192));
			for (int n1 = 0; n1 < 12; n1++)
			{

				ChVector<>csV1 = cs->armVerts[(int)(n1 / 4)][(n1) % 4];			//current smarticle vertex 1
				ChVector<>csV2 = cs->armVerts[(int)(n1 / 4)][(n1 + 1) % 4]; //current smarticle vertex 2
				overlap = overlap + overlaptest(csV1.x, csV1.y, csV2.x, csV2.y, lft, bot, rgt, bot);	//test if smart1 edge overlaps with box BOT edge
				overlap = overlap + overlaptest(csV1.x, csV1.y, csV2.x, csV2.y, lft, top, lft, bot);	//test if smart1 edge overlaps with box LFT edge
				overlap = overlap + overlaptest(csV1.x, csV1.y, csV2.x, csV2.y, rgt, bot, rgt, top);	//test if smart1 edge overlaps with box RGT edge
				overlap = overlap + overlaptest(csV1.x, csV1.y, csV2.x, csV2.y, rgt, top, lft, top);	//test if smart1 edge overlaps with box TOP edge
				if (overlap)
					goto retry;

				if (csV1.x > lft || csV1.x < rgt || csV2.x>lft || csV2.x<rgt) // pos side is on left
					goto retry;
				if (csV1.y > top || csV1.y < bot || csV2.y>top || csV2.y<bot) // pos side is on left
					goto retry;

				//now that we know smarticle is at least inside box, assign other smarticle and update
				auto & os = mySmarticlesVec[k]; // other smarticle
				os->SetEdges();

				for (int n2 = 0; n2 < 12; n2++)
				{
					ChVector<>osV1 = os->armVerts[(int)(n2 / 4)][n2 % 4]; //current smarticle vertex 1
					ChVector<>osV2 = os->armVerts[(int)(n2 / 4)][(n2 + 1) % 4]; //current smarticle vertex 2							
					overlap = overlap + overlaptest(csV1.x, csV1.y, csV2.x, csV2.y, osV1.x, osV1.y, osV2.x, osV2.y);			//tests if smart1's edge(n1) intersects any of smart2's edges (n2)
					if (overlap)
						goto retry;
				}

			}
		}
	}
	else//only 1 smarticle exists
	{
		if (bucketType == HOPPER)
		{
			cs->TransportSmarticle(ChVector<>(zposi, zposi, zposi));
		}
		else
		{
			cs->TransportSmarticle(ChVector<>(xposi, yposi, zposi));
		}
		//cs->RotateSmarticleBy(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(PI / 2 + Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x, genRand(-PI, PI), 0)));
		cs->SetEdges();

		application.DrawAll();

		application.GetVideoDriver()->endScene();
		application.GetVideoDriver()->beginScene(true, true,
			video::SColor(255, 140, 161, 192));
		for (int n1 = 0; n1 < 12; n1++)
		{
			ChVector<>csV1 = cs->armVerts[(int)(n1 / 4)][(n1) % 4];			//current smarticle vertex 1
			ChVector<>csV2 = cs->armVerts[(int)(n1 / 4)][(n1 + 1) % 4]; //current smarticle vertex 2
			overlap = overlap + overlaptest(csV1.x, csV1.y, csV2.x, csV2.y, lft, bot, rgt, bot);	//test if smart1 edge overlaps with box BOT edge
			overlap = overlap + overlaptest(csV1.x, csV1.y, csV2.x, csV2.y, lft, top, lft, bot);	//test if smart1 edge overlaps with box LFT edge
			overlap = overlap + overlaptest(csV1.x, csV1.y, csV2.x, csV2.y, rgt, bot, rgt, top);	//test if smart1 edge overlaps with box RGT edge
			overlap = overlap + overlaptest(csV1.x, csV1.y, csV2.x, csV2.y, rgt, top, lft, top);	//test if smart1 edge overlaps with box TOP edge
			if (overlap)
				goto retry;
		}

	}


}
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
	if (argc > 10){
		const char* text = argv[10];
		//percentToMoveToGlobal = atof(text);
		percentToChangeStressState = atof(text);
	}
	if (argc > 11){
		const char* text = argv[11];
		//percentToMoveToGlobal = atof(text);
		saveFrame = atoi(text);
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
	//const std::string simulationParams = out_dir + "/simulation_specific_parameters.txt";
	//simParams.open(simulationParams.c_str(), std::ios::app);
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
	simParams << "ang1:" << angle1<< std::endl<<
		"ang2:"<<angle2<<	std::endl <<
		"l_smarticle: " << l_smarticle << std::endl <<
		"l_smarticle mult for w (w = mult x l): " << l_smarticle / w_smarticle << std::endl <<
		"read from file: " << read_from_file << std::endl <<
		"dT: " << dT << std::endl << std::endl <<
		"tFinal: " << tFinal << std::endl <<
		"vibrate start: " << vibrateStart << std::endl <<
		"Active Percent: " << pctActive << std::endl <<
		"Start Angles: " << angle1 << " " << angle2<< std::endl;

	simParams << "Smarticle volume: " << vol << std::endl;
	simParams << "Smarticle rhos: arm: " <<rho_smarticleArm<< " mid: "<<rho_smarticleMid <<std::endl;

	//copy smarticle checkpoint if used to PostProcess folder
	if (read_from_file>=1)
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
	mphysicalSystem.SetSolverType(ChSystem::SOLVER_DEM);
	//mphysicalSystem.SetIntegrationType(ChSystem::INT_EULER_IMPLICIT_PROJECTED);
	mphysicalSystem.SetMaxItersSolverSpeed(50+.3*numPerLayer*numLayers);
  mphysicalSystem.SetMaxItersSolverStab(0);   // unuseful for Anitescu, only Tasora uses this
  mphysicalSystem.SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  mphysicalSystem.SetSolverWarmStarting(true);
  mphysicalSystem.SetUseSleeping(false);
  mphysicalSystem.Set_G_acc(ChVector<>(0, 0, gravity));
	vol = (t2_smarticle)* (t_smarticle)* (w_smarticle + 2 * (l_smarticle));
	//mphysicalSystem.SetTolForce(.0005);
	//mphysicalSystem.SetTol(.0001);
	//mphysicalSystem.SetMinBounceSpeed(.3);
	simParams.close();
}
// =============================================================================
int overlaptest(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
	double x = -((-x3 + x4)*(x2*y1 - x1*y2) + (x1 - x2)*(x4*y3 - x3*y4)) / ((x3 - x4)*(y1 - y2) + (x1 - x2)*(-y3 + y4));
	double y = -(x2*y1*y3 - x4*y1*y3 - x1*y2*y3 + x4*y2*y3 - x2*y1*y4 + x3*y1*y4 + x1*y2*y4 - x3*y2*y4) / (-x3*y1 + x4*y1 + x3*y2 - x4*y2 + x1*y3 - x2*y3 - x1*y4 + x2*y4);
	if (((x - x1)*(x - x2) <= 0) && ((x - x3)*(x - x4) <= 0))// % overlap
		return 1;
	else
		return 0;
}
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
		zpos = std::min(3 * sys->bucket_interior_halfDim.z, z) + w_smarticle;

		switch (bucketType){
		case DRUM:
			zpos = SaturateValue(z, sys->bucket_rad);
			myPos = sys->bucket_ctr + ChVector<>(genRand(sys->bucket_interior_halfDim.z / 2.5),
				genRand(sys->bucket_interior_halfDim.z / 2.5),
				zpos);
			break;
		case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
			if (!placeInMiddle)
			{
				myPos = sys->bucket_ctr + ChVector<>(sin(ang * i + phase) *(sys->bucket_rad / 2.2),
					cos(ang*i + phase)*(sys->bucket_rad / 2.2),
					std::max(sys->bucket_interior_halfDim.z * 2, zpos));
				dropSpeed = ChVector<>(0, 0, gravity*timeForDisp / 2.0 - 2 * w_smarticle / timeForDisp);
				myRot = ChQuaternion<>(genRand(-1, 1), genRand(-1, 1), genRand(-1, 1), genRand(-1, 1));
			}
			else////////////place in center of bucket on bucket bottom
			{
				myPos = sys->bucket_ctr + ChVector<>(0, -t_smarticle*1.45, sys->bucket_bott->GetPos().z + t_smarticle);
				dropSpeed = VNULL;
				myRot = Q_from_AngAxis(-PI_2, VECT_X);
			}
			////////////////////////////////////



			break;
		case HOPPER:
			myPos = sys->bucket_ctr + ChVector<>(sin(ang * i + phase) *(sys->bucket_rad / 2 + genRand(w)),
				cos(ang*i + phase)*(sys->bucket_rad / 2 + genRand(w)),
				zpos);

			myRot = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(genRand(-PI, PI), genRand(-PI, PI), genRand(-PI, PI)));
			break;
		case BOX: case FLATHOPPER:
		{
			//myPos = sys->bucket_ctr + ChVector<>((2*MyRand()-1)*.9*sys->bucket_interior_halfDim.x,
			//	(2*MyRand()-1)*.9*sys->bucket_interior_halfDim.y ,
			//	sys->bucket_interior_halfDim.z+w_smarticle);
			//myRot = ChQuaternion<>(2 * MyRand() - 1, 2 * MyRand() - 1, 2 * MyRand() - 1, 2 * MyRand() - 1);
			////dropSpeed = ChVector<>(0, 0, gravity*timeForDisp / 2.0 - 2 * w_smarticle / timeForDisp);
			//dropSpeed = VNULL;

			/////////////////////////////place in specific location///////////////////////////
			//dropSpeed = VNULL;
			//myRot = Q_from_AngAxis(PI_2, VECT_X);
			//ChQuaternion<> buckRot = sys->bucket->GetRot();
			//double buckRotAngx = Quat_to_Angle(ANGLESET_RXYZ, buckRot).x;
			////myRot = buckRot*Angle_to_Quat(ANGLESET_RXYZ,ChVector<>(PI/2,PI,0));
			////myRot = buckRot*Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(genRand(0, 2*PI), genRand(0, 2*PI), genRand(0, 2*PI)));

			////ypos genRand(-PI,PI)

			//myRot = buckRot*Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(PI / 2, genRand(-PI, PI), 0));

			//double xPos = genRand(-3, 3)*t2_smarticle / 1.25;
			//double yPos = (i - 4.2) * 2 * t2_smarticle;
			//myPos = sys->bucket_ctr + ChVector<>(xPos, yPos, (-yPos - 2 * sys->bucket_half_thick)*tan(buckRotAngx) + t_smarticle / 1.99);
			////myPos = sys->bucket_ctr + ChVector<>(xPos, yPos, (-yPos - 2 * sys->bucket_half_thick)*tan(buckRotAngx) + 3*t_smarticle) ;
			/////////////////////////////place in specific location///////////////////////////


			//////////////////////changed////////////////////
			dropSpeed = VNULL;
			myRot = Q_from_AngAxis(PI_2, VECT_X);
			ChQuaternion<> buckRot = sys->bucket->GetRot();
			double buckRotAngx = Quat_to_Angle(ANGLESET_RXYZ, buckRot).x;
			myRot = buckRot*Angle_to_Quat(ANGLESET_RXYZ,ChVector<>(PI/2,PI,0));
			//myRot = buckRot*Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(genRand(0, 2*PI), genRand(0, 2*PI), genRand(0, 2*PI)));

			//
			//double xPos = genRand(-3, 3)*t2_smarticle / 1.25;
			//double yPos = (i - 4.2) * 2 * t2_smarticle;

			double xPos = 0;// genRand(-3, 3)*t2_smarticle / 1.25;
			double yPos = (i)* 1.5 * t2_smarticle;

			myPos = sys->bucket_ctr + ChVector<>(xPos, yPos, (-yPos - 2 * sys->bucket_half_thick)*tan(buckRotAngx) + t_smarticle / 1.99);
			myPos = sys->bucket_ctr + ChVector<>(xPos, yPos, (-yPos - 2 * sys->bucket_half_thick)*tan(buckRotAngx) + 3*t_smarticle) ;
			//////////////////////changed///////////////////
			break;
		}

		default:

			myPos = sys->bucket_ctr + ChVector<>(sin(ang * i + phase) *(sys->bucket_rad / 2 + genRand(-w / 2, w / 2)),
				cos(ang*i + phase)*(sys->bucket_rad / 2 + genRand(-w / 2, w / 2)),
				zpos);
			myRot = ChQuaternion<>(genRand(-1, 1), genRand(-1, 1), genRand(-1, 1), genRand(-1, 1));
			break;
		}

		myRot.Normalize();
		/////////////////flat rot/////////////////
		Smarticle * smarticle0 = new Smarticle(&mphysicalSystem);
		smarticle0->Properties(mySmarticlesVec.size(),smartIdCounter * 4,
			rho_smarticleArm, rho_smarticleMid, mat_smarts,
			collisionEnvelope,
			//l_smarticle+t2_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
			l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
			sOmega,
			true,
			myPos,
			myRot,
			angle1*D2R, angle2*D2R);
		smartIdCounter++;
		if (genRand() < 1)//to reduce amount visualized amount
			smarticle0->visualize = true;
		smarticle0->populateMoveVector();
		smarticle0->SetAngles(angle1, angle2, true);
		smarticle0->SetInitialAngles();
		//if (i == 0)
		//{
		//	smarticle0->active = false;//##################
		//}
		smarticle0->Create();
		smarticle0->setCurrentMoveType((MoveType)Smarticle::global_GUI_value);
		smarticle0->vib.emplace_back(angle1*D2R, angle2*D2R);
		//must put this line because of how linspace is written;
		smarticle0->AssignState(VIB);
		smarticle0->GenerateVib(angle1*D2R, angle2*D2R);
		smarticle0->AssignState(Smarticle::global_GUI_value);
		smarticle0->activateStress = percentToChangeStressState;
		//smarticle0->ss.emplace_back(angle1, angle2);
		//smarticle0->midTorque.emplace_back(angle1*D2R + vibAmp, angle2*D2R + vibAmp);
		//smarticle0->midTorque.emplace_back(angle1*D2R + vibAmp, angle2*D2R + vibAmp);
		//GetLog()<< "\nMASS:"<<smarticle0->GetMass() <<"\n";
		double bucketX = sys->boxdim.x;
		double bucketY = sys->boxdim.y;

		

		//if (bucketType == BOX || bucketType == FLATHOPPER || bucketType == HOPPER)
		//{
		//	dropSmartsInLine()
		//}



		//if (bucketType == BOX || bucketType == FLATHOPPER || bucketType == HOPPER)
		//{
		//	placeSmarticles(mphysicalSystem, mySmarticlesVec, application, smarticle0);
		//}

		mySmarticlesVec.emplace_back((Smarticle*)smarticle0);
		GetLog() << "Smarticles: " << mySmarticlesVec.size() << "\n";
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
//std::shared_ptr<ChBody> create_flathopper(int num_boxes, int id, bool overlap, CH_SYSTEM* mphysicalSystem, std::shared_ptr<ChMaterialSurfaceBase> wallMat, int ridges = 5)
//{//essentially the same as create cyl container except made it bigger and added ridges
//	bucket = utils::CreateBoxContainer(mphysicalSystem, 1000000000, wallMat,
//		boxdim, sys->bucket_half_thick, sys->bucket_ctr, Q_from_AngAxis(-box_ang, VECT_X), true, false, true, false);
//	//bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
//	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));
//	bucket->AddAsset(bucketTexture);
//	bucket->SetCollide(true);
//	bucket->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
//	/*bucket_bott->GetCollisionModel()->SetFamily(1);
//	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);*/
//	//CreateBucket_bott(mphysicalSystem);
//	bucket_bott->SetCollide(false);
//	bucket_bott->SetPos(ChVector<>(5, 5, 5));
//	return bucket;
//}


void CreateMbdPhysicalSystemObjects(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec) {
	/////////////////
	// Ground body
	////////////////

	//sys->mat_wall->SetFriction(wall_fric);
	sys->create_Ground();
	sys->ground->GetCollisionModel()->SetFamily(sys->envFamily);
	sys->ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	sys->create_Container();
	sys->bucket->GetCollisionModel()->SetFamily(sys->envFamily);
	sys->bucket->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	sys->bucket_bott->GetCollisionModel()->SetFamily(sys->envFamily);
	sys->bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	///!@#$%
	//mat_wall->SetFriction(wall_fric); //steel- plexiglass   (plexiglass was outer cylinder material) // .6 for wall to staple using tan (theta) tested on 7/20
	
	
	mat_smarts->SetFriction(smart_fric);
	

	// 1: create bucket

		
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

    ++out_frame;
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
				zMax = bodyPtr->GetPos().z - sys->bucket_bott->GetPos().z;
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
void printFlowRate(double time,int count,bool end = false) //SAVE smarticle gaitType out for reference!
{

	flowRate_of << time << ", " << count << ", " << Smarticle::global_GUI_value<<std::endl;
}
// =============================================================================
void drawGlobalCoordinateFrame(CH_SYSTEM& mphysicalSystem)
{	
	double len = w_smarticle*2;
	double rad = t_smarticle/2;
	ChVector<> pos = sys->bucket_ctr + ChVector<>(2.5*sys->bucket_rad, 0, sys->bucket_interior_halfDim.z);

	auto xaxis = std::make_shared<ChBody>();
	auto yaxis = std::make_shared<ChBody>();
	auto zaxis = std::make_shared<ChBody>();


	//xaxis->SetPos(pos + ChVector<>(len - rad / 2, 0, 0));
	//yaxis->SetPos(pos + ChVector<>(0, len - rad / 2, 0));
	//zaxis->SetPos(pos + ChVector<>(0, 0, len - rad));
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
	xaxis->GetCollisionModel()->SyncPosition();
	yaxis->GetCollisionModel()->SyncPosition();
	zaxis->GetCollisionModel()->SyncPosition();

}



void recycleSmarticles(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mySmarticlesVec)
{
	double pos = -.75*sys->bucket_interior_halfDim.z;//z position below which smarticles are regenerated above pile inside container
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
				sPtr->TransportSmarticle(sys->bucket_ctr + ChVector<>(
					sin(ang*inc + rp)*(sys->bucket_rad / 2 + 4*w_smarticle*(genRand(-0.5,0.5))),
					cos(ang*inc + rp)*(sys->bucket_rad / 2 + w_smarticle*(genRand(-0.5, 0.5))),
					sys->bucket_interior_halfDim.z*2
					));

				//sPtr->SetSpeed(sPtr->GetArm(1)->GetPos_dt() / 4);
				sPtr->SetSpeed(ChVector<>(0, 0, -9.8*.01 / 2.0 - w_smarticle / .01));
			}
			else
			{

				sPtr->TransportSmarticle(sys->bucket_ctr + ChVector<>(
					sin(ang*inc + rp)*(sys->bucket_rad / 2 + w_smarticle*(genRand(-0.5, 0.5))),
					cos(ang*inc + rp)*(sys->bucket_rad / 2 + w_smarticle*(genRand(-0.5, 0.5))),
					sys->bucket_interior_halfDim.z*1.75
					));
			}


			++recycledSmarticles;
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
		if (bodyPtr->GetPos().z < -5.0*sys->bucket_interior_halfDim.z) {
			bodyPtr->SetBodyFixed(true);
		}
	}
}
// =============================================================================

void FixRotation(CH_SYSTEM& mphysicalSystem, Smarticle* sPtr) //reduces rotation speed by half 
{
	//if (sPtr->GetArm(0)->GetRot_dt().GetVector().Length2() > 10000)//added this because small arms can start to spin uncontrollably
	//{
	//	sPtr->GetArm(0)->SetRot_dt(sPtr->GetArm(0)->GetRot() / 2);
	//	sPtr->GetArm(1)->SetRot_dt(sPtr->GetArm(1)->GetRot() / 2);
	//	sPtr->GetArm(2)->SetRot_dt(sPtr->GetArm(2)->GetRot() / 2);
	//	GetLog() << "\n\n******WARNING******\n arm is rotating too fast and FixRotation method is running\n******WARNING******\n";
	//}
	
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
		if (bucketType != HOPPER || bucketType ==FLATHOPPER)
		{

		}

		switch (bucketType)
		{
		case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
			if (!IsInRadial(sPtr->Get_cm(), sys->bucket_bott->GetPos() + ChVector<>(0, 0, sys->bucket_interior_halfDim.z), ChVector<>(sys->bucket_rad*3, sys->bucket_bott->GetPos().z, sys->bucket_bott->GetPos().z + 4 * sys->bucket_interior_halfDim.z)))
			{
				EraseSmarticle(mphysicalSystem, myIter, *sPtr, mySmarticlesVec);
				GetLog() << "\nRemoving smarticle outside system \n";
				continue;
			}
			++myIter;
			break;


		case HOPPER:
			if (!IsInRadial(sPtr->GetArm(1)->GetPos(), sys->bucket->GetPos(), ChVector<>(2 * sys->bucket_rad, -4.0*sys->bucket_interior_halfDim.z, 4.0*sys->bucket_interior_halfDim.z)))
			{
				//double z = Find_Max_Z(mphysicalSystem, mySmarticlesVec);
				double x = sys->bucket_rad / 2 + w_smarticle;
				double y = sys->bucket_rad / 2 + w_smarticle;
				double z = genRand(.2, .6);
				//EraseSmarticle(mphysicalSystem, myIter, *sPtr, mySmarticlesVec);
				ChVector<> myPos = sys->bucket_ctr + ChVector<>(genRand(-x,x),
					genRand(-y, y),
					z);
				sPtr->TransportSmarticle(myPos);
				sPtr->SetSpeed(VNULL);
				GetLog() << "\nRemoving smarticle outside hopper \n";

			}
			else{ ++myIter; }
			break;
		case FLATHOPPER:
			if (sPtr->Get_cm().z < -3.2*sys->bucket_interior_halfDim.z)
			{
				EraseSmarticle(mphysicalSystem, myIter, *sPtr, mySmarticlesVec);
				smarticleHopperCount++;
				printFlowRate(mphysicalSystem.GetChTime(), smarticleHopperCount);
				GetLog() << "\nRemoving Smarticle far below container \n";
				continue;
			}
			else{
				++myIter;
			}
			break;
		default:
			++myIter;
			break;
		}
	}


}
void PrintRingPos(CH_SYSTEM* mphysicalSystem, int tstep, std::shared_ptr<ChBody>ring)
{
		ringPos << tstep << ", " << ring->GetPos().x << ", " << ring->GetPos().y << ", " << ring->GetPos().z << ", " << Smarticle::global_GUI_value << std::endl;
}
void PrintStress(CH_SYSTEM* mphysicalSystem, int tstep, double zmax,double cylrad) //TODO include knobs in calculation
{


	ChVector<> temp = bucket_bod_vec.at(1)->GetPos();
	double currBuckRad = sqrt(temp.x*temp.x + temp.y*temp.y) - sys->bucket_half_thick / 5.0;//sys->bucket_half_thick/5 is how wall thickness is defined!
	//GetLog() << sys->bucket_half_thick<< "thick\n";
	//showForce(mphysicalSystem)/(PI*2*cylrad*zmax)
	double force = showForce(mphysicalSystem);
	//GetLog() << "\nforce:" << force;
	stress_of << mphysicalSystem->GetChTime() << ", " << force <<","<< Smarticle::global_GUI_value <<", "<< currBuckRad<< std::endl;
	//stress_of.close();
}
void PrintStress2(CH_SYSTEM* mphysicalSystem, int tstep, double zmax, double cylrad, std::vector<Smarticle*> mySmarticlesVec) //TODO include knobs in calculation
{

	bool printAllSmarticleInfo = true;
	static int frame = 0;


	//else {
	//	stress_of.open(stress.c_str(), std::ios::app);
	//}
	//GetLog() << sys->bucket_half_thick<< "thick\n";
	//showForce(mphysicalSystem)/(PI*2*cylrad*zmax)
	ChVector<> temp;
	double currBuckRad;
	switch (bucketType)
	{
	case CYLINDER:
		temp = bucket_bod_vec.at(1)->GetPos();
		currBuckRad = sqrt(temp.x*temp.x + temp.y*temp.y) - sys->bucket_half_thick / 5.0;//sys->bucket_half_thick/5 is how wall thickness is defined!
		stress_of << mphysicalSystem->GetChTime() << ", " << 0 << ", " << Smarticle::global_GUI_value << ", " << currBuckRad << ", " << 0 << std::endl; //final 0 is a placeholder 
		break;
	case STRESSSTICK: case KNOBCYLINDER:
		temp = bucket_bod_vec.at(1)->GetPos();
		currBuckRad = sqrt(temp.x*temp.x + temp.y*temp.y) - sys->bucket_half_thick / 5.0;//sys->bucket_half_thick/5 is how wall thickness is defined!
		stress_of << mphysicalSystem->GetChTime() << ", " << showForce(mphysicalSystem) << ", " << Smarticle::global_GUI_value << ", " << currBuckRad << ", " << 0 << std::endl;
		break;
	case BOX:
		stress_of << mphysicalSystem->GetChTime() << ", " << angle1 << ", " << Smarticle::global_GUI_value << ", " << box_ang << ", " << angle2 << std::endl;
		break;

	}


	if (printAllSmarticleInfo)
	{
		for (size_t i = 0; i < mySmarticlesVec.size(); i++)
		{
			//works for plotlazy matlabfile
			stress_of << mySmarticlesVec[i]->GetAngle(0,true) << ", " << mySmarticlesVec[i]->GetAngle(1,true) << ", " << mySmarticlesVec[i]->moveType << ", " << mySmarticlesVec[i]->Get_cm().z << std::endl;
			//stress_of << mySmarticlesVec[i]->active << ", " << mySmarticlesVec[i]->Get_cm().x << ", " << mySmarticlesVec[i]->Get_cm().y << ", " << mySmarticlesVec[i]->Get_cm().z << std::endl;
		}
		stress_of << "#EF" <<frame<< std::endl;
	}
	//stress_of.close();
	frame = frame + 1;
}
void PrintFractions(CH_SYSTEM& mphysicalSystem, int tStep, std::vector<Smarticle*> mySmarticlesVec) {

	static int stepSave = 10;
	if (tStep % stepSave != 0) return;
	double zComz = 0;
	double totalTorque = 0;
	//static std::shared_ptr<ChBody> grid;  //uncomment to visualize vol frac boxes
	//static bool a = false;						//uncomment to visualize vol frac boxes


	//double sqSize = w_smarticle; // try increasing!
	//int rowSize = std::ceil(sys->bucket_rad*2/sqSize);
	double sqSizex = (w_smarticle+2*t2_smarticle); // try increasing!
	double sqSizey = (l_smarticle+2*t2_smarticle); // try increasing!
	int colSize = std::ceil(sys->bucket_rad*2/sqSizex);
	int rowSize = std::ceil(sys->bucket_rad*2/sqSizey);
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
	static int countInside =0;

	//double zMax = Find_Max_Z(mphysicalSystem,mySmarticlesVec);
	ChVector<> bucketMin = sys->bucket_bott->GetPos();


	// *** remember, 2 * sys->bucket_half_thick is needed since bucket is initialized inclusive. the half dims are extended 2*sys->bucket_half_thick from each side
	ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, sys->bucket_interior_halfDim.z);
	double totalVolume2 = 0;
	int countInside2 = 0;
	double volumeFraction = 0;
	

	switch (bucketType)
	{
	case BOX:
		zMax = Find_Max_Z(mphysicalSystem, mySmarticlesVec);
		zMax = std::min(zMax, bucketMin.z + 2 * sys->bucket_interior_halfDim.z);
		for (size_t i = 0; i < mySmarticlesVec.size(); i++) {
			Smarticle* sPtr = mySmarticlesVec[i];
			if (IsIn(sPtr->Get_cm(), bucketCtr - sys->bucket_interior_halfDim, bucketCtr + sys->bucket_interior_halfDim + ChVector<>(0, 0, 2.0 * sys->bucket_half_thick))) {
				countInside2++;
				totalVolume2 += sPtr->GetVolume();
			}
		}

		volumeFraction = totalVolume2 / (4.0 * sys->bucket_interior_halfDim.x * sys->bucket_interior_halfDim.y * (zMax - bucketMin.z));
		break;
	case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
		for (size_t i = 0; i < mySmarticlesVec.size(); i++) {
			Smarticle* sPtr = mySmarticlesVec[i];
			//isinradial rad parameter is Vector(bucketrad,zmin,zmax)
			if (IsInRadial(sPtr->Get_cm(), bucketCtr, ChVector<>(sys->bucket_rad, bucketMin.z, bucketMin.z + 2.0*sys->bucket_interior_halfDim.z))) {
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
				totalTorque += sPtr->GetReactTorqueVector(0).z +sPtr->GetReactTorqueVector(1).z;
				//zMax = std::max(zMax, sPtr->GetArm(1)->GetPos().z- bucketMin.z);

			}
		}
		volumeFraction = countInside2*vol / (max2*PI*sys->bucket_rad*sys->bucket_rad);
		//GetLog() << vol << " " << countInside2 << " " << sys->bucket_rad << " " << zMax << " " << volumeFraction << "\n";
		//GetLog() << "phi=" << volumeFraction << "\n";
		zComz = zComz / countInside2;
		totalTorque = totalTorque / (countInside2 * 2.0); //multiply by 2 (2 arms for each smarticle)
		break;
	case FLATHOPPER:
	{
		zComz = 0;
		zMax = Find_Max_Z(mphysicalSystem, mySmarticlesVec);
		zMax = std::min(zMax, bucketMin.z + 2 * sys->bucket_interior_halfDim.z);
		max2 = 0;
		for (size_t i = 0; i < mySmarticlesVec.size(); i++) 
		{
			Smarticle* sPtr = mySmarticlesVec[i];
			
			com = sPtr->Get_cm() - ChVector<>(0, 0, sys->bucket_bott->GetPos().z);
			max2 = std::max(max2, com.z);


			countInside2 = smarticleHopperCount;
			totalTorque += mySmarticlesVec[i]->GetTotalTorque();
			zComz += sPtr->GetArmTorque(1);
		}
		//volumeFraction = (countInside2*vol) / (max2*boxdim.x*abs(cos(box_ang)));
		volumeFraction = 0;
	}
	default:
		break;
	}
	//totalEnergy used to be meanOT
	vol_frac_of << mphysicalSystem.GetChTime() << ", " << countInside2 << ", " << volumeFraction << ", " << zMax << ", " << zComz << ", " << totalTorque << ", " << Smarticle::global_GUI_value << std::endl;
	//vol_frac_of.close();
	return;
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
	double pctglob[30] =
	{ 0.0009,
	0.0013,
	0.0016,
	0.0023,
	0.0027,
	0.0030,
	0.0041,
	0.0047,
	0.0052,
	0.0074,
	0.0085,
	0.0097,
	0.0151,
	0.0174,
	0.0197,
	0.0233,
	0.0247,
	0.0262,
	0.0343,
	0.0371,
	0.0399,
	0.0449,
	0.0473,
	0.0496,
	0.0620,
	0.0698,
	0.0776,
	0.0975,
	0.1199,
	0.1423 };

	double t = mphysicalSystem.GetChTime(); 
	
	for (size_t i = 0; i < mySmarticlesVec.size(); i++) {
		double tor1 = 0;
		double tor2 = 0; 
		if (mySmarticlesVec[i]->active)
		{
			mySmarticlesVec[i]->updateTorqueDeque();
			tor1 = std::get<0>(mySmarticlesVec[i]->torqueAvg);
			tor2 = std::get<1>(mySmarticlesVec[i]->torqueAvg);
		}
		int moveType=0;
		///////////////////random chance at current timestep for smarticle to not move to globalValue, models real life delay for smarticles to start motion to current state
		if (genRand() < 1)
		{
			//if (genRand() < pctglob[i] / 10)
			//{
			moveType = Smarticle::global_GUI_value;
			//}
		}
		else
		{
			moveType = mySmarticlesVec[i]->prevMoveType;
		}
		mySmarticlesVec[i]->ControllerMove(moveType, tor1, tor2);
		mySmarticlesVec[i]->steps = mySmarticlesVec[i]->steps + 1;

	}
}


void create_spring_cir(CH_SYSTEM& mphysicalSystem)
{
	int n_parts = 50;
	double rad = .0055;
	double k = 1000.0;
	double r = 1;
	double m = 0.005;
	double ropeRad = w_smarticle;
	ChVector<>p1 = ChVector<>(0, 0, .01);
	double ang = 2.0 * PPI / n_parts;
	ChVector<> pPos = (0, 0, 0);  //position of each plate
	ChQuaternion<> quat = QUNIT; //rotation of each plate
	std::vector<std::shared_ptr<ChBody>> partVec;
	std::shared_ptr<ChBody> partLast;
	for (int i = 0; i < n_parts; i++)
	{
		auto part = std::make_shared<ChBody>();

		pPos = p1 + ChVector<>(sin(ang * i) * (ropeRad),
			cos(ang*i)*(ropeRad),
			0);

		quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));

		part->AddAsset(sys->sphereTexture);

		part->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		part->SetPos(pPos);
		part->GetCollisionModel()->ClearModel();
		utils::AddSphereGeometry(part.get(), rad, pPos, QUNIT, true);
		part->SetMass(m);

		///////
		part->GetCollisionModel()->BuildModel();
		part->GetCollisionModel()->SetFamily(3); 
		part->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
		mphysicalSystem.Add(part);
		////////////
		
		partVec.emplace_back(part);
		part->GetCollisionModel()->SyncPosition();
		if (i > 0)
		{
			auto spring = std::make_shared<ChLinkSpring>();
			spring->Initialize(partVec.at(i), partVec.at(i - 1), false, partVec.at(i)->GetPos(), partVec.at(i - 1)->GetPos(), true, 0);
			//spring->Set_SpringRestLength(rad * 2.1);
			spring->Set_SpringK(k);
			spring->Set_SpringR(r);
			mphysicalSystem.Add(spring);
		}
			if (i == n_parts-1)
		{
			auto springLast = std::make_shared<ChLinkSpring>();
			springLast->Initialize(partVec.at(i), partVec.at(0), false, partVec.at(i)->GetPos(), partVec.at(0)->GetPos(), true, 0);
			//spring->Set_SpringRestLength(rad * 2.1);
			springLast->Set_SpringK(k);
			springLast->Set_SpringR(r);
			mphysicalSystem.Add(springLast);
			partLast = part;
		}

		//part->SetBodyFixed(true);
		part->SetCollide(true);
		

	}
	//partLast->Accumulate_force(ChVector<>(0, .4, 0), VNULL, true);
}
void create_spring(CH_SYSTEM& mphysicalSystem)
{
	double rad = .025;
	double k = 10000.0;
	double r = 4;
	double m = 0.001;
	ChVector<>p1 = ChVector<>(0, 0, .01);
	ChVector<>b1Pos = p1 + ChVector<>(rad, 0, 0);
	ChVector<>b2Pos = p1 + ChVector<>(-rad, 0, 0);
	auto b1 = std::make_shared<ChBody>();
	b1->SetMass(m);
	b1->AddAsset(sys->sphereTexture);
	//b1->SetPos(b1Pos);
	b1->GetCollisionModel()->ClearModel();
	b1->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddSphereGeometry(b1.get(), rad, b1Pos,QUNIT, true);
	b1->GetCollisionModel()->BuildModel();
	b1->SetCollide(true);


	auto b2 = std::make_shared<ChBody>();
	b2->SetMass(m);

	b2->GetCollisionModel()->ClearModel();
	b2->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	//b1->SetPos(b2Pos);
	utils::AddSphereGeometry(b2.get(), rad, b2Pos, QUNIT, true);
	b2->GetCollisionModel()->BuildModel();
	b2->SetCollide(true);
	b2->SyncCollisionModels();
	mphysicalSystem.Add(b1);
	mphysicalSystem.Add(b2);

	
	auto spring = std::make_shared<ChLinkSpring>();

	
	spring->Initialize(b1, b2, false, b1->GetPos(), b2->GetPos(), true, 0);
	//spring->Set_SpringRestLength(rad * 2.1);

	spring->Set_SpringK(k);
	spring->Set_SpringR(r);
	//mphysicalSystem.Add(spring);
	spring->AddAsset(sys->groundTexture);
	b1->Accumulate_force(ChVector<>(0, .5, 0), VNULL,true);
	mphysicalSystem.Add(spring);
}
// =============================================================================
bool SetGait(double time)
{
	//double tm = 3;
	//if (time <= tm*1)
	//	Smarticle::global_GUI_value = 1;
	//else if (time > tm*1 && time <= tm*2)
	//	Smarticle::global_GUI_value = 2;
	//else if (time > tm*2 && time <= tm*3)
	//	Smarticle::global_GUI_value = 1;
	//else if (time > tm*3 && time <= tm*4)
	//	Smarticle::global_GUI_value = 2;
	//else
	//	return true;


	//if (time <= .2)
	//	Smarticle::global_GUI_value = 2;
	//else if (time>.5)
	//	Smarticle::global_GUI_value = 0;
	if (time > 30)
		return true;
	/*else
		Smarticle::global_GUI_value = 1;
*/

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
	Smarticle::global_GUI_value = 0;

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
	//struct tm* timeinfoDiff;
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
	
	CH_SYSTEM mphysicalSystem;
	const std::string simulationParams = out_dir + "/simulation_specific_parameters.txt";
	simParams.open(simulationParams.c_str());

	sys = &SystemGeometry(&mphysicalSystem, bucketType, collisionEnvelope,
		l_smarticle, w_smarticle, t_smarticle, t2_smarticle);
	InitializeMbdPhysicalSystem_NonParallel(mphysicalSystem, argc, argv);
	
	
	
	

	

	simParams << "Job was submitted at date/time: " << asctime(timeinfo) << std::endl;
	//simParams.close();

	// define material property for everything
	//!@#$%
	//mat_wall = std::make_shared<ChMaterialSurface>();
	mat_smarts= std::make_shared<ChMaterialSurface>();
	// Create a ChronoENGINE physical system


	

	videoFrameInterval = 1 / (out_fps*dT); //dt = [sec/step], fps=[frames/sec] --> 1/(dt*fps)=[(sec*steps)/(sec*frames)]=[steps/frame]
	GetLog() << "\npctActive" << pctActive << "\n";
	Smarticle::pctActive = pctActive;
	MyBroadPhaseCallback mySmarticleBroadphaseCollisionCallback;
	mphysicalSystem.GetCollisionSystem()->SetBroadPhaseCallback(&mySmarticleBroadphaseCollisionCallback);

	//MyChCustomCollisionPointCallback myContactCallbackFriction;
	//mphysicalSystem.GetContactContainer()->SetAddContactCallback(myContactCallbackFriction.get());
	//mphysicalSystem.SetCustomCollisionPointCallback(&myContactCallbackFriction);
	//mphysicalSystem.GetContactContainer()->AddCollisionModelsToSystem();


	std::vector<Smarticle*> mySmarticlesVec;
	CreateMbdPhysicalSystemObjects(mphysicalSystem, mySmarticlesVec);

	//simParams.open(simulationParams.c_str(), std::ios::app);
	simParams << "SystemSize: " << sys->boxdim.x << sys->boxdim.y << sys->boxdim.z << std::endl;
	//simParams.close();

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
		core::vector3df(0, 0, 4*sys->bucket_rad)*sizeScale,
		core::vector3df(0, 0, 1*sys->bucket_rad)*sizeScale);

	ChIrrWizard::add_typical_Lights(application.GetDevice(),
		core::vector3df(.0139, -.39, -.0281)*sizeScale,
		core::vector3df(0.0139, .298, -.195)*sizeScale);

	ChIrrWizard::add_typical_Lights(application.GetDevice());
	

	RTSCamera* camera = new RTSCamera(application.GetDevice(), application.GetDevice()->getSceneManager()->getRootSceneNode(),
		application.GetDevice()->getSceneManager(), -1, -50.0f, 0.5f, 0.0005f);
	camera->setTranslateSpeed(0.005);
////camera->setPosition(core::vector3df(0, -.17, .07));
//	//camera->setPosition(core::vector3df(0, bucket_interior_halfDim.y, 4 * sys->bucket_rad / 2.0));
//	//camera->setTarget(core::vector3df(0, 2*bucket_interior_halfDim.y, -sys->bucket_rad / 2.1)); //	camera->setTarget(core::vector3df(0, 0, .01));

	//camera->setPosition(core::vector3df(.0139, -.445, .7));
	//camera->setTarget(core::vector3df(0.0139, .243, -.195)); //	camera->setTarget(core::vector3df(0, 0, .01));

	switch (bucketType)
	{
	case HOPPER:
	{
		if (stapleSize)
		{
			camera->setPosition(core::vector3df(0.00024, -0.181, .102));
			camera->setTarget(core::vector3df(0.00024, -.033, .0759)); //	camera->setTarget(core::vector3df(0, 0, .01));
		}
		else
		{
			camera->setPosition(core::vector3df(0.0139, -.45, .25));
			camera->setTarget(core::vector3df(0.0139, -.3, .24)); //	camera->setTarget(core::vector3df(0, 0, .01));
		}
	}
		break;
	case BOX:
		if (stapleSize)
		{
			camera->setPosition(core::vector3df(0.0139, -0.255, .045));
			camera->setTarget(core::vector3df(0.0139, -.1050, .03)); //	camera->setTarget(core::vector3df(0, 0, .01));
			if (box_ang == 0)
			{
				camera->setPosition(core::vector3df(0.0139, -0.255, .045));
				camera->setTarget(core::vector3df(0.0139, -.1050, .03)); //	camera->setTarget(core::vector3df(0, 0, .01));
			}
		}
		else
		{
			camera->setPosition(core::vector3df(0.0139, -0.65, -.180));
			camera->setTarget(core::vector3df(0.0139, -.50, -.195)); //	camera->setTarget(core::vector3df(0, 0, .01));
			if (box_ang == 0)
			{
				camera->setPosition(core::vector3df(-0.04, -0.0142, .943));
				camera->setTarget(core::vector3df(-0.0144, .019, -.497)); //	camera->setTarget(core::vector3df(0, 0, .01));
			}
		}
		break;
	case FLATHOPPER:
		camera->setPosition(core::vector3df(0.014, -1.34, -.4338));
		camera->setTarget(core::vector3df(0.0139, -.50, -.495)); //	camera->setTarget(core::vector3df(0, 0, .01));

		if (box_ang == 0)
		{
			camera->setPosition(core::vector3df(-0.04, -0.0142, .943));
			camera->setTarget(core::vector3df(-0.0144, .019, -.497)); //	camera->setTarget(core::vector3df(0, 0, .01));
		}

		break;

	case DRUM:

		if (stapleSize)
		{
			camera->setPosition(core::vector3df(-0.0011, -0.115, 0.015));
			camera->setTarget(core::vector3df(-0.0011, 0.035, 1e-8)); //	camera->setTarget(core::vector3df(0, 0, .01));
		}
		else
		{
			camera->setPosition(core::vector3df(0.0139, -0.65, -.180));
			camera->setTarget(core::vector3df(0.0139, -.50, -.195)); //	camera->setTarget(core::vector3df(0, 0, .01));
		}
		
		break;

	case KNOBCYLINDER: case CYLINDER: case STRESSSTICK: case HOOKRAISE:

		if (stapleSize)
		{
			camera->setPosition(core::vector3df(-0.0061, -0.095, 0.03));
			camera->setTarget(core::vector3df(-0.0061, 0.055, 0.015)); //	camera->setTarget(core::vector3df(0, 0, .01));
		}
		else
		{
			camera->setPosition(core::vector3df(0.0139, -0.65, -.180));
			camera->setTarget(core::vector3df(0.0139, -.50, -.195)); //	camera->setTarget(core::vector3df(0, 0, .01));
		}

		break;

	default:
		camera->setPosition(core::vector3df(0.0139, -0.65, -.180));
		camera->setTarget(core::vector3df(0.0139, -.50, -.195)); //	camera->setTarget(core::vector3df(0, 0, .01));
		break;
	}



	ChIrrWizard::add_typical_Lights(application.GetDevice(),
		camera->getPosition(),
		camera->getTarget());

	camera->setNearValue(0.0005f);
	camera->setMinZoom(0.01f);
	camera->setZoomSpeed(0.1f);
	camera->setUpVector(core::vector3df(0, -1, 0));
	auto x = camera->getAbsoluteTransformation().getRotationDegrees();
	x.X = 90;
	auto c= camera->getAbsoluteTransformation();
	c.setRotationDegrees(x);
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
	IrrGui receiver(&application, &mySmarticlesVec,sys);
	// scan all contact
	// note how to add the custom event receiver to the default interface:
	application.SetUserEventReceiver(&receiver);



	receiver.drawAngle();//initialize draw angle
	application.SetShowInfos(true);

#endif


	int stepEnd = int(tFinal / dT);  // 1.0e6;//2.4e6;//600000;//2.4e6 * (.02 * paramsH.sizeScale) /


	switch (bucketType)
	{
	case STRESSSTICK: case HOOKRAISE:// case KNOBCYLINDER:
	{
		double rodLen = sys->bucket_interior_halfDim.z*1.5;
		
		sys->create_CentralColumn(rodLen);
		sys->create_Truss();
		sys->create_Prismatic(sys->stick);

		break;
	}

	case DRUM:
	{
		sys->setUpBucketActuator(sys->bucket->GetRot());
		sys->bucket->SetBodyFixed(true);
		break;
	}
	case BOX:
	{

		sys->setUpBucketActuator(Q_from_AngAxis(PI_2, VECT_Y));
		sys->bucket->SetBodyFixed(true);
		break;
	}

	
	case HOPPER:
	{
		sys->create_Truss();
		sys->create_VibrateLink(omega_bucket, vibration_amp, vibrateStart, sys->bucket);
		break;

	}
	case CYLINDER:
	{
		sys->create_Truss();
		sys->create_VibrateLink(omega_bucket, vibration_amp, vibrateStart, sys->bucket_bott);

		break;
	}
	case KNOBCYLINDER:
	{	
		unsigned int kpr = 4;//knobs per row
		unsigned int rows = 15; //knob per z
		double rodLen = sys->bucket_interior_halfDim.z*2.0;
		sys->create_CentralColumn(rodLen);
		sys->create_Knobs(kpr, rows, rodLen);
		sys->setUpBucketActuator();
		break;

	}
		default:
			break;
	}
	double timeForVerticalDisplacement = 0.015;
	if (bucketType == DRUM)
		timeForVerticalDisplacement = 0.095; // 1.5 for safety proximity .015
	if (bucketType == BOX)
		timeForVerticalDisplacement = .15;
	int numGeneratedLayers = 0;


	if (read_from_file>=1)
	{
		CheckPointSmarticlesDynamic_Read(mphysicalSystem, mySmarticlesVec,application);
		application.AssetBindAll();
		application.AssetUpdateAll();
		numGeneratedLayers = numLayers;
	}

	//create_spring(mphysicalSystem);
	//create_spring_cir(mphysicalSystem);
	double xPos = 0;
	double yPos = 1.5 * t2_smarticle;
	ChVector<> pos2 = sys->bucket_ctr + ChVector<>(xPos, yPos, t2_smarticle);

	std::shared_ptr<ChBody> ring = sys->create_EmptyCylinder(25, true, false, t2_smarticle, sys->bucket_half_thick, ringRad, pos2, false, sys->groundTexture);
	ring->SetIdentifier(455465);
	ring->SetCollide(true);
	ring->SetMass(.001);
	mphysicalSystem.AddBody(ring);
	
	
	PrintRingPos(&mphysicalSystem, mphysicalSystem.GetChTime(),ring);
	
	const std::string stress = out_dir + "/Stress.txt";
	const std::string flowRate = out_dir + "/flowrate.txt";
	const std::string vol_frac = out_dir + "/volumeFraction.txt";


	ringPos.open(out_dir + "/RingPos.txt");
	stress_of.open(stress.c_str());
	flowRate_of.open(flowRate.c_str());
	vol_frac_of.open(vol_frac.c_str());

	ringPos << "# ring rad = " << ringRad << "tstep, x, y, z, globalGUI" << std::endl;;
	stress_of << dT << ", " << out_fps << ", " << videoFrameInterval << ", " << sys->bucket_rad << ", " << bucketType << std::endl;
	flowRate_of << 0 << ", " << 0 << ", " << Smarticle::global_GUI_value << std::endl;

	


//  for (int tStep = 0; tStep < 1; tStep++) {
	//START OF LOOP 
	application.DrawAll();
	
	for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
		double t = mphysicalSystem.GetChTime();
		if (read_from_file < 1)
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
		if (bucketType == DRUM)
		{
			sys->bucket_bott->SetBodyFixed(true);
			sys->rotate_body_sp(t, sys->bucket, sys->bucket_actuator, drum_omega);
		}
		if (bucketType == BOX)
		{
			sys->bucket_bott->SetBodyFixed(true);
			sys->rotate_body_rot(t, sys->bucket, sys->bucket_actuator, Quat_to_Angle(ANGLESET_RXYZ, sys->bucket->GetRot()).x);
		}
		//vibration movement
		if (t > vibrateStart && t < vibrateStart + 3)
		{
			switch (bucketType)
			{
			case HOOKRAISE: case STRESSSTICK:
			{
				sys->stick->SetBodyFixed(false);
				sys->pris_link->SetDisabled(false);
				break;

				//if (sys->pris_engine->IsDisabled())
				//{
				//	sys->stick->SetBodyFixed(false);
				//	sys->pris_engine->SetDisabled(false);
				//
				//}
				//sys->pris_engine->GetDist_dt();
				//break;
			}
			case KNOBCYLINDER:
			{
				double rotSpeed = 2; //rads/sec
				sys->bucket_actuator->SetDisabled(false);
				sys->stick->SetBodyFixed(false);				
				sys->rotate_body_sp(t, sys->stick, sys->bucket_actuator, PPI);
				break;
			}

			case CYLINDER:
			{
				sys->bucket_bott->SetBodyFixed(false);
				sys->vibrate_link->SetDisabled(false);
				break;
			}
			case HOPPER:
			{
				sys->bucket->SetBodyFixed(false);
				sys->vibrate_link->SetDisabled(false);
				break;
			}
			case FLATHOPPER:
			{
				sys->bucket_bott->SetPos(ChVector<>(1, 0, 0));
				bucket_exist = false;
				break;
			}
			default:
				break;
			}
		}


		if ((fmod(t, timeForVerticalDisplacement) < dT) && (mySmarticlesVec.size() < numPerLayer*numLayers) && (numGeneratedLayers == numLayers))
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
		//mphysicalSystem.DoStepDynamics(dT);
		UpdateSmarticles(mphysicalSystem, mySmarticlesVec);



		receiver.drawSmarticleAmt(numGeneratedLayers);
		receiver.drawSuccessful();
		/////////////////////////////////////////////////////////////
		//	double jointClearance = .0065;
		//	double armt = t_smarticle;
		//	double armt2 = .00806 / 2 * sizeScale; //8.06 mm with solar 3.2 without
		//	double l_mod = l_smarticle + 2 * t2_smarticle/2 - jointClearance;
		//	//application.GetVideoDriver()->draw3DLine(vector3df(mySmarticlesVec[0]->GetArm(1)->GetPos().x, mySmarticlesVec[0]->GetArm(1)->GetPos().y, mySmarticlesVec[0]->GetArm(1)->GetPos().z),
		//	//	vector3df(mySmarticlesVec[0]->GetArm(1)->GetPos().x, mySmarticlesVec[0]->GetArm(1)->GetPos().y, mySmarticlesVec[0]->GetArm(1)->GetPos().z + 1), irr::video::SColor(70, 30, 200, 200));
		//	ChVector<>pos(mySmarticlesVec[0]->GetArm(1)->TransformPointLocalToParent(ChVector<>(-w_smarticle/2, 0,0)));
		//	ChVector<>pos2(mySmarticlesVec[0]->GetArm(1)->TransformPointLocalToParent(ChVector<>(-w_smarticle / 2, 0, t2_smarticle / 2)));
		//	ChVector<>pos3(mySmarticlesVec[0]->GetArm(1)->TransformPointLocalToParent(ChVector<>(w_smarticle / 2, 0, t2_smarticle/2)));

		//	vector3df armpos(pos.x+pos2.x, pos.y+pos2.y, pos.z+pos2.z);
		//	vector3df armpos2(-pos2.y, pos2.x, pos2.z);
		//	vector3df armpos3(-pos3.y, pos3.x, pos3.z);
		//	//l + 2 * r2 - jointClearance
		//	application.GetVideoDriver()->setTransform(irr::video::ETS_WORLD, core::IdentityMatrix);
		//	application.GetVideoDriver()->draw3DLine(armpos2,
		//		armpos3, irr::video::SColor(70, 255, 0, 0));

		//	application.GetVideoDriver()->draw3DLine(armpos2 - vector3df(0, 0, 10),
		//		armpos2 + vector3df(0, 0, 10), irr::video::SColor(70, 0, 255, 0));
		//	
		///*	application.GetVideoDriver()->draw3DLine(armpos3 - vector3df(0, 0, 10),
		//		armpos3 + vector3df(0, 0, 10), irr::video::SColor(70, 0, 255, 0));*/

		//	//////////////////////////////////////////////////////////////////////

		application.DrawAll();

		application.GetVideoDriver()->endScene();

#else

		mphysicalSystem.DoStepDynamics(dT);
		UpdateSmarticles(mphysicalSystem, mySmarticlesVec);
#endif
#endif

		if (SetGait(t) == true)
			break;


		if (bucketType == STRESSSTICK || bucketType == KNOBCYLINDER || bucketType == CYLINDER || bucketType == BOX)
		{
			double zmax = Find_Max_Z(mphysicalSystem, mySmarticlesVec);
			PrintStress2(&mphysicalSystem, tStep, zmax, sys->rad, mySmarticlesVec);
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
		
		if (read_from_file == -1|| read_from_file == 2)
		{
			CheckPointSmarticlesDynamic_Write(mySmarticlesVec,
			 		tStep,
			 		mat_smarts,
			 		l_smarticle,
			 		w_smarticle,
			 		t_smarticle,
			 		t2_smarticle,
			 		collisionEnvelope,
					rho_smarticleArm,
					rho_smarticleMid);
		}
		PrintRingPos(&mphysicalSystem, mphysicalSystem.GetChTime(), ring);
  }
	//simParams.open(simulationParams.c_str(), std::ios::app);
	simParams << "Smarticle OT: " << mySmarticlesVec.at(0)->OTThresh << std::endl;

	//need to print out last step at 15 secs
	//if (bucketType == FLATHOPPER)
	//{
	//	recycleSmarticles(mphysicalSystem, mySmarticlesVec);
	//}
  for (int i = 0; i < mySmarticlesVec.size(); i++) {
	  delete mySmarticlesVec[i];
  }
	if (saveFrame)
	{
		receiver.SaveToMovie();
		receiver.DeleteJpegs();
	}
  mySmarticlesVec.clear();

	simParams << "completed"<<std::endl;

	flowRate_of.close();
	simParams.close();
	stress_of.close();
	vol_frac_of.close();
	ringPos.close();
  simParams.close();
  return 0;
}
