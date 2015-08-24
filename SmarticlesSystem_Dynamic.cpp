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
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "utils/ChUtilsCreators.h"     //Arman: why is this
#include "utils/ChUtilsInputOutput.h"  //Arman: Why is this
#include "utils/ChUtilsGenerators.h"

#include <ctime>
#include <stdlib.h>  // system, rand, srand, RAND_MAX
#include "core/ChFileutils.h" // for MakeDirectory
#include "Smarticle.h"
#include "SmarticleU.h"
#include "CheckPointSmarticles.h"

#include <memory>



#if irrlichtVisualization

#ifdef CHRONO_OPENGL
#undef CHRONO_OPENGL
#endif

//#include "unit_IRRLICHT/ChIrrApp.h"
#include "unit_IRRLICHT/ChBodySceneNode.h"
#include "unit_IRRLICHT/ChBodySceneNodeTools.h"
//#include "unit_IRRLICHT/ChIrrTools.h"
#include "unit_IRRLICHT/ChIrrWizard.h"
#include "core/ChRealtimeStep.h"
//#include <irrlicht.h>
#include "assets/ChTexture.h"

using namespace irr;
using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;
#endif


#ifdef CHRONO_OPENGL

#include "chrono_opengl/ChOpenGLWindow.h"
#endif

	using namespace chrono;

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

enum SmarticleType {SMART_ARMS , SMART_U};
enum BucketType { CYLINDER, BOX };
SmarticleType smarticleType = SMART_ARMS;//SMART_U;
BucketType bucketType = CYLINDER;
// =============================================================================

class MyBroadPhaseCallback : public collision::ChBroadPhaseCallback {
  public:
    /// Callback used to report 'near enough' pairs of models.
    /// This must be implemented by a child class of ChBroadPhaseCallback.
    /// Return false to skip narrow-phase contact generation for this pair of bodies.
   virtual bool BroadCallback( collision::ChCollisionModel* mmodelA,  ///< pass 1st model
    							collision::ChCollisionModel* mmodelB   ///< pass 2nd model
                               ) {
    	return (!(abs(mmodelA->GetPhysicsItem()->GetIdentifier() - mmodelB->GetPhysicsItem()->GetIdentifier()) < 3));
    }
};

// =============================================================================


double Find_Max_Z(CH_SYSTEM& mphysicalSystem);
std::ofstream simParams;
ChSharedPtr<ChBody> bucket;
ChSharedPtr<ChBody> bucket_bott;

	int appWidth = 1280;
	int appHeight = 720;
	double sizeScale = 1;
	double gravity = -9.81 * sizeScale;
	
	double vibration_freq = 30;
	
	double omega_bucket = 2 * CH_C_PI * vibration_freq;  // 30 Hz vibration similar to Gravish 2012, PRL
	//double vibration_amp = sizeScale * 0.00055;
	double mGamma = 2.0 * gravity;
	double vibration_amp = mGamma / (omega_bucket*omega_bucket);
	


	//double dT = std::min(0.001, 1.0 / vibration_freq / 200);;//std::min(0.0005, 1.0 / vibration_freq / 200);
	double dT = 0.001;//std::min(0.0005, 1.0 / vibration_freq / 200);
	double contact_recovery_speed = 0.5 * sizeScale;
	double tFinal = 1000;
	double vibrateStart= tFinal-5.0;

	double rho_smarticle = 7850.0 / (sizeScale * sizeScale * sizeScale);
	double rho_cylinder = 1180.0 / (sizeScale * sizeScale * sizeScale);
	ChSharedPtr<ChMaterialSurface> mat_g;
	int numLayers = 100;
	double armAngle = 90;
	double sOmega = 5;  // smarticle omega
	

		// staple smarticle geometry
		double w_smarticle 	= sizeScale * 0.0117;
		double l_smarticle 	= 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	//	double t_smarticle 	= sizeScale * .00127;
	//	double t2_smarticle	= sizeScale * .0005;
		double t_smarticle 	= sizeScale * .00254;
		double t2_smarticle	= sizeScale * .001;

	////robot smarticle geometry
	//double w_smarticle = 0.046; //4.6cm
	//double l_smarticle = 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	//double t_smarticle = .03;
	//double t2_smarticle = .021;

	double collisionEnvelope = .1 * t2_smarticle;
	int global_GUI_value = 0;
	bool bucket_exist = true;


	bool povray_output = true;
	int out_fps = 120;
	const std::string out_dir = "PostProcess";
	const std::string pov_dir_mbd = out_dir + "/povFilesSmarticles";

	ChVector<> bucket_ctr = ChVector<>(0,0,0);
	//ChVector<> Cbucket_interior_halfDim = sizeScale * ChVector<>(.05, .05, .025);
	double bucket_rad = sizeScale*0.034;
	ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, .030);

	
	//ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(.1, .1, .05);
	double bucket_half_thick = sizeScale * .005;

	ChSharedPtr<ChTexture> bucketTexture(new ChTexture());
	
	ChSharedPtr<ChTexture> groundTexture(new ChTexture());
	ChSharedPtr<ChTexture> floorTexture(new ChTexture());

// =============================================================================
#if irrlichtVisualization
	class MyEventReceiver : public IEventReceiver {
	public:

		MyEventReceiver(ChIrrAppInterface* myapp, std::vector<Smarticle*> *mySmarticlesVec) {
			sv = mySmarticlesVec;
			// store pointer applicaiton
			app = myapp;
			// ..add a GUI slider to control friction

			text_SmarticleAmt = app->GetIGUIEnvironment()->addStaticText(L"Layers: 0, Smarticles: 0",
				rect<s32>(850, 65, 1050, 80), true);

			text_Q = app->GetIGUIEnvironment()->addStaticText(L"Press Q to activate GUI1",
				rect<s32>(850, 85, 1050, 100), true);
			text_W = app->GetIGUIEnvironment()->addStaticText(L"Press W to activate GUI2",
				rect<s32>(850, 105, 1050, 120), true);
			text_E = app->GetIGUIEnvironment()->addStaticText(L"Press E to activate GUI3",
				rect<s32>(850, 125, 1050, 140), true);
			text_R = app->GetIGUIEnvironment()->addStaticText(L"Press R to vibrate around current angle",
				rect<s32>(850, 145, 1050, 160), true);
			text_T = app->GetIGUIEnvironment()->addStaticText(L"Press T to vibrate around specified angles",
				rect<s32>(850, 165, 1050, 180), true);
			
			angle1Input = app->GetIGUIEnvironment()->addEditBox(L"75",
				rect<s32>(1050, 165, 1100, 180), true);
			angle2Input = app->GetIGUIEnvironment()->addEditBox(L"75",
				rect<s32>(1100, 165, 1150, 180), true);

			text_Y = app->GetIGUIEnvironment()->addStaticText(L"Press Y to move cylinder away",
				rect<s32>(850, 185, 1050, 200), true);


		}
		
		bool OnEvent(const SEvent& event) {
			// check if user moved the sliders with mouse..
			if (event.EventType == irr::EET_KEY_INPUT_EVENT && !event.KeyInput.PressedDown) {
				switch (event.KeyInput.Key)
				{
				case irr::KEY_KEY_Q:
					if (global_GUI_value != 1)
						global_GUI_value = 1;
					else
						global_GUI_value = 0;
					return true;

				case irr::KEY_KEY_W:
					if (global_GUI_value != 2)
						global_GUI_value = 2;
					else
						global_GUI_value = 0;
					return true;
				case irr::KEY_KEY_E:
					if (global_GUI_value != 3)
						global_GUI_value = 3;
					else
						global_GUI_value = 0;
					return true;
				case irr::KEY_KEY_R: //vibrate around current theta
					if (global_GUI_value != 4)
					{
						double CurrTheta01;
						double CurrTheta12;
						global_GUI_value = 4;
						std::pair<double, double> angPair;
						for (int i = 0; i < sv->size(); i++) //get each particles current theta
						{
							Smarticle* sPtr = sv->at(i);

							MoveType currMoveType = sPtr->moveType;
							std::vector<std::pair<double, double>> *v;



							switch (currMoveType) //TODO fix this and put this in a function inside smarticle class
							{
							case 0:
								v = &sPtr->global;
								break;
							case 1:
								v = &sPtr->gui1;
								break;
							case 2:
								v = &sPtr->gui2;
								break;
							case 3:
								v = &sPtr->gui3;
								break;
							case 4:
								v = &sPtr->vib;
								break;
							case 5:
								v = &sPtr->ot;
								break;
							default:
								v = &sPtr->global;
								break;
							}

							CurrTheta01 = v->at(sPtr->moveTypeIdxs.at(currMoveType)).first;
							CurrTheta12 = v->at(sPtr->moveTypeIdxs.at(currMoveType)).second;
							sPtr->vib.clear();


							angPair.first = CurrTheta01;
							angPair.second = CurrTheta12;
							sPtr->vib.push_back(angPair);

							angPair.first = CurrTheta01 - vibAmp;
							angPair.second = CurrTheta12 - vibAmp;
							sPtr->vib.push_back(angPair);

							angPair.first = CurrTheta01;
							angPair.second = CurrTheta12;
							sPtr->vib.push_back(angPair);

							angPair.first = CurrTheta01 + vibAmp;
							angPair.second = CurrTheta12 + vibAmp;
							sPtr->vib.push_back(angPair);
						}

					}
					else
						global_GUI_value = 0;
					return true;

				case irr::KEY_KEY_T: //TODO vibrate around theta specified in boxes  
					if (global_GUI_value != 5)
					{

						std::pair<double, double> angPair;
						double ang1;
						double ang2;
						global_GUI_value = 5;
						for (int i = 0; i < sv->size(); i++) //get each particles current theta
						{
							Smarticle* sPtr = sv->at(i);

							ang1 = wcstod(angle1Input->getText(), NULL)*CH_C_PI / 180;
							ang2 = wcstod(angle2Input->getText(), NULL)*CH_C_PI / 180;
							sPtr->vib.clear();

							//in case strange values are written
							if (ang2 > CH_C_PI || ang2 < -CH_C_PI)
							{
								global_GUI_value = 0;
								return true;
							}
							if (ang2 > CH_C_PI || ang2 < -CH_C_PI)
							{
								global_GUI_value = 0;
								return true;
							}

							angPair.first = ang1;
							angPair.second = ang2;
							sPtr->vib.push_back(angPair);
							//sPtr->vib.assign(0, angPair);


							angPair.first = ang1 - vibAmp;
							angPair.second = ang2 - vibAmp;
							sPtr->vib.push_back(angPair);
							//sPtr->vib.assign(1, angPair);

							angPair.first = ang1;
							angPair.second = ang2;
							sPtr->vib.push_back(angPair);
							//sPtr->vib.assign(2, angPair);

							angPair.first = ang1 + vibAmp;
							angPair.second = ang2 + vibAmp;
							sPtr->vib.push_back(angPair);
							//sPtr->vib.assign(3, angPair);


						}
					}
					else
						global_GUI_value = 0;
					return true;

				case irr::KEY_KEY_Y:
					if (bucket_exist)
					{
						bucket->SetPos(ChVector<>(100, 0, 0));
						bucket_exist = false;
					}

					return true;

				}


			}
			return false;
		}
		void drawSmarticleAmt(int numLayers)
		{
			char message[100]; sprintf(message, "Layers: %d, Smarticles: %d", numLayers, sv->size());
			this->text_SmarticleAmt->setText(core::stringw(message).c_str());
		}
	private:
		double vibAmp = 2 * CH_C_PI / 180; //vibrate by some amount of degrees back and forth
		std::vector<Smarticle*> *sv;
		ChIrrAppInterface* app;
		IGUIScrollBar* scrollbar_friction;
		IGUIStaticText* text_Q;
		IGUIStaticText* text_W;
		IGUIStaticText* text_E;
		IGUIStaticText* text_R;
		IGUIStaticText* text_T;
		IGUIStaticText* text_Y;
		IGUIStaticText* text_SmarticleAmt;
		IGUIScrollBar* scrollbar_cohesion;
		IGUIStaticText* text_cohesion;
		IGUIScrollBar* scrollbar_compliance;
		IGUIStaticText* text_compliance;
		IGUIStaticText* text_angle1;
		IGUIStaticText* text_angle2;
		IGUIEditBox* angle1Input;
		IGUIEditBox* angle2Input;
		

	};
#endif
// =============================================================================
void MySeed(double s = time(NULL)) { srand(s); }
double MyRand() { return float(rand()) / RAND_MAX; }
// =============================================================================
void SetArgumentsForMbdFromInput(int argc, char* argv[], int& threads, int& max_iteration_sliding, int& max_iteration_bilateral, double& dt, int& num_layers, double& mangle) {
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
		mangle = atof(text);
	}

	/// if parallel, get solver setting
  if (USE_PARALLEL) {
	  if (argc > 5) {
		const char* text = argv[5];
		threads = atoi(text);
	  }
	  if (argc > 6) {
		const char* text = argv[6];
		max_iteration_sliding = atoi(text);
	  }
	  if (argc > 7) {
		const char* text = argv[7];
		max_iteration_bilateral = atoi(text);
	  }
  }
}
// =============================================================================
void InitializeMbdPhysicalSystem_NonParallel(ChSystem& mphysicalSystem, int argc, char* argv[]) {
	// initializd random seeder
	MySeed(923);


  // ---------------------
  // Print the rest of parameters
  // ---------------------

	int dummyNumber0;
	int dummyNumber1;
	int dummyNumber2;
  SetArgumentsForMbdFromInput(argc, argv, dummyNumber0, dummyNumber1, dummyNumber2, dT,numLayers, armAngle);

  simParams << std::endl <<
		  " l_smarticle: " << l_smarticle << std::endl <<
		  " l_smarticle mult for w (w = mult x l): " << l_smarticle /  w_smarticle << std::endl <<
		  " dT: " << dT << std::endl << std::endl;

  // ---------------------
  // Edit mphysicalSystem settings.
  // ---------------------

  // Modify some setting of the physical system for the simulation, if you want
  mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR); // LCP_ITERATIVE_SOR_MULTITHREAD , LCP_ITERATIVE_SOR  (LCP_ITERATIVE_SOR_MULTITHREAD does not work)
  mphysicalSystem.SetIterLCPmaxItersSpeed(50);
  mphysicalSystem.SetIterLCPmaxItersStab(5);   // unuseful for Anitescu, only Tasora uses this
  mphysicalSystem.SetParallelThreadNumber(1);  //TODO figure out if this can increase speed
  mphysicalSystem.SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  mphysicalSystem.SetIterLCPwarmStarting(true);
  mphysicalSystem.SetUseSleeping(false);
  mphysicalSystem.Set_G_acc(ChVector<>(0, 0, gravity));
}
// =============================================================================
void InitializeMbdPhysicalSystem_Parallel(ChSystemParallelDVI& mphysicalSystem, int argc, char* argv[]) {
	// initializd random seeder
	MySeed(964);
  // Desired number of OpenMP threads (will be clamped to maximum available)
  int threads = 1;
  // Perform dynamic tuning of number of threads?
  bool thread_tuning = true;

  //	uint max_iteration = 20;//10000;
  int max_iteration_normal = 50;
  int max_iteration_sliding = 50;
  int max_iteration_spinning = 0;
  int max_iteration_bilateral = 50;

  // ----------------------
  // Set params from input
  // ----------------------

  SetArgumentsForMbdFromInput(argc, argv, threads, max_iteration_sliding, max_iteration_bilateral, dT,numLayers, armAngle);

  // ----------------------
  // Set number of threads.
  // ----------------------

  //  
	int max_threads = omp_get_num_procs();
	if (threads > max_threads)
	threads = max_threads;
	mphysicalSystem.SetParallelThreadNumber(threads);
	omp_set_num_threads(threads);

	mphysicalSystem.GetSettings()->perform_thread_tuning = thread_tuning;
	mphysicalSystem.GetSettings()->min_threads = std::max(1, threads/2);
	mphysicalSystem.GetSettings()->max_threads = std::min(max_threads, int(3.0 * threads / 2));
  // ---------------------
  // Print the rest of parameters
  // ---------------------
	simParams << std::endl <<
		" number of threads: " << threads << std::endl <<
		" max_iteration_normal: " << max_iteration_normal << std::endl <<
		" max_iteration_sliding: " << max_iteration_sliding << std::endl <<
		" max_iteration_spinning: " << max_iteration_spinning << std::endl <<
		" max_iteration_bilateral: " << max_iteration_bilateral << std::endl <<
		" l_smarticle: " << l_smarticle << std::endl <<
		" l_smarticle mult for w (w = mult x l): " << l_smarticle / w_smarticle << std::endl <<
		" dT: " << dT << std::endl <<
		" tFinal: " << tFinal << std::endl <<
		" vibrate start: " << vibrateStart << std::endl <<
		" arm angle: " << armAngle << std::endl << std::endl;
	

  // ---------------------
  // Edit mphysicalSystem settings.
  // ---------------------

  double tolerance = 0.001;  // 1e-3;  // Arman, move it to paramsH
  mphysicalSystem.Set_G_acc(ChVector<>(0, 0, gravity));

  mphysicalSystem.GetSettings()->solver.solver_mode = SLIDING;                              // NORMAL, SPINNING
  mphysicalSystem.GetSettings()->solver.max_iteration_normal = max_iteration_normal;        // max_iteration / 3
  mphysicalSystem.GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;      // max_iteration / 3
  mphysicalSystem.GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;    // 0
  mphysicalSystem.GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;  // max_iteration / 3
  mphysicalSystem.GetSettings()->solver.tolerance = tolerance;
  mphysicalSystem.GetSettings()->solver.alpha = 0;  // Arman, find out what is this
  mphysicalSystem.GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
  mphysicalSystem.ChangeSolverType(APGD);  // Arman check this APGD APGDBLAZE
  //  mphysicalSystem.GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

  mphysicalSystem.GetSettings()->collision.collision_envelope = collisionEnvelope;
  mphysicalSystem.GetSettings()->collision.bins_per_axis = _make_int3(40, 40, 40);  // Arman check
}
#if irrlichtVisualization
void AddParticlesLayer1(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec, ChIrrApp& application) {
#else
void AddParticlesLayer1(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec) {
#endif
	double z;
	int numPerLayer =1;
	int smarticleCount = mySmarticlesVec.size();
	double ang = 2*CH_C_PI / numPerLayer;
	double w = w_smarticle;
	if (smarticleCount < numPerLayer){ z = 0; }
	else{ z = Find_Max_Z(mphysicalSystem); }
	double phase = MyRand()*CH_C_PI / 2;
	for (int i = 0; i < numPerLayer; i++)
	{
		phase = MyRand()*CH_C_PI / 2;
		ChVector<> myPos = bucket_ctr + ChVector<>(sin(ang * i+phase) *(bucket_rad/2 + w*MyRand()-w/2),
			cos(ang*i + phase)*(bucket_rad / 2 + w*MyRand() - w / 2),
			std::min(4 * bucket_interior_halfDim.z, z) + (i + 1)*w_smarticle / 4);
			
		ChQuaternion<> myRot = ChQuaternion<>(MyRand(), MyRand(), MyRand(), MyRand());
			myRot.Normalize();
		
			Smarticle * smarticle0 = new Smarticle(&mphysicalSystem);
			smarticle0->Properties(mySmarticlesVec.size(), mySmarticlesVec.size() * 4,
				rho_smarticle, mat_g,
				collisionEnvelope,
				l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
				sOmega,
				true,
				myPos,
				myRot
				);

			smarticle0->populateMoveVector();

			//TODO figure out why I cannot start at initial position of input file(doesn't move correctly if done)
			smarticle0->SetAngle(0, 0, true);
			smarticle0->Create();
			//smarticle0->AddMotion(myMotionDefault);
			//smarticle0->AddMotion(myMotion);
			mySmarticlesVec.push_back((Smarticle*)smarticle0);
#if irrlichtVisualization
			application.AssetBindAll();
			application.AssetUpdateAll();
#endif
	}


}
// =============================================================================
#if irrlichtVisualization
void AddParticlesLayer(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec, ChIrrApp& application) {
#else
void AddParticlesLayer(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec) {
#endif

	/////////////////
	// Smarticle body
	/////////////////
	ChVector<> smarticleLengths(l_smarticle, w_smarticle, t_smarticle); // l, w, t
	ChVector<> sLenghWithTol = 1.3 * ChVector<>(smarticleLengths.x, smarticleLengths.y, smarticleLengths.z);
	ChVector<> myPos;
	double z;
	double maxDim = 1.3 * std::max(sLenghWithTol.x, sLenghWithTol.y);
	int nX = bucket_interior_halfDim.x / maxDim;
	int nY = bucket_interior_halfDim.y / maxDim;

	int smarticleCount = mySmarticlesVec.size();
	if (smarticleCount < 9){ z = 0; }
	else{ z = Find_Max_Z(mphysicalSystem); }
	int numPerLayer = 4;
	//this filling method works better for smaller diameter where diameter < 3*width of staple	
		//for (int i = -nX + 1; i < nX; i++) {
			//for (int j = -nY + 1; j < nY; j++) {
		for (int i = 0; i < numPerLayer; i++)
		{
				ChQuaternion<> myRot = ChQuaternion<>(MyRand(), MyRand(), MyRand(), MyRand());
				myRot.Normalize();


				//ChVector<> myPos = ChVector<>(i * maxDim + bucket_ctr.x + MyRand()*w_smarticle - 0.5 * w_smarticle
				//	, j * maxDim + bucket_ctr.y + MyRand() * w_smarticle - 0.5 * w_smarticle
				//	, z + maxDim);

				ChVector<> myPos = ChVector<>(bucket_ctr.x + MyRand()* (bucket_interior_halfDim.x - MyRand()*bucket_interior_halfDim.x / 2.0),
					bucket_ctr.y + (MyRand()*bucket_interior_halfDim.y - MyRand()*bucket_interior_halfDim.y / 2.0),
					std::min(4 * bucket_interior_halfDim.z, z) + (i+1)*w_smarticle / 4);

				//ChVector<> myPos = ChVector<>(i * maxDim + bucket_ctr.x, j * maxDim + bucket_ctr.y, bucket_ctr.z + 6.0 * bucket_interior_halfDim.z + 2 * bucket_half_thick);


				//ChVector<> myPos = ChVector<>(i * maxDim, j * maxDim, bucket_ctr.z + 6.0 * bucket_interior_halfDim.z + 2 * bucket_half_thick);
				// ***  added 2*bucket_half_thick to make sure stuff are initialized above bucket. Remember, bucket is inclusive, i.e. the sizes are extende 2*t from each side



//				ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));



				/*ChSharedPtr<SmarticleMotionPiece> myMotionDefault(new SmarticleMotionPiece);
				myMotionDefault->joint_01.theta1 = -.00001;
				myMotionDefault->joint_01.theta2 =  .00001;
				myMotionDefault->joint_01.omega = 0;
				myMotionDefault->joint_12.theta1 = -.00001;
				myMotionDefault->joint_12.theta2 =  .00001;
				myMotionDefault->joint_12.omega = 0;
				myMotionDefault->timeInterval = 0.5;
				myMotionDefault->startTime = 0;
				myMotionDefault->SetMotionType(RELEASE_G);


				ChSharedPtr<SmarticleMotionPiece> myMotion(new SmarticleMotionPiece);
				myMotion->joint_01.theta1 = -0.5 * CH_C_PI;
				myMotion->joint_01.theta2 =  0.5 * CH_C_PI;
				myMotion->joint_01.omega = sOmega;
				myMotion->joint_12.theta1 = -0.5 * CH_C_PI;
				myMotion->joint_12.theta2 =  0.5 * CH_C_PI;
				myMotion->joint_12.omega = sOmega;
				myMotion->timeInterval = 0.5;
				myMotion->startTime = 0;
				myMotion->SetMotionType(SQUARE_G);*/

				
				if (smarticleType == SMART_ARMS) {
					Smarticle * smarticle0 = new Smarticle(&mphysicalSystem);
					smarticle0->Properties(smarticleCount, smarticleCount * 4,
						rho_smarticle, mat_g,
						collisionEnvelope,
						l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
						sOmega,
						true,
						myPos,
						myRot
						);

					smarticle0->populateMoveVector();
					
					//TODO figure out why I cannot start at initial position of input file(doesn't move correctly if done)
					smarticle0->SetAngle(0,0, true);
					smarticle0->Create();
					//smarticle0->AddMotion(myMotionDefault);
					//smarticle0->AddMotion(myMotion);
					mySmarticlesVec.push_back((Smarticle*)smarticle0);
					#if irrlichtVisualization
										application.AssetBindAll();
										application.AssetUpdateAll();
					#endif
				}

				else if (smarticleType == SMART_U) {
					SmarticleU * smarticle0 = new SmarticleU(&mphysicalSystem);
					smarticle0->Properties(smarticleCount,
						rho_smarticle, mat_g,
						collisionEnvelope,
						l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
						myPos,
						myRot);
					smarticle0->Create();
					mySmarticlesVec.push_back(smarticle0);


					#if irrlichtVisualization
										application.AssetBindAll();
										application.AssetUpdateAll();
					#endif
				}
				else {
					std::cout << "Error! Smarticle type is not set correctly" << std::endl;
				}
				smarticleCount++;
			//}
		}
}
// =============================================================================
//creates an approximate cylinder from a n-sided regular polygon
//num_boxes = number of boxes to use
//bucket_rad = radius of cylinder, center point to midpoint of side a side

ChSharedPtr<ChBody> create_cylinder_from_blocks(int num_boxes, int id, bool overlap, CH_SYSTEM* mphysicalSystem, ChSharedPtr<ChMaterialSurfaceBase> wallMat)
{
	ChSharedPtr<ChBody> cyl_container;
	if (USE_PARALLEL) {
		cyl_container = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	}
	else {
		cyl_container = ChSharedPtr<ChBody>(new ChBody);
	}
	cyl_container->SetIdentifier(id);
	//cyl_container->SetMass(mass);
	cyl_container->SetPos(bucket_ctr);
	cyl_container->SetRot(QUNIT);
	cyl_container->SetBodyFixed(false);
	cyl_container->SetCollide(true);
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double wallt = t / 5; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
	double half_height = bucket_interior_halfDim.z;
	double box_side = bucket_rad * 2.0 * tan(CH_C_PI / num_boxes);//side length of cyl
	double o_lap = 0;
	if (overlap){ o_lap = t * 2; }
	double ang = 2.0 * CH_C_PI / num_boxes;
	ChVector<> box_size = (0,0,0); //size of plates
	ChVector<> pPos = (0,0,0);  //position of each plate
	ChQuaternion<> quat=QUNIT; //rotation of each plate
	ChSharedPtr<ChBoxShape> box(new ChBoxShape);
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
		if (ang*i < CH_C_PI  || ang*i > 3.0 * CH_C_PI / 2.0) 
		{
			m_visualization = true;
			cyl_container->AddAsset(bucketTexture);
		}
		cyl_container->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		utils::AddBoxGeometry(cyl_container.get_ptr(), box_size, pPos, quat, m_visualization);
		
	}
	//Add ground piece
	//
	//utils::AddBoxGeometry(cyl_container.get_ptr(), Vector(bucket_rad, bucket_rad + t, t), Vector(0, 0, -t), QUNIT, true);

	//checks top,bottom, and middle location
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos() + Vector(0,0,2 * bucket_interior_halfDim.z + 2 * bucket_half_thick), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos(), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos() + Vector(0, 0, bucket_interior_halfDim.z), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	
	//ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, bucket_interior_halfDim.z);
	
	
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad + 2 * t, t, ChVector<>(0, 0, -t), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//add up volume of bucket and multiply by rho to get mass;
	double cyl_volume = CH_C_PI*(2 * box_size.z - 2 * t)*(2 * box_size.z - 2 * t)*((2 * bucket_rad + 2 * t)*(2 * bucket_rad + 2 * t) - bucket_rad*bucket_rad) + (CH_C_PI)*(bucket_rad + 2 * t)*(bucket_rad + 2 * t) * 2 * t;
	cyl_container->SetMass(rho_cylinder*cyl_volume);

	//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	cyl_container->GetCollisionModel()->BuildModel();

	mphysicalSystem->AddBody(cyl_container);
	return cyl_container;
}

// =============================================================================
//creates an approximate cylinder from a n-sided regular polygon
//num_boxes = number of boxes to use
//bucket_rad = radius of cylinder, center point to midpoint of side a side
ChSharedPtr<ChBody> Create_hopper(CH_SYSTEM* mphysicalSystem, ChSharedPtr<ChMaterialSurfaceBase> wallMat, double w1, double w2, double w3, double h1, double h2,  bool overlap)
{
	ChSharedPtr<ChBody> cyl_container;
	if (USE_PARALLEL) {
		cyl_container = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	}
	else {
		cyl_container = ChSharedPtr<ChBody>(new ChBody);
	}


	double hw1 = 0.5 * w1;
	double hw2 = 0.5 * w2;
	double hw3 = 0.5 * w3;
	double hh1 = 0.5 * h1;
	double hh2 = 0.5 * h2;
	double ht = bucket_half_thick;

	//cyl_container->SetIdentifier(id);
	//cyl_container->SetMass(mass);
	cyl_container->SetPos(bucket_ctr);
	cyl_container->SetRot(QUNIT);
	cyl_container->SetBodyFixed(true);
	cyl_container->SetCollide(true);


	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double o_lap = 0;
	if (overlap){ o_lap = 2 * t; }

	cyl_container->GetCollisionModel()->ClearModel();
	cyl_container->SetMaterialSurface(wallMat);
	double mtheta = atan((hw1 - hw3) / h1);


	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));
	
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(hw1 + ht, 0, h1 + hh2), QUNIT, true); // upper part, max_x plate

	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(-hw1 - ht, 0, h1 + hh2), QUNIT, true); // upper part, min_x plate
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, hw2 + ht, h1 + hh2), QUNIT, true); // upper part, min_x plate
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, -hw2 - ht, h1 + hh2), QUNIT, true); // upper part, min_x plate

	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, -hw2 - ht, hh1), QUNIT, true); // upper part, min_x plate
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, hw2 + ht, hh1), QUNIT, true); // upper part, min_x plate

	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(hw3 + hh1 * tan(mtheta) + ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(mtheta, VECT_Y), true); // upper part, min_x plate
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(-hw3 - hh1 * tan(mtheta) - ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(-mtheta, VECT_Y), true); // upper part, min_x plate
	cyl_container->AddAsset(bucketTexture);

	double estimated_volume = 8 * (w1 * t * h1); // Arman : fix this
	cyl_container->SetMass(rho_cylinder*estimated_volume);
	cyl_container->GetCollisionModel()->BuildModel();
	mphysicalSystem->AddBody(cyl_container);
	return cyl_container;
}

// =============================================================================
void CreateMbdPhysicalSystemObjects(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec) {
	/////////////////
	// Ground body
	/////////////////

	// ground
	ChVector<> boxDim = sizeScale * ChVector<>(0.1, 0.1, .002);
	ChVector<> boxLoc = sizeScale * ChVector<>(0, 0, -5.0*bucket_interior_halfDim.z);
	ChSharedPtr<ChBody> ground;
	if (USE_PARALLEL) {
		ground = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
		bucket = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
		bucket_bott = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	} else {
		ground = ChSharedPtr<ChBody>(new ChBody);
		bucket = ChSharedPtr<ChBody>(new ChBody);
		bucket_bott = ChSharedPtr<ChBody>(new ChBody);
	}
	ground->SetMaterialSurface(mat_g);
	ground->SetPos(boxLoc);

	// ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(true);

	ground->GetCollisionModel()->ClearModel();
	ground->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddCylinderGeometry(ground.get_ptr(), boxDim.x, boxDim.z, ChVector<>(0,0,0), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	ground->GetCollisionModel()->SetFamily(1);
	ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	ground->GetCollisionModel()->BuildModel();
	mphysicalSystem.AddBody(ground);
	
	groundTexture->SetTextureFilename(GetChronoDataFile("greenwhite.png"));
	ground->AddAsset(groundTexture);

	// 1: create bucket
		mat_g->SetFriction(0.4); //steel- plexiglass   (plexiglass was outer cylinder material)
	if (bucketType == BOX){
	//	bucket = utils::CreateBoxContainer(&mphysicalSystem, 1, mat_g, bucket_interior_halfDim, bucket_half_thick, bucket_ctr, QUNIT, true, false, true, false);
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));
		//bucket->AddAsset(mtexture);
		bucket = Create_hopper(&mphysicalSystem, mat_g, bucket_interior_halfDim.x, bucket_interior_halfDim.y, 0.5 * bucket_interior_halfDim.x, bucket_interior_halfDim.z, 2 * bucket_interior_halfDim.z,  true);

	}
	if (bucketType == CYLINDER){
		//http://www.engineeringtoolbox.com/friction-coefficients-d_778.html to get coefficients

		bucket = create_cylinder_from_blocks(25, 1, true, &mphysicalSystem, mat_g);
		

		bucket_bott->SetBodyFixed(true);
		bucket_bott->SetCollide(true);
		bucket_bott->GetCollisionModel()->ClearModel();
		bucket_bott->SetPos(bucket_ctr);
		bucket_bott->SetMaterialSurface(mat_g);
		floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));//custom file
		bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		utils::AddBoxGeometry(bucket_bott.get_ptr(), Vector(bucket_rad + 2 * bucket_half_thick, bucket_rad + 2 * bucket_half_thick, bucket_half_thick), Vector(0, 0, -bucket_half_thick), QUNIT, true);
		bucket_bott->AddAsset(floorTexture);
		bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
		
		bucket_bott->GetCollisionModel()->BuildModel();



	
	
		
	
		mphysicalSystem.AddBody(bucket_bott);
		mat_g->SetFriction(0.5); //steel - steel
	}

	bucket->SetBodyFixed(false);
	bucket->SetCollide(true);
	bucket->GetCollisionModel()->SetFamily(1);
	bucket->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	// 2: create plate
//	bucket->SetMaterialSurface(mat_g);
//	bucket->SetBodyFixed(true);
//	bucket->SetCollide(true);
//	bucket->GetCollisionModel()->ClearModel();
//	utils::AddBoxGeometry(bucket.get_ptr(), ChVector<>(bucket_interior_halfDim.x, bucket_interior_halfDim.y, bucket_half_thick), bucket_ctr - ChVector<>(0, 0, bucket_interior_halfDim.z));
//	bucket->GetCollisionModel()->BuildModel();
//	mphysicalSystem.AddBody(bucket);



//	//	/////////////////
//	// Smarticle body
//	/////////////////
//
//	MySeed(964);
//	ChVector<> smarticleLengths(l_smarticle, w_smarticle, t_smarticle); // l, w, t
//	ChVector<> sLenghWithTol = 1.3 * ChVector<>(smarticleLengths.x, smarticleLengths.y, 2 * smarticleLengths.z);
//
//	double maxDim = 1.3 * std::max(sLenghWithTol.x, sLenghWithTol.y);
//	int nX = bucket_interior_halfDim.x / maxDim;
//	int nY = bucket_interior_halfDim.y / maxDim;
//	// test one smarticle
//
//	ChQuaternion<> myRot = Q_from_AngAxis(CH_C_PI / 2, VECT_X);// ChQuaternion<>(MyRand(), MyRand(), MyRand(), MyRand());
//	myRot.Normalize();
//	ChVector<> myPos = ChVector<>(nX / 2 * maxDim, nY / 2 * maxDim , 2 * maxDim);
//	SmarticleU * smarticle0  = new SmarticleU(&mphysicalSystem);
//	smarticle0->Properties( 3 /* 1 and 2 are the first two objects */,
//					  rho_smarticle, mat_g, l_smarticle, w_smarticle, t_smarticle, t2_smarticle,
//					  myPos,
//					  myRot);
//	smarticle0->Create();
//
//	ChVector<> inertiaS =1e6 * smarticle0->GetSmarticleBodyPointer()->GetInertiaXX();
//	printf("e inertia %f %f %f \n",inertiaS.x, inertiaS.y, inertiaS.z );
//
//	mySmarticlesVec.push_back(smarticle0);
//
//////*** stuff needed to be printed
////	ChVector<> IXX = smarticle0->GetSmarticleBodyPointer()->GetInertiaXX();
////	ChVector<> IXY = smarticle0->GetSmarticleBodyPointer()->GetInertiaXY();
////	ChVector<> posB = smarticle0->GetSmarticleBodyPointer()->GetPos();
////	ChVector<> posS = smarticle0->Get_InitPos();
////	double mass = smarticle0->GetSmarticleBodyPointer()->GetMass();






/////////////////////////////////////////////
//    ChSharedPtr<ChLinkLockRevolute> bucketGroundPrismatic(new ChLinkLockRevolute);
//    bucketGroundPrismatic->Initialize(
//    		smarticle0.GetArm(0), ground, true, ChCoordsys<>(ChVector<>(0, l/2, 0)), ChCoordsys<>(posRel + ChVector<>(0, l/2, 0)));
//    bucketGroundPrismatic->SetName("ship_ground_prismatic");
//    mphysicalSystem.AddLink(bucketGroundPrismatic);


//  ChSharedPtr<ChLinkLockRevolute> bucketGroundPrismatic(new ChLinkLockRevolute);
//  bucketGroundPrismatic->Initialize(
//  		smarticle0.GetArm(0), smarticle0.GetArm(1), true, ChCoordsys<>(ChVector<>(0, l/2, 0)), ChCoordsys<>(posRel + ChVector<>(0, l/2, 0)));
//  bucketGroundPrismatic->SetName("ship_ground_prismatic");
//  mphysicalSystem.AddLink(bucketGroundPrismatic);

}
// =============================================================================

void SavePovFilesMBD(CH_SYSTEM& mphysicalSystem,
                     int tStep) {
  int out_steps = std::ceil((1.0 / dT) / out_fps);
  printf("tStep %d , outstep %d, num bodies %d chrono_time %f\n", tStep, out_steps, mphysicalSystem.Get_bodylist()->size(), mphysicalSystem.GetChTime());

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
double Find_Max_Z(CH_SYSTEM& mphysicalSystem) {
	std::string smarticleTypeName;
	if (smarticleType == SMART_ARMS) {
		smarticleTypeName = "smarticle_arm";
	} else if (smarticleType == SMART_U) {
		smarticleTypeName = "smarticle_u";
	} else {
		std::cout << "Error! Smarticle type is not set correctly" << std::endl;
	}
	double zMax = -999999999;
	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
	for (int i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
		ChBody* bodyPtr = *(myIter + i);
		if ( strcmp(bodyPtr->GetName(), smarticleTypeName.c_str()) == 0 ) {
			if (zMax < bodyPtr->GetPos().z) {
				//zMax = bodyPtr->GetPos().z;
				zMax = bodyPtr->GetPos().z - bucket->GetPos().z;
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
	double xydist = (std::sqrt(dist.x * dist.x + dist.y * dist.y));
	xydist = (std::sqrt((pt.x - centralPt.x)*(pt.x - centralPt.x) + (pt.y - centralPt.y)*(pt.y - centralPt.y)));
//	if ((xydist < rad.x) && (pt.z >= rad.y) && (pt.z < rad.z)) {
//		return true;
//	}
//	return false;
	if (xydist >= rad.x) { return false; } // if outside radius
	if (pt.z < rad.y || pt.z >rad.z){ return false; }
	return true;
}
// =============================================================================
void FixBodies(CH_SYSTEM& mphysicalSystem, int tStep) {
	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
	for (int i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
		ChBody* bodyPtr = *(myIter + i);
		if (bodyPtr->GetPos().z < -1.5 * bucket_interior_halfDim.z) {
			bodyPtr->SetBodyFixed(true);
		}
	}
}
// =============================================================================
void PrintFractions(CH_SYSTEM& mphysicalSystem, int tStep, std::vector<Smarticle*> mySmarticlesVec) {
	const std::string vol_frac = out_dir + "/volumeFraction.txt";
	int stepSave = 10;
	if (tStep % stepSave != 0) return;

	std::ofstream vol_frac_of;
	if (tStep == 0) {
	  vol_frac_of.open(vol_frac.c_str());
	} else {
	  vol_frac_of.open(vol_frac.c_str(), std::ios::app);
	}

	double zMax = Find_Max_Z(mphysicalSystem);
	ChVector<> bucketMin = bucket->GetPos();

	zMax = std::min(zMax, bucketMin.z + 2 * bucket_interior_halfDim.z);
	// *** remember, 2 * bucket_half_thick is needed since bucket is initialized inclusive. the half dims are extended 2*bucket_half_thick from each side

	ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, bucket_interior_halfDim.z);
//	const std::string smarticleTypeName;
//	if (smarticleType == SMART_ARMS) {
//		smarticleTypeName = "smarticle_arm";
//	} else if (smarticleType == SMART_U) {
//		smarticleTypeName = "smarticle_u";
//	} else {
//		std::cout << "Error! Smarticle type is not set correctly" << std::endl;
//	}
//	int countInside = 0;
//	double totalVolume1 = 0;
//	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
//	for (int i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
//		ChBody* bodyPtr = *(myIter + i);
//		if ( strcmp(bodyPtr->GetName(), smarticleTypeName.c_str()) == 0 ) {
//			if ( IsIn(bodyPtr->GetPos(), bucketCtr - bucket_interior_halfDim, bucketCtr + bucket_interior_halfDim + 2 * ChVector<>(0, 0, 2 * bucket_half_thick)) ) {
//				countInside ++;
//				totalVolume1 += bodyPtr->GetMass() / bodyPtr->GetDensity();
//			}
//		}
//	}

	double totalVolume2 = 0;
	int countInside2 = 0;
	double volumeFraction = 0;
	if (bucketType == BOX)
	{
		for (int i = 0; i < mySmarticlesVec.size(); i++) {
			Smarticle* sPtr = mySmarticlesVec[i];
			if (IsIn(sPtr->Get_cm(), bucketCtr - bucket_interior_halfDim, bucketCtr + bucket_interior_halfDim + ChVector<>(0, 0, 2 * bucket_half_thick))) {
				countInside2++;
				totalVolume2 += sPtr->GetVolume();
			}
		}

		volumeFraction = totalVolume2 / (4 * bucket_interior_halfDim.x * bucket_interior_halfDim.y * (zMax - bucketMin.z));
	}
	if (bucketType == CYLINDER)
	{
		for (int i = 0; i < mySmarticlesVec.size(); i++) {
			Smarticle* sPtr = mySmarticlesVec[i];
			//isinradial rad parameter is Vector(bucketrad,zmin,zmax)
			if (IsInRadial(sPtr->Get_cm(), bucketCtr, ChVector<>(bucket_rad, bucketMin.z, bucketMin.z+2*bucket_interior_halfDim.z))) {
				countInside2++;
				totalVolume2 += sPtr->GetVolume();
			}
		}

		volumeFraction = totalVolume2 / (CH_C_PI * bucket_rad * bucket_rad * 2 * bucket_interior_halfDim.z);
	}

	vol_frac_of << mphysicalSystem.GetChTime() << ", " << countInside2  << ", " << volumeFraction << ", " << zMax << std::endl;

	vol_frac_of.close();
}
// =============================================================================
//void SetEnvelopeForSystemObjects(ChSystem& mphysicalSystem) {
//	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
//	for (int i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
//		(*myIter)->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
//		myIter++;
//	}
//
//}
// =============================================================================
// move bucket
void vibrate_bucket(double t) {
	double x_bucket = vibration_amp*sin(omega_bucket * t);
	double xDot_bucket = vibration_amp*omega_bucket*cos(omega_bucket * t);
	double xDDot_bucket = vibration_amp*omega_bucket*omega_bucket*-1 * sin(omega_bucket * t);
	bucket->SetPos(ChVector<>(0, 0, x_bucket));
	bucket->SetPos_dt(ChVector<>(0, 0, xDot_bucket));
	bucket->SetPos_dtdt(ChVector<>(0, 0, xDDot_bucket));
	bucket->SetRot(QUNIT);
}
// =============================================================================
void UpdateSmarticles(
		CH_SYSTEM& mphysicalSystem,
		std::vector<Smarticle*> mySmarticlesVec) {

	double current_time = mphysicalSystem.GetChTime();
	for (int i = 0; i < mySmarticlesVec.size(); i++) {
		//mySmarticlesVec[i]->MoveLoop();

		mySmarticlesVec[i]->MoveLoop2(global_GUI_value);
		//TODO moveloop2(guistate)
//		mySmarticlesVec[i]->UpdateMySmarticleMotion();
//
//		if (current_time > 0.4 && current_time < 0.8) {
//					mySmarticlesVec[i]->MoveToAngle(CH_C_PI/3, -CH_C_PI/3);
//		} else if (current_time > 0.8 && current_time < 1.2) {
//			mySmarticlesVec[i]->MoveToAngle(CH_C_PI/2, CH_C_PI/2);
//		} else if (current_time >= 1.2 && current_time < 2.0) {
//			mySmarticlesVec[i]->MoveToAngle(CH_C_PI/2, -CH_C_PI/2);
//		} else if (current_time >= 2.0 && current_time < 2.4) {
//			mySmarticlesVec[i]->MoveToAngle(0, 0);
//		} else {
//			mySmarticlesVec[i]->MoveToAngle(CH_C_PI/3, CH_C_PI/3);
//		}

	}
}
// =============================================================================
//TODO write gui (dont forget ifopengl in it)


//bool screenshot(char *fileName){
//	int Xres = 1280;
//	int Yres = 720;
//	static unsigned char header[54] = {
//		0x42, 0x4D, 0x36, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x36, 0x00, 0x00, 0x00, 0x28, 0x00,
//		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x18, 0x00, 0x00, 0x00,
//		0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0xC4, 0x0E, 0x00, 0x00, 0xC4, 0x0E, 0x00, 0x00, 0x00, 0x00,
//		0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
//
//	unsigned char *pixels = (unsigned char *)malloc(Xres * Yres * 3);
//	((unsigned __int16 *)header)[9] = Xres;
//	((unsigned __int16 *)header)[11] = Yres;
//
//	glReadPixels(0, 0, Xres, Yres, GL_RGB, GL_UNSIGNED_BYTE, pixels);
//
//	unsigned char temp;
//	for (unsigned int i = 0; i < Xres * Yres * 3; i += 3){
//		temp = pixels[i];
//		pixels[i] = pixels[i + 2];
//		pixels[i + 2] = temp;
//	}
//
//	HANDLE FileHandle;
//	unsigned long Size;
//
//	if (fileName == NULL){
//		char file[256];
//		unsigned int i = 0;
//		do {
//			sprintf(file, "Screenshot%d.bmp", i);
//			FileHandle = CreateFile(file, GENERIC_WRITE, 0, NULL, CREATE_NEW, FILE_ATTRIBUTE_NORMAL, NULL);
//			i++;
//		} while (FileHandle == INVALID_HANDLE_VALUE);
//	}
//	else {
//		FileHandle = CreateFile(fileName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
//		if (FileHandle == INVALID_HANDLE_VALUE)	return false;
//	}
//	DWORD NumberOfBytesWritten;
//	WriteFile(FileHandle, header, sizeof(header), &NumberOfBytesWritten, NULL);
//	WriteFile(FileHandle, pixels, Xres * Yres * 3, &NumberOfBytesWritten, NULL);
//
//	CloseHandle(FileHandle);
//
//	free(pixels);
//	return true;
//}
int main(int argc, char* argv[]) {
	  time_t rawtime;
	  struct tm* timeinfo;
	  time(&rawtime);
	  timeinfo = localtime(&rawtime);
	  ChTimerParallel step_timer;

		//set chrono dataPath to data folder placed in smarticle directory so we can share created files
		#ifdef _WIN64
			std::string fp = "\\..\\data\\";
			fp = __FILE__ + fp;
			SetChronoDataPath(fp);
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

	  const std::string simulationParams = out_dir + "/simulation_specific_parameters.txt";
	  simParams.open(simulationParams.c_str());
	  simParams << " Job was submitted at date/time: " << asctime(timeinfo) << std::endl;


	  // define material property for everything
		mat_g = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
		mat_g->SetFriction(0.5); // .6 for wall to staple using tan (theta) tested on 7/20

  // Create a ChronoENGINE physical system
  CH_SYSTEM mphysicalSystem;
	
#if (USE_PARALLEL)
	  InitializeMbdPhysicalSystem_Parallel(mphysicalSystem, argc, argv);
#else
	  InitializeMbdPhysicalSystem_NonParallel(mphysicalSystem, argc, argv);
#endif

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
	gl_window.viewer->render_camera.camera_scale = 2.0/(1000.0)*sizeScale;
	gl_window.viewer->render_camera.near_clip = .001;
	gl_window.SetRenderMode(opengl::WIREFRAME);
	//TODO study ChOpenGlWindow.cpp to figure out how to make buttons

// Uncomment the following two lines for the OpenGL manager to automatically
// run the simulation in an infinite loop.
// gl_window.StartDrawLoop(time_step);
// return 0;
#endif

#if irrlichtVisualization
  std::cout << "@@@@@@@@@@@@@@@@  irrlicht stuff  @@@@@@@@@@@@@@@@" << std::endl;
  // Create the Irrlicht visualization (open the Irrlicht device,
  // bind a simple user interface, etc. etc.)
	ChIrrApp application(&mphysicalSystem, L"Dynamic Smarticles",
		core::dimension2d<u32>(appWidth, appHeight), false, true);
  
  // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
  ChIrrWizard::add_typical_Logo(application.GetDevice());
  ChIrrWizard::add_typical_Sky(application.GetDevice());
  ChIrrWizard::add_typical_Lights(application.GetDevice(),
                                  core::vector3df(-.1, -.06, .1),
                                  core::vector3df(0, 0, -.01));
	ChIrrWizard::add_typical_Lights(application.GetDevice());


	scene::RTSCamera* camera = new scene::RTSCamera(application.GetDevice(), application.GetDevice()->getSceneManager()->getRootSceneNode(),
		application.GetDevice()->getSceneManager(), -1, -50.0f, 0.5f, 0.0005f);
	camera->setUpVector(core::vector3df(0, 0, 1));//TODO ask arman why up vector isn't changing camera orientation at beginning
	camera->setPosition(core::vector3df(-.1, -.06, .1));
	camera->setTarget(core::vector3df(0, 0, -.01));
	camera->setNearValue(0.01f);
	camera->setMinZoom(0.6f);



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
	mphysicalSystem.SetIterLCPmaxItersSpeed(120);
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
	MyEventReceiver receiver(&application,&mySmarticlesVec);
	// note how to add the custom event receiver to the default interface:
	application.SetUserEventReceiver(&receiver);

#endif


  int stepEnd = int(tFinal / dT);  // 1.0e6;//2.4e6;//600000;//2.4e6 * (.02 * paramsH.sizeScale) /
  // ***************************** Simulation loop ********************************************

  ChSharedPtr<ChFunction_Const> fun1 = ChSharedPtr<ChFunction_Const>(new ChFunction_Const(0));
  ChSharedPtr<ChFunction_Ramp> fun2 = ChSharedPtr<ChFunction_Ramp>(new ChFunction_Ramp(0,5));
  ChSharedPtr<ChFunction_Ramp> fun3 = ChSharedPtr<ChFunction_Ramp>(new ChFunction_Ramp(0,-1));

  ChSharedPtr<ChFunction_Const> fun4 = ChSharedPtr<ChFunction_Const>(new ChFunction_Const(CH_C_PI / 2));
  ChSharedPtr<ChFunction_Const> fun5 = ChSharedPtr<ChFunction_Const>(new ChFunction_Const(-CH_C_PI / 2));


  //AddParticlesLayer(mphysicalSystem, mySmarticlesVec);


//  for (int i = 0; i < mySmarticlesVec.size(); i++) {
//	  mySmarticlesVec[i]->SetActuatorFunction(0, fun2);
//	  mySmarticlesVec[i]->SetActuatorFunction(1, fun1);
//
//  }



  //double timeForVerticalDisplcement = 1.0 * sqrt(2 * w_smarticle / mphysicalSystem.Get_G_acc().Length()); // 1.5 for safety proximity
	//removed length since unnecessary sqrt in that calc

	double timeForVerticalDisplcement = 0.1; // 1.5 for safety proximity
 
	int numGeneratedLayers = 0;

	  //int sSize1 = mySmarticlesVec.size();
	  //if (  (fmod(mphysicalSystem.GetChTime(), timeForVerticalDisplcement) < dT)  &&
			//  (numGeneratedLayers < numLayers) ){
		 // AddParticlesLayer(mphysicalSystem, mySmarticlesVec);
		 // numGeneratedLayers ++;
	  //}

//  CheckPointSmarticles_Read(mphysicalSystem, mySmarticlesVec);

  printf("************** size sys %d \n", mySmarticlesVec.size());
//  for (int tStep = 0; tStep < 1; tStep++) {
	
  for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
	  double t = mphysicalSystem.GetChTime();

	  int sSize1 = mySmarticlesVec.size();
	  if (  (fmod(mphysicalSystem.GetChTime(), timeForVerticalDisplcement) < dT)  &&
			  (numGeneratedLayers < numLayers) ){
#if irrlichtVisualization
			//AddParticlesLayer(mphysicalSystem, mySmarticlesVec,application);
			AddParticlesLayer1(mphysicalSystem, mySmarticlesVec,application);
#else
			//AddParticlesLayer(mphysicalSystem, mySmarticlesVec);
			AddParticlesLayer1(mphysicalSystem, mySmarticlesVec);
#endif 

			

		  numGeneratedLayers ++;
	  }



		//if (numGeneratedLayers == numLayers)
		//{
		//	//start shaking
		//}


//		if (smarticleType == SMART_ARMS)
//		{
//
//			for (int i = 0; i < mySmarticlesVec.size(); i++) {
//				double omega = 10;
//				if (tStep < 500) {
//					mySmarticlesVec[i]->SetActuatorFunction(0, -omega, dT);
//				}
//				else {
//					mySmarticlesVec[i]->SetActuatorFunction(0, omega, dT);
//				}
//
//			}
//
//		}

	   printf("\n");

		 if (t > vibrateStart){
			 bucket->SetBodyFixed(false);
			 vibrate_bucket(t);
		 }
		 else{ bucket->SetBodyFixed(true);}
		 receiver.drawSmarticleAmt(numGeneratedLayers);
		
	



//	  int stage = int(t / (CH_C_PI/2));
//	  printf("yo %d \n", stage%4);
//	  switch (stage % 4) {
//	  case 0: {
//		  smarticle0->SetActuatorFunction(0, fun2);
//	  } break;
//	  case 1: {
////		  smarticle0->SetActuatorFunction(0, fun1);
//	  } break;
//	  case 2: {
//		  smarticle0->SetActuatorFunction(0, fun3);
//	  } break;
//	  case 3: {
////		  smarticle0->SetActuatorFunction(0, fun1);
//	  } break;
//	  }
	  SavePovFilesMBD(mphysicalSystem, tStep);
	  step_timer.start("step time");

#ifdef CHRONO_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(dT);
      gl_window.Render();
    }
#else
#if irrlichtVisualization
    if (!(application.GetDevice()->run())) break;
    application.GetVideoDriver()->beginScene(true, true,
    		video::SColor(255, 140, 161, 192));
//    ChIrrTools::drawGrid(
//        application.GetVideoDriver(), .2, .2, 150, 150,
//        ChCoordsys<>(ChVector<>(0.5 * hdim.x, boxMin.y, 0.5 * hdim.z),
//                     Q_from_AngAxis(CH_C_PI / 2, VECT_X)),
//        video::SColor(50, 90, 90, 150), true);
		//application.AssetBindAll();
		//application.AssetUpdateAll();
		application.DrawAll();
    application.DoStep();
    application.GetVideoDriver()->endScene();
#else
    mphysicalSystem.DoStepDynamics(dT);
#endif
#endif



    UpdateSmarticles(mphysicalSystem, mySmarticlesVec);
	  time(&rawtimeCurrent);
	  double timeDiff = difftime(rawtimeCurrent, rawtime);
		char filename[100];
		sprintf(filename, "screenshot%d.bmp", tStep);
		//screenshot(filename);
	  step_timer.stop("step time");
	  std::cout << "step time: " << step_timer.GetTime("step time") << ", time passed: " << int(timeDiff)/3600 <<":"<< (int(timeDiff) % 3600) / 60 << ":" << (int(timeDiff) % 60) <<std::endl;

	  FixBodies(mphysicalSystem, tStep);
	  PrintFractions(mphysicalSystem, tStep, mySmarticlesVec);

	  std::cout.flush();


//	  CheckPointSmarticles_Write(mySmarticlesVec,
//	  		tStep,
//	  		mat_g,
//	  		l_smarticle,
//	  		w_smarticle,
//	  		t_smarticle,
//	  		t2_smarticle,
//	  		collisionEnvelope,
//	  		rho_smarticle);

  }
  for (int i = 0; i < mySmarticlesVec.size(); i++) {
	  delete mySmarticlesVec[i];

  }
  mySmarticlesVec.clear();
	simParams << "completed"<<std::endl;
  simParams.close();
  return 0;
}


