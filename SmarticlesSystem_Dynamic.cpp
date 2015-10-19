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
#include "chrono_irrlicht/ChBodySceneNode.h"  //changed path from unit to chrono to reflect changes in updated chrono
#include "chrono_irrlicht/ChBodySceneNodeTools.h"
//#include "unit_IRRLICHT/ChIrrTools.h"
#include "chrono_irrlicht/ChIrrWizard.h"
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
enum BucketType { CYLINDER, BOX, HULL,RAMP,HOPPER,DRUM};
SmarticleType smarticleType = SMART_ARMS;//SMART_U;
BucketType bucketType = CYLINDER;
// =============================================================================

class MyBroadPhaseCallback : public collision::ChBroadPhaseCallback {
  public:
    /// Callback used to report 'near enough' pairs of models.
    /// This must be implemented by a child class of ChBroadPhaseCallback.
    /// Return false to skip narrow-phase contact generation for this pair of bodies.
   virtual bool BroadCallback( collision::ChCollisionModel* mmodelA,  ///< pass 1st model
    							collision::ChCollisionModel* mmodelB)   ///< pass 2nd model
		{return (!(abs(mmodelA->GetPhysicsItem()->GetIdentifier() - mmodelB->GetPhysicsItem()->GetIdentifier()) < 3));}
};

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


// =============================================================================
double Find_Max_Z(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mSmartVec);
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
	double dT = 0.0005;//std::min(0.0005, 1.0 / vibration_freq / 200);
	double contact_recovery_speed = 2* sizeScale;
	double tFinal = 6;
	double vibrateStart= .750+100;

	double rho_smarticle = 7850.0 / (sizeScale * sizeScale * sizeScale);
	double rho_cylinder = 1180.0 / (sizeScale * sizeScale * sizeScale);
	ChSharedPtr<ChMaterialSurface> mat_g;
	int numLayers = 100;
	double armAngle = 90;
	double sOmega = 5;  // smarticle omega

	ChSharedPtr<ChLinkEngine> drum_actuator;

		// staple smarticle geometry
	double w_smarticle 	= sizeScale * 0.0117;
	double l_smarticle 	= 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	double t_smarticle 	= sizeScale * .00127;
	double t2_smarticle	= sizeScale * .0005;
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

	bool read_from_file = false;
	bool povray_output = false;
	int out_fps = 120;
	const std::string out_dir = "PostProcess";
	const std::string pov_dir_mbd = out_dir + "/povFilesSmarticles";
	int numPerLayer = 5;
	ChVector<> bucket_ctr = ChVector<>(0,0,0);
	//ChVector<> Cbucket_interior_halfDim = sizeScale * ChVector<>(.05, .05, .025);
	//double bucket_rad = sizeScale*0.034;
	//double bucket_rad = sizeScale*0.02;
	double bucket_rad = sizeScale*0.022;
	ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, .030);


	//ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(.1, .1, .05);
	double bucket_half_thick = sizeScale * .005;

	double max_z = 0;
	double rampAngle = 10 * CH_C_PI / 180;
	double rampInc = 1.0/60.0;
	double drum_freq = 1;
	double drum_omega = drum_freq*2*CH_C_PI;
	double pctActive = 1.0;
	double angle1 = 90;
	double angle2 = 90;
	double vibAmp = 5 * CH_C_PI / 180; //vibrate by some amount of degrees back and forth
	ChSharedPtr<ChTexture> bucketTexture(new ChTexture());
	ChSharedPtr<ChTexture> groundTexture(new ChTexture());
	ChSharedPtr<ChTexture> floorTexture(new ChTexture());

// =============================================================================
#if irrlichtVisualization
	class MyEventReceiver : public IEventReceiver {
	public:

		MyEventReceiver(ChIrrApp* myapp, std::vector<Smarticle*> *mySmarticlesVec) {

			sv = mySmarticlesVec;
			// store pointer applicaiton
			app = myapp;
			// ..add a GUI slider to control friction


			text_Angle = app->GetIGUIEnvironment()->addStaticText(L"Angle: 0, Increment: 0",
				rect<s32>(850, 45, 1050, 60), true);
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
			text_successful = app->GetIGUIEnvironment()->addStaticText(L"Successfully moved smarticles",
				rect<s32>(850, 205, 1050, 220), true);
			resetSuccessfulCount();

		}

		bool OnEvent(const SEvent& event) {
			// check if user moved the sliders with mouse..
			if (event.EventType == irr::EET_KEY_INPUT_EVENT && !event.KeyInput.PressedDown) {
				switch (event.KeyInput.Key)
				{
				case irr::KEY_KEY_Q:
					if (Smarticle::global_GUI_value != 1)
						Smarticle::global_GUI_value = 1;
					else
						Smarticle::global_GUI_value = 0;
					return true;
					break;

				case irr::KEY_KEY_W:
					if (Smarticle::global_GUI_value != 2)
						Smarticle::global_GUI_value = 2;
					else
						Smarticle::global_GUI_value = 0;
					return true;
					break;
				case irr::KEY_KEY_E:
					if (Smarticle::global_GUI_value != 3)
						Smarticle::global_GUI_value = 3;
					else
						Smarticle::global_GUI_value = 0;
					return true;
					break;
				case irr::KEY_KEY_R: //vibrate around current theta
					if (Smarticle::global_GUI_value != 4) //TODO create a boolean and then call a method later which performs this so we don't have to run through smarticle vec here!
					{
						
						double CurrTheta01;
						double CurrTheta12;
						Smarticle::global_GUI_value = 4;
						std::pair<double, double> angPair;
						for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
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

								sPtr->vib.emplace_back(CurrTheta01, CurrTheta12);
								sPtr->vib.emplace_back(CurrTheta01 - vibAmp, CurrTheta12 - vibAmp);
								sPtr->vib.emplace_back(CurrTheta01, CurrTheta12);
								sPtr->vib.emplace_back(CurrTheta01 + vibAmp, CurrTheta12 + vibAmp);
							}
						}

					else
						Smarticle::global_GUI_value = 0;
					return true;
					break;

				case irr::KEY_KEY_T:
					if (Smarticle::global_GUI_value != 5)//TODO create a boolean and then call a method later which performs this so we don't have to run through smarticle vec here!
					{
						std::pair<double, double> angPair;
						double ang1;
						double ang2;
						Smarticle::global_GUI_value = 5;
						for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
						{
							Smarticle* sPtr = sv->at(i);

							ang1 = wcstod(angle1Input->getText(), NULL)*CH_C_PI / 180;
							ang2 = wcstod(angle2Input->getText(), NULL)*CH_C_PI / 180;
							sPtr->vib.clear();

							//in case strange values are written
							if (ang2 > CH_C_PI || ang2 < -CH_C_PI)
							{
								Smarticle::global_GUI_value = 0;
								return true;
								break;
							}
							if (ang2 > CH_C_PI || ang2 < -CH_C_PI)
							{
								Smarticle::global_GUI_value = 0;
								return true;
								break;
							}

							sPtr->vib.emplace_back(ang1, ang2);

							sPtr->vib.emplace_back(ang1 - vibAmp, ang2 - vibAmp);

							sPtr->vib.emplace_back(ang1, ang2);

							sPtr->vib.emplace_back(ang1 + vibAmp, ang2 + vibAmp);

						}
					}
					else
						Smarticle::global_GUI_value = 0;
					return true;
					break;

				case irr::KEY_KEY_Y: //remove conainer or floor
					if (bucket_exist)
					{
						switch (bucketType)
						{
						case CYLINDER:
							bucket->SetPos(ChVector<>(100, 0, 0));
							break;
						case HOPPER:
							bucket_bott->SetPos(ChVector<>(100, 0, 0));
							break;
						case RAMP:
							bucket_bott->SetPos(ChVector<>(100, 0, 0));
							break;
						}

						bucket_exist = false;
					}

					return true;
					break;


				case irr::KEY_KEY_1:
					switch (bucketType)
					{
					case DRUM:
						drum_freq = drum_freq - rampInc;
						break;
					case RAMP:

						rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0;
						bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
							Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0
							, 0, 0)));
						bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
							Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0
							, 0, 0)));
						break;
					default:
						rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0;
						bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
							Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x - rampInc * CH_C_PI / 180.0
							, 0, 0)));
						break;
					}
					drawAngle();
					return true;
					break;
				case irr::KEY_KEY_2:			//increase angle of bucket by rampInc
					switch (bucketType)
					{
					case DRUM:
						drum_freq = drum_freq + rampInc;
						break;
					case RAMP:

						rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0;
						bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
							Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0
							, 0, 0)));
						bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
							Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0
							, 0, 0)));
						break;
					default:
						rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0;
						bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
							Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x + rampInc * CH_C_PI / 180.0
							, 0, 0)));
						break;
					}
					drawAngle();
					return true;
					break;
				case irr::KEY_KEY_3:			//decrease rampInc
					rampInc = rampInc-1.0/60.0;
					drawAngle();
					return true;
					break;
				case irr::KEY_KEY_4:			//increase rampInc
					rampInc = rampInc + 1.0/60.0;
					drawAngle();
					return true;
					break;
				}


			}
			return false;
		}
		void drawSmarticleAmt(int numLayers)//nu
		{
			char message[100]; sprintf(message, "Layers: %d, Smarticles: %d, GUI: %d", numLayers, sv->size(),Smarticle::global_GUI_value);
			this->text_SmarticleAmt->setText(core::stringw(message).c_str());
		}
		void drawAngle()
		{

			if (bucketType == DRUM)
			{
				char message[100]; sprintf(message, "AngVel: %g rpm, Increment: %g", drum_freq*60, rampInc*60);
				this->text_Angle->setText(core::stringw(message).c_str());
			}
			else{
				char message[100]; sprintf(message, "Angle: %1.1g, Increment: %1.3g", Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x * 180 / CH_C_PI, rampInc);
				this->text_Angle->setText(core::stringw(message).c_str());
			}

		}
		void drawSuccessful()
		{
			int count = 0;
			for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
			{
				if(sv->at(i)->successfulMotion) count++;
			}
			char message[100]; sprintf(message, "Successfully Moving: %d/%d", count, sv->size());
			this->text_successful->setText(core::stringw(message).c_str());
		}
		void drawSuccessful2()
		{
			char message[100]; sprintf(message, "Successfully Moving: %d/%d", successfulCount, sv->size());
			this->text_successful->setText(core::stringw(message).c_str());
		}
		void addSuccessful(Smarticle &sPtr)
		{
			if (sPtr.successfulMotion)
				successfulCount++;

		}
		void drawOTArms(Smarticle &sPtr)
		{
			if (sPtr.GetArm0OT())
			{
				sPtr.GetArm(0)->AddAsset(Smarticle::mtextureOT);
			}
			else
			{
				sPtr.GetArm(0)->AddAsset(Smarticle::mtextureArm);
			}
			if (sPtr.GetArm2OT())
			{
				sPtr.GetArm(2)->AddAsset(Smarticle::mtextureOT);
			}
			else
			{
				sPtr.GetArm(2)->AddAsset(Smarticle::mtextureArm);
			}

			app->AssetBind(sPtr.GetArm(0));
			app->AssetBind(sPtr.GetArm(2));
			app->AssetUpdate(sPtr.GetArm(0));
			app->AssetUpdate(sPtr.GetArm(2));
		}
		void drawOTArms()
		{
			for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
			{
				Smarticle* sPtr = sv->at(i);

				if (sPtr->GetArm0OT())
				{
					sPtr->GetArm(0)->AddAsset(Smarticle::mtextureOT);
				}
				else
				{
					sPtr->GetArm(0)->AddAsset(Smarticle::mtextureArm);
				}
				if (sPtr->GetArm2OT())
				{
					sPtr->GetArm(2)->AddAsset(Smarticle::mtextureOT);
				}
				else
				{
					sPtr->GetArm(2)->AddAsset(Smarticle::mtextureArm);
				}

				app->AssetBind(sPtr->GetArm(0));
				app->AssetBind(sPtr->GetArm(2));
				app->AssetUpdate(sPtr->GetArm(0));
				app->AssetUpdate(sPtr->GetArm(2));
			}
		}
		void resetSuccessfulCount()
		{
			successfulCount = 0;
		}
	private:

		std::vector<Smarticle*> *sv;
		ChIrrApp* app;
		IGUIScrollBar* scrollbar_friction;
		IGUIStaticText* text_Q;
		IGUIStaticText* text_W;
		IGUIStaticText* text_E;
		IGUIStaticText* text_R;
		IGUIStaticText* text_T;
		IGUIStaticText* text_Y;
		IGUIStaticText* text_SmarticleAmt;
		IGUIStaticText* text_Angle;
		IGUIScrollBar* scrollbar_cohesion;
		IGUIStaticText* text_cohesion;
		IGUIScrollBar* scrollbar_compliance;
		IGUIStaticText* text_compliance;
		IGUIStaticText* text_angle1;
		IGUIStaticText* text_angle2;
		IGUIEditBox* angle1Input;
		IGUIEditBox* angle2Input;
		IGUIStaticText* text_successful;
		static int successfulCount;


	};
	int MyEventReceiver::successfulCount = 0;
#endif
// =============================================================================
void MySeed(double s = time(NULL)) { srand(s); }
double MyRand() { return float(rand()) / RAND_MAX; }
// =============================================================================
void SetArgumentsForMbdFromInput(int argc, char* argv[], int& threads, int& max_iteration_sliding, int& max_iteration_bilateral, double& dt, int& num_layers, double& mangle,bool& readFile,double& mpctActive, double& mangle1, double& mangle2) {
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
	}
	if (argc > 7){
		const char* text = argv[7];
		angle2 = atof(text);
	}
	/// if parallel, get solver setting
  if (USE_PARALLEL) {
	  if (argc > 8) {
		const char* text = argv[8];
		threads = atoi(text);
	  }
	  if (argc > 9) {
		const char* text = argv[9];
		max_iteration_sliding = atoi(text);
	  }
	  if (argc > 10) {
		const char* text = argv[10];
		max_iteration_bilateral = atoi(text);
	  }
  }
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
	if (read_from_file)
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
	mphysicalSystem.SetIterLCPmaxItersSpeed(int(2.85*numLayers*numPerLayer));
  mphysicalSystem.SetIterLCPmaxItersStab(0);   // unuseful for Anitescu, only Tasora uses this
  mphysicalSystem.SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  mphysicalSystem.SetIterLCPwarmStarting(true);
  mphysicalSystem.SetUseSleeping(false);
  mphysicalSystem.Set_G_acc(ChVector<>(0, 0, gravity));

	//mphysicalSystem.SetTolForce(.0005);
	//mphysicalSystem.SetTol(.0001);
	//mphysicalSystem.SetMinBounceSpeed(.3);
	simParams.close();
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

  SetArgumentsForMbdFromInput(argc, argv, threads, max_iteration_sliding, max_iteration_bilateral, dT,numLayers, armAngle,read_from_file,pctActive,angle1,angle2);

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
	const std::string simulationParams = out_dir + "/simulation_specific_parameters.txt";
	simParams.open(simulationParams.c_str(), std::ios::app);
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
		" read from file: " << read_from_file << std::endl <<
		" arm angle: " << angle1<< " " << angle2 << std::endl << std::endl;


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
void AddParticlesLayer1(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec, ChIrrApp& application,double timeForDisp) {
#else
void AddParticlesLayer1(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec,double timeForDisp) {
#endif


	double z;
	double zpos;
	int smarticleCount = mySmarticlesVec.size();
	double ang = 2*CH_C_PI / numPerLayer;
	double w = w_smarticle;
	if (smarticleCount < numPerLayer){ z = w_smarticle / 2; }
	//else{ z = max_z; }
	else{ z = Find_Max_Z(mphysicalSystem,mySmarticlesVec); }
	double phase = MyRand()*CH_C_PI / 2;
	for (int i = 0; i < numPerLayer; i++)
	{
		phase = MyRand()*CH_C_PI / 2 - MyRand()*CH_C_PI / 4 ;
		zpos = std::min(3 * bucket_interior_halfDim.z, z) + w_smarticle / 2;
		if (bucketType==DRUM)
			zpos = std::min(bucket_interior_halfDim.z/2, z) + w_smarticle / 4;

		ChVector<> myPos = bucket_ctr + ChVector<>(sin(ang * i + phase) *(bucket_rad / 2 + w*MyRand()-w/2), //TODO for hopper no -w/2.0
			cos(ang*i + phase)*(bucket_rad / 2 + w*MyRand() - w/2.0),
			zpos);
		if (bucketType == CYLINDER)
		{
				myPos = bucket_ctr + ChVector<>(sin(ang * i + phase) *(bucket_rad / 2.2), //TODO for hopper no -w/2.0
				cos(ang*i + phase)*(bucket_rad / 2.2),
				zpos);

				//ChVector<> myPos = bucket_ctr + ChVector<>((MyRand()-.5)*(bucket_rad / 2), //TODO for hopper no -w/2.0
				//	(MyRand() - .5)*(bucket_rad / 2),
				//	zpos);
		}

			//std::min(4 * bucket_interior_halfDim.z, z) + (i + 1)*w_smarticle / 4);
			//std::min(3 * bucket_interior_halfDim.z, z) + w_smarticle / 4);
			// z + w_smarticle / 2);

		ChQuaternion<> myRot = ChQuaternion<>(2*MyRand()-1, 2*MyRand()-1, 2*MyRand()-1, 2*MyRand()-1);
			myRot.Normalize();

			Smarticle * smarticle0 = new Smarticle(&mphysicalSystem);
			smarticle0->Properties(mySmarticlesVec.size(), mySmarticlesVec.size() * 4,
				rho_smarticle, mat_g,
				collisionEnvelope,
				//l_smarticle+t2_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
				l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
				sOmega,
				true,
				myPos,
				myRot
				);

			if (MyRand()<1)
				smarticle0->visualize = true;
			smarticle0->populateMoveVector();
			smarticle0->SetAngle(angle1, angle2, true);
			smarticle0->Create();

			smarticle0->vib.emplace_back(angle1*CH_C_PI / 180, angle2*CH_C_PI / 180);
			smarticle0->vib.emplace_back(angle1*CH_C_PI / 180 - vibAmp, angle1*CH_C_PI / 180 - vibAmp);
			smarticle0->vib.emplace_back(angle1*CH_C_PI / 180, angle2*CH_C_PI / 180);
			smarticle0->vib.emplace_back(angle1*CH_C_PI / 180 + vibAmp, angle1*CH_C_PI / 180 + vibAmp);

			mySmarticlesVec.emplace_back((Smarticle*)smarticle0);
			smarticle0->SetSpeed(ChVector<>(0, 0, -9.8*timeForDisp / 2.0 - w_smarticle / timeForDisp));
#if irrlichtVisualization
			application.AssetBindAll();
			application.AssetUpdateAll();
#endif
	}


}

ChSharedPtr<ChBody> create_drum(int num_boxes, int id, bool overlap, CH_SYSTEM* mphysicalSystem, ChSharedPtr<ChMaterialSurfaceBase> wallMat,int ridges = 5)
{//essentially the same as create cyl container except made it bigger and added ridges
	ChSharedPtr<ChBody> drum;
	if (USE_PARALLEL) {
		drum = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	}
	else {
		drum = ChSharedPtr<ChBody>(new ChBody);
	}
	double radMult =1.5;
	drum->SetIdentifier(id);
	drum->SetPos(bucket_ctr);
	drum->SetRot(QUNIT);
	drum->SetBodyFixed(false);
	drum->SetCollide(true);
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double wallt = t / 5.0; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
	double half_height = bucket_interior_halfDim.z/(radMult*2);
	double box_side = bucket_rad * radMult*2 * tan(CH_C_PI / num_boxes);//side length of cyl
	double o_lap = 0;
	if (overlap){ o_lap = t * 2; }
	double ang = 2.0 * CH_C_PI / num_boxes;
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
		utils::AddBoxGeometry(drum.get_ptr(), box_size, pPos, quat);



		if (i%ridgeNum == 0)
		{
			ridge_size = ChVector<>((box_side) / 4.0,
				w_smarticle/8.0,
				half_height + o_lap);
			pPos = bucket_ctr + ChVector<>(sin(ang * i) * (-wallt + bucket_rad*radMult),
				cos(ang*i)*(-wallt + bucket_rad*radMult),
				0);

			drum->GetCollisionModel()->SetEnvelope(collisionEnvelope);
			//utils::AddBoxGeometry(drum.get_ptr(), ridge_size, pPos, quat);
		}
		drum->SetRot(Q_from_AngAxis(CH_C_PI / 2.0, VECT_X));


	}
//TODO add bucketVolume as global variable and set it in each function to calculate for each shape volumefraction seamlessly
	//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);

	drum->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	//front wall made invisible so we can see inside
	utils::AddBoxGeometry(drum.get_ptr(), ChVector<>(wallt + bucket_rad*radMult, wallt + bucket_rad*radMult, wallt), bucket_ctr + VECT_Z*(half_height + 2 * o_lap - 2 * t), QUNIT,false);
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
	utils::AddBoxGeometry(bucket_bott.get_ptr(), ChVector<>(wallt + bucket_rad*radMult, wallt + bucket_rad*radMult, wallt), bucket_ctr - VECT_Z*(half_height + 2 * o_lap - 2 * t), QUNIT);


	bucket_bott->GetCollisionModel()->BuildModel();
	mphysicalSystem->AddBody(bucket_bott);
	bucket_bott->SetRot(Q_from_AngAxis(CH_C_PI / 2.0, VECT_X));
	bucket_bott->GetCollisionModel()->SetFamily(1);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	return drum;
}
ChSharedPtr<ChBody> create_complex_convex_hull(CH_SYSTEM* mphysicalSystem, ChSharedPtr<ChMaterialSurfaceBase> wallMat, double numBoxes) //TODO finish this method, currently not being used
{
	ChSharedPtr<ChBody> convexShape;
	if (USE_PARALLEL) { convexShape = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel)); }
	else{ convexShape = ChSharedPtr<ChBody>(new ChBody); }
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code

	//cyl_container->SetMass(mass);
	convexShape->SetPos(bucket_ctr);
	convexShape->SetRot(QUNIT);
	convexShape->SetBodyFixed(false);
	convexShape->SetCollide(true);

	std::vector<ChVector<>> points;

	convexShape->GetCollisionModel()->ClearModel();
	double ang = 2 * CH_C_PI / numBoxes;
	for (size_t i = 0; i < numBoxes; i++)
	{
		convexShape->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	}

	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));

	ChVector<> rampSize(w_smarticle * 5, w_smarticle * 5, t);


	//utils::AddBoxGeometry(convexShape.get_ptr(), rampSize, rampPos, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(rampAngle, 0, 0)), true);



	convexShape->GetCollisionModel()->BuildModel();
	convexShape->AddAsset(floorTexture);
	mphysicalSystem->AddBody(convexShape);
	return convexShape;
}
ChSharedPtr<ChBody> create_ramp(int id, CH_SYSTEM* mphysicalSystem, ChSharedPtr<ChMaterialSurfaceBase> wallMat)
{
	ChSharedPtr<ChBody> ramp;
	if (USE_PARALLEL) {ramp = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));}
	else{ramp = ChSharedPtr<ChBody>(new ChBody);}
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

	ChVector<> rampPos(0, 0, -sin(rampAngle)*w - t);

	ramp->GetCollisionModel()->ClearModel();
	ramp->SetMaterialSurface(wallMat);



	ramp->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	ramp->GetCollisionModel()->SetEnvelope(collisionEnvelope);

	utils::AddBoxGeometry(ramp.get_ptr(), ChVector<>(w, w, t), rampPos, QUNIT, true);//bottom
	utils::AddBoxGeometry(ramp.get_ptr(), ChVector<>(t, w, h), rampPos - VECT_X*(w - t) + VECT_Z*(h - t), QUNIT, true);// -x
	utils::AddBoxGeometry(ramp.get_ptr(), ChVector<>(t, w, h), rampPos + VECT_X*(w - t) + VECT_Z*(h - t), QUNIT, true);// +x
	utils::AddBoxGeometry(ramp.get_ptr(), ChVector<>(w, t, h), rampPos - VECT_Y*(w - t) + VECT_Z*(h - t), QUNIT, true);//high side -y
	//utils::AddBoxGeometry(ramp.get_ptr(), ChVector<>(w, t, h), rampPos + VECT_Y*(w - t) + VECT_Z*(h - t), QUNIT, true); //down side +y


	///bucket_bott in ramp is the +y side!
	bucket_bott->SetBodyFixed(true);
	bucket_bott->SetCollide(true);
	bucket_bott->GetCollisionModel()->ClearModel();
	bucket_bott->SetPos(bucket_ctr);
	bucket_bott->SetMaterialSurface(mat_g);
	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));//custom file
	bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(bucket_bott.get_ptr(), ChVector<>(w, t, h), rampPos + VECT_Y*(w - t) + VECT_Z*(h - t), QUNIT, true); //down side +y
	bucket_bott->AddAsset(floorTexture);
	bucket_bott->GetCollisionModel()->SetFamily(1);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	bucket_bott->GetCollisionModel()->BuildModel();
	bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(rampAngle, 0, 0)));
	mphysicalSystem->AddBody(bucket_bott);
	bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(rampAngle, 0, 0)));

	ramp->GetCollisionModel()->BuildModel();
	ramp->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(rampAngle, 0, 0)));
	ramp->AddAsset(floorTexture);
	mphysicalSystem->AddBody(ramp);
	return ramp;
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
		if (ang*i < 3*CH_C_PI/4  || ang*i > 5* CH_C_PI/4)
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
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos() + Vector(0,0,2 * bucket_interior_halfDim.z + 2.0 * bucket_half_thick), Q_from_AngAxis(CH_C_PI/2.0, VECT_Y),true);
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos(), Q_from_AngAxis(CH_C_PI / 2.0, VECT_Y));
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
ChSharedPtr<ChBody> Create_hopper2(CH_SYSTEM* mphysicalSystem, ChSharedPtr<ChMaterialSurfaceBase> wallMat, double theta, double holeSize, bool overlap)
{
	ChSharedPtr<ChBody> hopper;
	if (USE_PARALLEL) {
		hopper = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	}
	else {
		hopper = ChSharedPtr<ChBody>(new ChBody);
	}
	holeSize = holeSize / 2;
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double r = bucket_rad;
	double ang = theta*CH_C_PI / 180.0;
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
	//utils::AddBoxGeometry(ramp.get_ptr(), ChVector<>(w, w, t), rampPos, QUNIT, true);//bucket_bottom


	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(t, r + t, h), hop_pos - VECT_X*(sin(ang)*h + holeSize+cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, ang, 0)), true);// -x negAngle
	//utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(t, r + t, h), hop_pos - VECT_X*(sin(ang)*(h)+holeSize+cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, ang, 0)), true);// -x negAngle +theta
	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(t, r + t, h), hop_pos + VECT_X*(sin(ang)*h + holeSize+cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, -ang, 0)), true);// +x  -theta = /


	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(2 * sin(ang)*h + holeSize+2*cos(ang)*t, t, h), hop_pos - VECT_Y*(r + o_lap) + VECT_Z*(h - o_lap), QUNIT, false);//camera pos -y
	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(2 * sin(ang)*h + holeSize+2*cos(ang)*t, t, h), hop_pos + VECT_Y*(r + o_lap) + VECT_Z*(h - o_lap), QUNIT, true); //camera target +y

	hopper->GetCollisionModel()->BuildModel();
	mphysicalSystem->AddBody(hopper);
	//utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(r, t, h + o_lap), ChVector<>(-r/2, 0, h+o_lap), QUNIT, true); // front plate, max_x plate

	//utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(-hw1 - ht, 0, h1 + hh2), QUNIT, true); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, hw2 + ht, h1 + hh2), QUNIT, true); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, -hw2 - ht, h1 + hh2), QUNIT, false); // upper part, min_x plate

	//utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, -hw2 - ht, hh1), QUNIT, false); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, hw2 + ht, hh1), QUNIT, true); // upper part, min_x plate

	//utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(hw3 + hh1 * tan(mtheta) + ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(mtheta, VECT_Y), true); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(-hw3 - hh1 * tan(mtheta) - ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(-mtheta, VECT_Y), true); // upper part, min_x plate
	//hopper->AddAsset(bucketTexture);

	//double estimated_volume = 8 * (w1 * t * h1); // Arman : fix this
	//hopper->SetMass(rho_cylinder*estimated_volume);
	//hopper->GetCollisionModel()->BuildModel();
	//mphysicalSystem->AddBody(hopper);
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
	utils::AddBoxGeometry(bucket_bott.get_ptr(), Vector(bucket_rad + 2 * bucket_half_thick, bucket_rad + 2 * bucket_half_thick, bucket_half_thick), Vector(0, 0, -bucket_half_thick), QUNIT, true);
	bucket_bott->AddAsset(floorTexture);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	bucket_bott->GetCollisionModel()->BuildModel();
	mphysicalSystem.AddBody(bucket_bott);
}
ChSharedPtr<ChBody> Create_hopper(CH_SYSTEM* mphysicalSystem, ChSharedPtr<ChMaterialSurfaceBase> wallMat, double w1, double w2, double w3, double h1, double h2,  bool overlap)
{
	ChSharedPtr<ChBody> hopper;
	if (USE_PARALLEL) {
		hopper = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	}
	else {
		hopper = ChSharedPtr<ChBody>(new ChBody);
	}


	double hw1 = w1;
	double hw2 = w2;
	double hw3 = w3;
	double hh1 = h1*.5;
	double hh2 = h2*.5;
	double ht = bucket_half_thick;

	hopper->SetPos(bucket_ctr);
	hopper->SetRot(QUNIT);
	hopper->SetBodyFixed(true);
	hopper->SetCollide(true);


	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double o_lap = 0;
	if (overlap){ o_lap = 2 * t; }

	hopper->GetCollisionModel()->ClearModel();
	hopper->SetMaterialSurface(wallMat);
	double mtheta = atan((hw1 - hw3) / h1);


	bucketTexture->SetTextureFilename(GetChronoDataFile("greenwhite.png"));
	hopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(hw1 + ht, 0, h1 + hh2), QUNIT, true); // upper part, max_x plate

	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(-hw1 - ht, 0, h1 + hh2), QUNIT, true); // upper part, min_x plate
	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, hw2 + ht, h1 + hh2), QUNIT, true); // upper part, min_x plate
	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, -hw2 - ht, h1 + hh2), QUNIT, false); // upper part, min_x plate

	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, -hw2 - ht, hh1), QUNIT, false); // upper part, min_x plate
	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, hw2 + ht, hh1), QUNIT, true); // upper part, min_x plate

	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(hw3 + hh1 * tan(mtheta) + ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(mtheta, VECT_Y), true); // upper part, min_x plate
	utils::AddBoxGeometry(hopper.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(-hw3 - hh1 * tan(mtheta) - ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(-mtheta, VECT_Y), true); // upper part, min_x plate
	hopper->AddAsset(bucketTexture);

	double estimated_volume = 8 * (w1 * t * h1); // Arman : fix this
	hopper->SetMass(rho_cylinder*estimated_volume);
	hopper->GetCollisionModel()->BuildModel();
	mphysicalSystem->AddBody(hopper);
	return hopper;
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
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));
		switch (bucketType)		//http://www.engineeringtoolbox.com/friction-coefficients-d_778.html to get coefficients
		{
		case BOX:
			bucket = utils::CreateBoxContainer(&mphysicalSystem, 1, mat_g, bucket_interior_halfDim, bucket_half_thick, bucket_ctr, QUNIT, true, false, true, false);
			bucket->AddAsset(bucketTexture);
			break;

		case CYLINDER:
			bucket = create_cylinder_from_blocks(25, 1, true, &mphysicalSystem, mat_g);
			CreateBucket_bott(mphysicalSystem);
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
		mat_g->SetFriction(0.5); //steel - steel


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
	for (size_t i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
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
	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
	for (size_t i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
		ChBody* bodyPtr = *(myIter + i);
		if (strcmp(bodyPtr->GetName(), smarticleTypeName.c_str()) == 0) {
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
void drawGlobalCoordinateFrame(CH_SYSTEM& mphysicalSystem,	///< the chrono::engine physical system
	double len = w_smarticle,
	double rad = t_smarticle,
	ChVector<> pos = bucket_ctr + ChVector<>(2.5*bucket_rad,0,bucket_interior_halfDim.z))
{
	ChSharedPtr<ChBody> xaxis, yaxis, zaxis;
	if (USE_PARALLEL) {
		xaxis = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
		yaxis = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
		zaxis = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	}
	else{
		xaxis = ChSharedPtr<ChBody>(new ChBody);
		yaxis = ChSharedPtr<ChBody>(new ChBody);
		zaxis = ChSharedPtr<ChBody>(new ChBody);
	}

	xaxis->SetPos(pos);						yaxis->SetPos(pos);							zaxis->SetPos(pos);
	xaxis->SetCollide(false);			yaxis->SetCollide(false);				zaxis->SetCollide(false);
	xaxis->SetBodyFixed(true);		yaxis->SetBodyFixed(true);			zaxis->SetBodyFixed(true);
	xaxis->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	yaxis->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	zaxis->GetCollisionModel()->SetEnvelope(collisionEnvelope);

	utils::AddCylinderGeometry(xaxis.get_ptr(), rad, len, ChVector<>(len-rad/2, 0, 0) + pos, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, CH_C_PI / 2)), true);//bottom
	utils::AddCylinderGeometry(yaxis.get_ptr(), rad, len, ChVector<>(0, len-rad/2, 0) + pos, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, 0)), true);//bottom, true);//bottom
	utils::AddCylinderGeometry(zaxis.get_ptr(), rad, len, ChVector<>(0, 0, len-rad) + pos, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(CH_C_PI / 2, 0, 0)), true);//bottom

	xaxis->AddAsset(ChSharedPtr<ChColorAsset>(new ChColorAsset(1.0f, 0, 0)));
	yaxis->AddAsset(ChSharedPtr<ChColorAsset>(new ChColorAsset(0, 1.0f, 0)));
	zaxis->AddAsset(ChSharedPtr<ChColorAsset>(new ChColorAsset(0, 0, 1.0f)));
	mphysicalSystem.AddBody(xaxis); mphysicalSystem.AddBody(yaxis); mphysicalSystem.AddBody(zaxis);
}


void recycleSmarticles(CH_SYSTEM& mphysicalSystem, Smarticle& sPtr)
{
	double pos = -.75*bucket_interior_halfDim.z;//z position below which smarticles are regenerated above pile inside container
	double ang = 2 * CH_C_PI / 5;
	double rp = MyRand()*ang / 4; //add slight offset to angInc to allow particles not always fall in nearly same position
	static int recycledSmarticles = 0;
	static int inc = 0;
	//ChVector<> myPos = bucket_ctr + ChVector<>(sin(ang * i + phase) *(bucket_rad / 2 + w*MyRand()), //TODO for hopper no -w/2.0
	//	cos(ang*i + phase)*(bucket_rad / 2 + w*MyRand() - w / 2.0),
	//	zpos);
		if (sPtr.GetArm(1)->GetPos().z < pos)
		{
			if (bucketType == HOPPER)
			{
				sPtr.TransportSmarticle(bucket_ctr + ChVector<>(
					sin(ang*inc + rp)*(bucket_rad / 2 + 4 * w_smarticle*(MyRand() - 1 / 2.0)),
					cos(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(MyRand() - 1 / 2.0)),
					bucket_interior_halfDim.z * 2
					));

				//sPtr->SetSpeed(sPtr->GetArm(1)->GetPos_dt() / 4);
				sPtr.SetSpeed(ChVector<>(0, 0, -9.8*.01 / 2.0 - w_smarticle / .01));
			}
			else
			{

				sPtr.TransportSmarticle(bucket_ctr + ChVector<>(
					sin(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(MyRand() - 1 / 2.0)),
					cos(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(MyRand() - 1 / 2.0)),
					bucket_interior_halfDim.z*1.75
					));
				//sPtr->TransportSmarticle(ChVector<>
				//	(ChVector<>(sPtr->GetArm(1)->GetPos().x,
				//	sPtr->GetArm(1)->GetPos().y,
				//	bucket_interior_halfDim.z*1.75)));
				//sPtr->SetSpeed(sPtr->GetArm(1)->GetPos_dt()/2);
			}


			recycledSmarticles++;
			inc = (inc + 1) % 5;
		}
	printFlowRate(mphysicalSystem.GetChTime(), recycledSmarticles);
}
void recycleSmarticles(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mySmarticlesVec)
{
	double pos = -.75*bucket_interior_halfDim.z;//z position below which smarticles are regenerated above pile inside container
	double ang = 2 * CH_C_PI / 5;
	double rp = MyRand()*ang/4 ; //add slight offset to angInc to allow particles not always fall in nearly same position
	static int recycledSmarticles = 0;
	static int inc = 0;
	//ChVector<> myPos = bucket_ctr + ChVector<>(sin(ang * i + phase) *(bucket_rad / 2 + w*MyRand()), //TODO for hopper no -w/2.0
	//	cos(ang*i + phase)*(bucket_rad / 2 + w*MyRand() - w / 2.0),
	//	zpos);
	for (size_t i = 0; i < mySmarticlesVec.size(); i++)
	{
		Smarticle* sPtr = mySmarticlesVec[i];
		if (sPtr->GetArm(1)->GetPos().z < pos)
		{


			if (bucketType == HOPPER)
			{
				sPtr->TransportSmarticle(bucket_ctr + ChVector<>(
					sin(ang*inc + rp)*(bucket_rad / 2 + 4*w_smarticle*(MyRand() - 1/2.0)),
					cos(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(MyRand() - 1 / 2.0)),
					bucket_interior_halfDim.z*2
					));

				//sPtr->SetSpeed(sPtr->GetArm(1)->GetPos_dt() / 4);
				sPtr->SetSpeed(ChVector<>(0, 0, -9.8*.01 / 2.0 - w_smarticle / .01));
			}
			else
			{

				sPtr->TransportSmarticle(bucket_ctr + ChVector<>(
					sin(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(MyRand() - 1 / 2.0)),
					cos(ang*inc + rp)*(bucket_rad / 2 + w_smarticle*(MyRand() - 1 / 2.0)),
					bucket_interior_halfDim.z*1.75
					));
				//sPtr->TransportSmarticle(ChVector<>
				//	(ChVector<>(sPtr->GetArm(1)->GetPos().x,
				//	sPtr->GetArm(1)->GetPos().y,
				//	bucket_interior_halfDim.z*1.75)));
				//sPtr->SetSpeed(sPtr->GetArm(1)->GetPos_dt()/2);
			}


			recycledSmarticles++;
			inc = (inc+1)%5;
		}
	}
	printFlowRate(mphysicalSystem.GetChTime(), recycledSmarticles);
}
// =============================================================================
void FixBodies(CH_SYSTEM& mphysicalSystem, int tStep) {
	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
	for (size_t i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
		ChBody* bodyPtr = *(myIter + i);
		if (bodyPtr->GetPos().z < -5.0*bucket_interior_halfDim.z) {
			bodyPtr->SetBodyFixed(true);
		}
	}
}
// =============================================================================
bool FixSmarticles(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mySmarticlesVec, Smarticle &sPtr, int tstep, int idx) //true if removed, false if not fixed
{
	//recycle smarticles if bucket or hopper
	if (bucketType == HOPPER && bucket_exist == false) //if hopper, put smarticles back inside after reaching below hopper if bucket_bott still exists delete
	{
		recycleSmarticles(mphysicalSystem, sPtr);
		return false;
	}
	if (sPtr.armBroken)
	{
		sPtr.~Smarticle();
		mySmarticlesVec.erase(mySmarticlesVec.begin()+idx);
		GetLog() << "\narm broken removing smarticle\n";
		return true;
	}
	if (bucketType != HOPPER)
	{
		if (sPtr.GetArm(1)->GetPos().z < -3.0*bucket_interior_halfDim.z) //if far below bucket
		{
			sPtr.~Smarticle();
			mySmarticlesVec.erase(mySmarticlesVec.begin() + idx);
			GetLog() << "\nbelow bucket\n";
			return true;
		}
		if (bucketType == CYLINDER && !IsInRadial(sPtr.Get_cm(), bucket->GetPos() + ChVector<>(0, 0, bucket_interior_halfDim.z), ChVector<>(bucket_rad, bucket->GetPos().z, bucket->GetPos().z + 2 * bucket_interior_halfDim.z))) //if outside radius
		{
			//sPtr->~Smarticle();
			//mySmarticlesVec.erase(mySmarticlesVec.begin() + idx);
			GetLog() << "\noutside radius removing HERE!\n";
			return true;
		}
	}
	else
	{

		if (!IsInRadial(sPtr.Get_cm(), bucket->GetPos(), ChVector<>(2 * bucket_rad, -4.0*bucket_interior_halfDim.z, 4.0*bucket_interior_halfDim.z)))
		{
			sPtr.~Smarticle();
			mySmarticlesVec.erase(mySmarticlesVec.begin() + idx);
			GetLog() << "\noutside radius\n";
			return true;
		}
	}
	return false;

}
void FixSmarticles(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> &mySmarticlesVec, double tstep) { ///remove all traces of smarticle from system
	if (bucketType == HOPPER && bucket_exist==false) //if hopper, put smarticles back inside after reaching below hopper if bucket_bott still exists delete
	{
		recycleSmarticles(mphysicalSystem,mySmarticlesVec);
	}

	std::vector<Smarticle*>::iterator myIter;
	//for(myIter = mySmarticlesVec.begin(); myIter != mySmarticlesVec.end(); ++myIter)
		for (myIter = mySmarticlesVec.begin(); myIter != mySmarticlesVec.end();)
	{
		Smarticle* sPtr = *(myIter);
		if (sPtr->armBroken)
		{
			sPtr->~Smarticle();
			myIter = mySmarticlesVec.erase(myIter);
			GetLog() << "\narm broken removing smarticle\n";
			continue;
		}
		if (sPtr->GetArm(0)->GetRot_dt().GetVector().Length2() > 10000)
		{
			sPtr->GetArm(0)->SetRot_dt(sPtr->GetArm(0)->GetRot() / 2); //added this because small arms can start to spin uncontrollably
			sPtr->GetArm(1)->SetRot_dt(sPtr->GetArm(1)->GetRot() / 2);
			sPtr->GetArm(2)->SetRot_dt(sPtr->GetArm(2)->GetRot() / 2);
		}
		if (bucketType != HOPPER)
		{
			if (sPtr->GetArm(1)->GetPos().z < -3.0*bucket_interior_halfDim.z)
			{
				sPtr->~Smarticle();
				myIter = mySmarticlesVec.erase(myIter);
				continue;
			}
			if (bucketType == CYLINDER && !IsInRadial(sPtr->Get_cm(), bucket->GetPos() + ChVector<>(0, 0, bucket_interior_halfDim.z), ChVector<>(bucket_rad, bucket->GetPos().z, bucket->GetPos().z + 2 * bucket_interior_halfDim.z)))
			{
				sPtr->~Smarticle();
				myIter = mySmarticlesVec.erase(myIter);
				continue;
			}
			 ++myIter;
		}
		else
		{

			if (!IsInRadial(sPtr->GetArm(1)->GetPos(), bucket->GetPos(), ChVector<>(2*bucket_rad, -4.0*bucket_interior_halfDim.z, 4.0*bucket_interior_halfDim.z)))
			{
				sPtr->~Smarticle();
				myIter = mySmarticlesVec.erase(myIter);
			}
			else{ ++myIter; }

		}
	}


}
void PrintFractions(CH_SYSTEM& mphysicalSystem, int tStep, std::vector<Smarticle*> mySmarticlesVec) {
	const static std::string vol_frac = out_dir + "/volumeFraction.txt";
	static int stepSave = 10;
	if (tStep % stepSave != 0) return;
	double zComz = 0;
	double meanOT = 0;
	//static ChSharedPtr<ChBody> grid;  //uncomment to visualize vol frac boxes
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
	ChVector<> bucketMin = bucket->GetPos();


	// *** remember, 2 * bucket_half_thick is needed since bucket is initialized inclusive. the half dims are extended 2*bucket_half_thick from each side

	ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, bucket_interior_halfDim.z);
	double totalVolume2 = 0;
	int countInside2 = 0;
	double volumeFraction = 0;
	if (bucketType == BOX)
	{
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
	}
	if (bucketType == CYLINDER)
	{
		//if (a)
		//	mphysicalSystem.RemoveBody(grid);			//uncomment to visualize vol frac boxes
		//grid = ChSharedPtr<ChBody>(new ChBody); //uncomment to visualize vol frac boxes
		//grid->RecomputeCollisionModel();				//uncomment to visualize vol frac boxes
		//grid->GetCollisionModel()->ClearModel();//uncomment to visualize vol frac boxes
		for (size_t i = 0; i < mySmarticlesVec.size(); i++) {
			Smarticle* sPtr = mySmarticlesVec[i];
			//isinradial rad parameter is Vector(bucketrad,zmin,zmax)
			if (IsInRadial(sPtr->Get_cm(), bucketCtr, ChVector<>(bucket_rad, bucketMin.z, bucketMin.z+2.0*bucket_interior_halfDim.z))) {
				countInside2++;
				//com = sPtr->Get_cm()-ChVector<>(0,0,bucketMin.z);
				com = sPtr->Get_cm() - ChVector<>(0, 0, bucketMin.z);
				zComz += com.z;
				max2 = std::max(max2, com.z);
				if(max2>zMax)
				{
					double temp = zMax;
					zMax = max2;
					max2 = temp;
				}

				//pos = com;
				//xpos = int((pos.x+bucket_rad)/sqSizex);
				//ypos = int((pos.y+bucket_rad)/sqSizey);
				//vecPos = (colSize)*ypos+xpos;
				//zHeights[vecPos]=std::pair<int,double>(zHeights[vecPos].first+1,std::max(pos.z,zHeights[vecPos].second));
				
				//utils::AddBoxGeometry(grid.get_ptr(), ChVector<>(sqSizex / 2, sqSizey / 2, pos.z / 2), ChVector<>((4)*sqSizex + sqSizex / 2 - bucket_rad, (6)*sqSizey + sqSizey / 2 - bucket_rad, bucket->GetPos().z + pos.z / 2), QUNIT, true); // upper part, min_x plate
				//GetLog() << "\n (smartx,smarty): (" << com.x << "," << com.y << ")\n";
				//GetLog() << "\n(xpos,ypos): (" << (xpos)*sqSizex - bucket_rad << "," << (ypos)*sqSizex - bucket_rad << ")\n";
				meanOT += sPtr->GetReactTorqueLen01() + sPtr->GetReactTorqueLen12();
				//zMax = std::max(zMax, sPtr->GetArm(1)->GetPos().z- bucketMin.z);

			}
		}
		
		//uncomment to visualize vol frac boxes
		//grid->AddAsset(Smarticle::mtextureOT);
		//grid->GetCollisionModel()->BuildModel();
		//a = true;
		//mphysicalSystem.AddBody(grid);

		
		//for (size_t i = 0; i < zHeights.size(); i++)
		//{
		//	volumeFraction += zHeights[i].second*sqSizey*sqSizex;
		//}
		//volumeFraction = countInside2*vol/volumeFraction;
	
		volumeFraction = countInside2*vol / (max2*CH_C_PI*bucket_rad*bucket_rad);
		//GetLog() << vol << " " << countInside2 << " " << bucket_rad << " " << zMax << " " << volumeFraction << "\n";
		GetLog() <<"phi="<< volumeFraction << "\n";
		zComz = zComz / countInside2;
		meanOT = meanOT / (countInside2 * 2.0); //multiply by 2 (2 arms for each smarticle)
	}

	vol_frac_of << mphysicalSystem.GetChTime() << ", " << countInside2 << ", " << volumeFraction << ", " << zMax << ", " << zComz << ", " << meanOT<< ", " << Smarticle::global_GUI_value << std::endl;
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
	if (bucketType == CYLINDER)
	{
		double phase = -omega_bucket*vibrateStart;
		double x_bucket = vibration_amp*sin(omega_bucket * t + phase);
		double xDot_bucket = vibration_amp*omega_bucket*cos(omega_bucket * t + phase);
		double xDDot_bucket = vibration_amp*omega_bucket*omega_bucket*-1 * sin(omega_bucket * t + phase);
		bucket->SetPos(ChVector<>(0, 0, x_bucket));
		bucket->SetPos_dt(ChVector<>(0, 0, xDot_bucket));
		bucket->SetPos_dtdt(ChVector<>(0, 0, xDDot_bucket));
		bucket->SetRot(QUNIT);
		bucket_bott->SetPos(ChVector<>(0, 0, x_bucket));
		bucket_bott->SetPos_dt(ChVector<>(0, 0, xDot_bucket));
		bucket_bott->SetPos_dtdt(ChVector<>(0, 0, xDDot_bucket));
		bucket_bott->SetRot(QUNIT);
	}
}
void rotate_drum(double t)//method is called on each iteration to rotate drum at an angular velocity of drum_omega
{

	////bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()) + ChVector<>(0, 0, bucket_omega*t)));
	//bucket->SetRot_dt(Angle_to_Quat(ANGLESET_RXYZ,ChVector<>(0,0,bucket_omega)));
	//bucket->SetRot_dtdt(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, 0)));
	bucket->SetBodyFixed(false);
	bucket_bott->SetBodyFixed(true);


	static ChSharedPtr<ChFunction_Const> mfun2 = drum_actuator->Get_spe_funct().DynamicCastTo<ChFunction_Const>();
	drum_omega = drum_freq*CH_C_PI * 2;
	mfun2->Set_yconst(drum_omega);
	//ChQuaternion<> rspeed = (bucket->GetRot()*Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, drum_omega)));
	//rspeed.Normalize();
	//bucket->SetRot_dt(rspeed);
	//bucket->SetPos(VNULL);

}
void setUpDrumActuator(CH_SYSTEM& mphysicalSystem)
{
	drum_actuator = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
	ChVector<> pR01(0, 0, 0);
	ChQuaternion<> qx = Q_from_AngAxis(CH_C_PI / 2.0, VECT_Z);
	drum_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), bucket->GetRot()));
	drum_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
	mphysicalSystem.AddLink(drum_actuator);
	ChSharedPtr<ChFunction_Const> mfun2 = drum_actuator->Get_spe_funct().DynamicCastTo<ChFunction_Const>();
	drum_omega = drum_freq*CH_C_PI * 2;
	mfun2->Set_yconst(drum_omega);
}

// =============================================================================
void UpdateSmarticles(
	CH_SYSTEM& mphysicalSystem,
	Smarticle& sPtr) {

	double current_time = mphysicalSystem.GetChTime();
	sPtr.MoveLoop2(Smarticle::global_GUI_value);
}
void UpdateSmarticles(
		CH_SYSTEM& mphysicalSystem,
		std::vector<Smarticle*> mySmarticlesVec) {

	double current_time = mphysicalSystem.GetChTime();
	for (size_t i = 0; i < mySmarticlesVec.size(); i++) {
		mySmarticlesVec[i]->MoveLoop2(Smarticle::global_GUI_value);
	}
}
// =============================================================================

int main(int argc, char* argv[]) {
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	ChTimerParallel step_timer;
	Smarticle::global_GUI_value = 0;
	//set chrono dataPath to data folder placed in smarticle directory so we can share created files
#if defined(_WIN64)
	std::string fp = "D:\\ChronoCode\\chronoPkgs\\Smarticles\\data\\";
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

	const std::string simulationParams = out_dir + "/simulation_specific_parameters.txt";
	simParams.open(simulationParams.c_str());
	simParams << "Job was submitted at date/time: " << asctime(timeinfo) << std::endl;
	simParams.close();
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
	GetLog() << "\npctActive" << pctActive<<"\n";
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
		core::vector3df(-.1, 0, .1),
		core::vector3df(-.1, 0, 0));
	ChIrrWizard::add_typical_Lights(application.GetDevice());
	ChIrrWizard::add_typical_Lights(application.GetDevice(),
		core::vector3df(0, -.1, 0),
		core::vector3df(0, 0, -.01));

	scene::RTSCamera* camera = new scene::RTSCamera(application.GetDevice(), application.GetDevice()->getSceneManager()->getRootSceneNode(),
		application.GetDevice()->getSceneManager(), -1, -50.0f, 0.5f, 0.0005f);
	camera->setUpVector(core::vector3df(0, 0, 1));
	camera->setPosition(core::vector3df(0, -.1, 0));
	camera->setTarget(core::vector3df(0, 0, -.01));
	camera->setNearValue(0.0005f);
	camera->setMinZoom(0.1f);
	camera->setZoomSpeed(0.1f);

	drawGlobalCoordinateFrame(mphysicalSystem);


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
	MyEventReceiver receiver(&application, &mySmarticlesVec);

	// note how to add the custom event receiver to the default interface:
	application.SetUserEventReceiver(&receiver);
	receiver.drawAngle();//initialize draw angle

	//ool OnEvent(const SEvent& event) {
	//	// check if user moved the sliders with mouse..
	//	if (event.EventType == irr::EET_KEY_INPUT_EVENT && !event.KeyInput.PressedDown) {
	//		switch (event.KeyInput.Key)
	//		{
	//		case irr::KEY_KEY_Q:
	//			if (Smarticle::global_GUI_value != 1)


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

	double timeForVerticalDisplcement = 0.01; // 1.5 for safety proximity

	int numGeneratedLayers = 0;

	  //int sSize1 = mySmarticlesVec.size();
	  //if (  (fmod(mphysicalSystem.GetChTime(), timeForVerticalDisplcement) < dT)  &&
			//  (numGeneratedLayers < numLayers) ){
		 // AddParticlesLayer(mphysicalSystem, mySmarticlesVec);
		 // numGeneratedLayers ++;
	  //}
	if (bucketType==DRUM)
		setUpDrumActuator(mphysicalSystem);

	if (read_from_file)
	{
		CheckPointSmarticlesDynamic_Read(mphysicalSystem, mySmarticlesVec);
		application.AssetBindAll();
		application.AssetUpdateAll();
		numGeneratedLayers = numLayers;
	}
//  for (int tStep = 0; tStep < 1; tStep++) {
	Smarticle::global_GUI_value = 1;

	for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
		double t = mphysicalSystem.GetChTime();
		if (!read_from_file)
		{
			if ((fmod(t, timeForVerticalDisplcement) < dT) &&
				(numGeneratedLayers < numLayers)){
#if irrlichtVisualization
				AddParticlesLayer1(mphysicalSystem, mySmarticlesVec, application, timeForVerticalDisplcement);
#else
				AddParticlesLayer1(mphysicalSystem, mySmarticlesVec);
#endif
				numGeneratedLayers++;
			}
		}
		if (t > vibrateStart && t<vibrateStart+3){
			//if (t > vibrateStart || bucketType==HOPPER){
			bucket->SetBodyFixed(false);
			vibrate_bucket(t);
		}
		else{ bucket->SetBodyFixed(true); }

		if (bucketType == DRUM)
		{
			rotate_drum(t);
		}
		size_t vecSize = mySmarticlesVec.size();
		////////////////////////////////////////////////////////////////////////////////////
		//receiver.resetSuccessfulCount();
		//bool removedSmart = false;
		//max_z = 0;//reset max_z
		//for (int i = vecSize - 1; i >= 0; --i)
		//{
		//	Smarticle& sPtr = *mySmarticlesVec.at(i);
		//	if (FixSmarticles(mphysicalSystem, mySmarticlesVec, sPtr, tStep, i)) //if removed smarticle
		//		continue;
		//	UpdateSmarticles(mphysicalSystem, sPtr);

		//	i == 0 ? max_z = PrintFractions(mphysicalSystem, tStep, sPtr, true) : max_z = PrintFractions(mphysicalSystem, tStep, sPtr, false);
		//
		//receiver.addSuccessful(sPtr);
		//}
		//receiver.drawSuccessful2();
		//receiver.drawSmarticleAmt(numGeneratedLayers);

		////////FixSmarticles(mphysicalSystem, mySmarticlesVec, tStep);
		////////UpdateSmarticles(mphysicalSystem, mySmarticlesVec);
		////////PrintFractions(mphysicalSystem, tStep, mySmarticlesVec);
		////////receiver.drawOTArms();

		//////





		/////////////////////////////////////////////////////////////////////////////////////



		if (fmod(t, timeForVerticalDisplcement) < dT
			&&mySmarticlesVec.size()< numPerLayer*numLayers && (numGeneratedLayers == numLayers))
			AddParticlesLayer1(mphysicalSystem, mySmarticlesVec, application, timeForVerticalDisplcement);

	  //SavePovFilesMBD(mphysicalSystem, tStep);
	  //step_timer.start("step time");
		//receiver.drawOTArms();
		///
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

		//framerecord
		application.SetVideoframeSaveInterval(20);//only save every 2 frames
		application.DrawAll();
		application.AssetBindAll();  //uncomment to visualize vol frac boxes
		application.AssetUpdateAll();//uncomment to visualize vol frac boxes
    application.DoStep();//
    application.GetVideoDriver()->endScene();
	#else
			mphysicalSystem.DoStepDynamics(dT);
	#endif
#endif

		//Smarticle::global_GUI_value = 1;
		//if (t < .75)
		//	Smarticle::global_GUI_value = 1;
		//else if (t > .75 && t < 1.5)
		//	Smarticle::global_GUI_value = 2;
		//else if (t > 1.5 && t < 2.25)
		//	Smarticle::global_GUI_value = 3;
		//else if (t > 2.25 && t < 4.6)
		//	break;
		//else if (t > 4.6 &&t < 5.6)
		//	Smarticle::global_GUI_value = 1;
		//else
		//	break;


			//if (t < vibrateStart+3.1)
			//	Smarticle::global_GUI_value = 1;
			//else if (t > 3.1)
			//	Smarticle::global_GUI_value = 4;
			//else
			//break;

			//if (bucket_exist && t>.1)
			//{
			//	bucket_bott->SetPos(ChVector<>(100, 0, 0));
			//	bucket_exist = false;
			//}

		receiver.drawSuccessful();
		FixSmarticles(mphysicalSystem, mySmarticlesVec, tStep);
		UpdateSmarticles(mphysicalSystem, mySmarticlesVec);
	  time(&rawtimeCurrent);
	  double timeDiff = difftime(rawtimeCurrent, rawtime);
	  step_timer.stop("step time");
		//max_z = 0;
		PrintFractions(mphysicalSystem, tStep, mySmarticlesVec);
	  std::cout.flush();
		receiver.drawSmarticleAmt(numGeneratedLayers);
		CheckPointSmarticlesDynamic_Write(mySmarticlesVec,
	  		tStep,
	  		mat_g,
	  		l_smarticle,
	  		w_smarticle,
	  		t_smarticle,
	  		t2_smarticle,
	  		collisionEnvelope,
	  		rho_smarticle);







  }
	simParams.open(simulationParams.c_str(), std::ios::app);
	simParams << "Smarticle OT: " <<	 mySmarticlesVec.at(0)->torqueThresh2 << std::endl;
  for (int i = 0; i < mySmarticlesVec.size(); i++) {
	  delete mySmarticlesVec[i];

  }
  mySmarticlesVec.clear();

	simParams << "completed"<<std::endl;
  simParams.close();
  return 0;
}



double PrintFractions(CH_SYSTEM& mphysicalSystem, int tStep, Smarticle &sPtr,bool lastIdx) { //TODO make zmax proper!
	const static std::string vol_frac = out_dir + "/volumeFraction.txt";
	static int stepSave = 10;
	static double zCom = 0;
	static double meanOT = 0;
	static double zMax = 0;
	static double totalVolume2 = 0;
	static int countInside2 = 0;
	static double volumeFraction = 0;
	ChVector<> bucketMin = bucket->GetPos();



	// *** remember, 2 * bucket_half_thick is needed since bucket is initialized inclusive. the half dims are extended 2*bucket_half_thick from each side

	ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, bucket_interior_halfDim.z);

	if (bucketType == CYLINDER)
	{
			//isinradial rad parameter is Vector(bucketrad,zmin,zmax)
			if (IsInRadial(sPtr.Get_cm(), bucketCtr, ChVector<>(bucket_rad, bucketMin.z, bucketMin.z + 2.0*bucket_interior_halfDim.z))) {
				countInside2++;
				totalVolume2 += sPtr.GetVolume();
				zCom += sPtr.Get_cm().z + bucketMin.z;
				meanOT += sPtr.GetReactTorqueLen01() + sPtr.GetReactTorqueLen12();
				zMax = std::max(zMax, sPtr.GetArm(1)->GetPos().z);
			}
	}

	if (lastIdx&& tStep % stepSave == 0)
	{

		std::ofstream vol_frac_of;
		if (tStep == 0) {
			vol_frac_of.open(vol_frac.c_str());
		}
		else {
			vol_frac_of.open(vol_frac.c_str(), std::ios::app);
		}
		volumeFraction = totalVolume2 / (CH_C_PI * bucket_rad * bucket_rad * 2.0 * bucket_interior_halfDim.z);
		zCom = zCom / countInside2;
		meanOT = meanOT / (countInside2 * 2.0); //multiply by 2 (2 arms for each smarticle)
		vol_frac_of << mphysicalSystem.GetChTime() << ", " << countInside2 << ", " << volumeFraction << ", " << zMax << ", " << zCom << ", " << meanOT << ", " << Smarticle::global_GUI_value << std::endl;
		vol_frac_of.close();
	}
	return zMax;
}
