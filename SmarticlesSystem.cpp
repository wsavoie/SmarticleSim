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
//     - collisions and contacts
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
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"     //Arman: why is this
#include "chrono_utils/ChUtilsInputOutput.h"  //Arman: Why is this
#include "chrono_utils/ChUtilsGenerators.h"

#include <ctime>
#include <stdlib.h>  // system, rand, srand, RAND_MAX
#include "core/ChFileutils.h" // for MakeDirectory
#include "Smarticle.h"


//#include <fstream> // Arman: I don't know why this is not required
//#include <cstring> // Arman: I don't know why this is not required
//#include <stdio.h> // Arman: I don't know why this is not required



//#ifdef CHRONO_PARALLEL_HAS_OPENGL
//#undef CHRONO_PARALLEL_HAS_OPENGL
//#endif

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif
//***********************************
// Use the namespace of Chrono

using namespace chrono;
//using namespace chrono::collision;
//using namespace std;

// =============================================================================
std::ofstream simParams;
ChSharedPtr<ChBody> bucket;

double gravity = -9.81;
double dT = .01;
double contact_recovery_speed = .3;
double tFinal = 30;
double sizeScale = 1000;
double rho_smarticle = 7850 / (sizeScale * sizeScale * sizeScale);


bool povray_output = true;
int out_fps = 25;
const std::string out_dir = "PostProcess";
const std::string pov_dir_mbd = out_dir + "/povFilesSmarticles";

ChVector<> bucket_ctr = ChVector<>(0,0,0);
ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(.05, .05, .025);
double bucket_thick = sizeScale * .001;

// =============================================================================
void MySeed(double s = time(NULL)) { srand(s); }
double MyRand() { return float(rand()) / RAND_MAX; }
// =============================================================================
void SetArgumentsForMbdFromInput(int argc, char* argv[], int& threads, int& max_iteration_sliding, int& max_iteration_bilateral) {
  if (argc > 1) {
    const char* text = argv[1];
    threads = atoi(text);
  }
  if (argc > 2) {
    const char* text = argv[2];
    max_iteration_sliding = atoi(text);
  }
  if (argc > 3) {
    const char* text = argv[3];
    max_iteration_bilateral = atoi(text);
  }
}
// =============================================================================
void InitializeMbdPhysicalSystem(ChSystemParallelDVI& mphysicalSystem, int argc, char* argv[]) {
  // Desired number of OpenMP threads (will be clamped to maximum available)
  int threads = 1;
  // Perform dynamic tuning of number of threads?
  bool thread_tuning = true;

  //	uint max_iteration = 20;//10000;
  int max_iteration_normal = 200;
  int max_iteration_sliding = 200;
  int max_iteration_spinning = 0;
  int max_iteration_bilateral = 100;

  // ----------------------
  // Set params from input
  // ----------------------

  SetArgumentsForMbdFromInput(argc, argv, threads, max_iteration_sliding, max_iteration_bilateral);

  // ----------------------
  // Set number of threads.
  // ----------------------

  //  omp_get_num_procs();
  int max_threads = mphysicalSystem.GetParallelThreadNumber();
  if (threads > max_threads)
    threads = max_threads;
  mphysicalSystem.SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);

  mphysicalSystem.GetSettings()->perform_thread_tuning = thread_tuning;
  mphysicalSystem.GetSettings()->min_threads = std::max(1, threads/2);
  mphysicalSystem.GetSettings()->max_threads = int(3.0 * threads / 2);

  // ---------------------
  // Print the rest of parameters
  // ---------------------
  simParams << std::endl <<
		  " number of threads: " << threads << std::endl <<
		  " max_iteration_normal: " << max_iteration_normal << std::endl <<
		  " max_iteration_sliding: " << max_iteration_sliding << std::endl <<
		  " max_iteration_spinning: " << max_iteration_spinning << std::endl <<
		  " max_iteration_bilateral: " << max_iteration_bilateral << std::endl << std::endl;

  // ---------------------
  // Edit mphysicalSystem settings.
  // ---------------------

  double tolerance = 0.1;  // 1e-3;  // Arman, move it to paramsH
  //double collisionEnvelop = 0.04 * paramsH.HSML;
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

//    mphysicalSystem.GetSettings()->collision.collision_envelope = collisionEnvelop;   // global collisionEnvelop does not work. Maybe due to sph-tire size mismatch
  mphysicalSystem.GetSettings()->collision.bins_per_axis = _make_int3(40, 40, 40);  // Arman check
}

// =============================================================================

void CreateMbdPhysicalSystemObjects(ChSystemParallelDVI& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec) {
	  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
		mat_g->SetFriction(0.1);
		mat_g->SetCohesion(0);
		mat_g->SetCompliance(0.0);
		mat_g->SetComplianceT(0.0);
		mat_g->SetDampingF(0.2);

	/////////////////
	// Ground body
	/////////////////

	// ground
	ChVector<> boxDim = sizeScale * ChVector<>(0.1, 0.1, .002);
	ChVector<> boxLoc = sizeScale * ChVector<>(0, 0, -.0027);
	ChSharedPtr<ChBody> ground = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	ground->SetMaterialSurface(mat_g);
	ground->SetPos(boxLoc);

	//  ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(true);

	ground->GetCollisionModel()->ClearModel();
	utils::AddCylinderGeometry(ground.get_ptr(), boxDim.x, boxDim.z, ChVector<>(0,0,0), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	ground->GetCollisionModel()->BuildModel();
	mphysicalSystem.AddBody(ground);

	// bucket
	bucket = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	bucket = utils::CreateBoxContainer(&mphysicalSystem, 1, mat_g, bucket_interior_halfDim, bucket_thick, bucket_ctr);


	/////////////////
	// Smarticle body
	/////////////////
	MySeed(964);
	double w_smarticle 	= sizeScale * 0.0117;
	double l_smarticle 	= 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	double t_smarticle 	= sizeScale * .00127;
	double t2_smarticle	= sizeScale * .0005;
	ChVector<> smarticleLengths(l_smarticle, w_smarticle, t_smarticle); // l, w, t
	ChVector<> sLenghWithTol = 1.3 * ChVector<>(smarticleLengths.x, smarticleLengths.y, 2 * smarticleLengths.z);

//	int nX = bucket_interior_halfDim.x / sLenghWithTol.x;
//	int nY = bucket_interior_halfDim.y / sLenghWithTol.z;
//	int nZ = 1;
	double maxDim = 1.3 * std::max(sLenghWithTol.x, sLenghWithTol.y);
	int nX = bucket_interior_halfDim.x / maxDim;
	int nY = bucket_interior_halfDim.y / maxDim;
	int nZ = 1;

	int smarticleCount = 0;
	for (int k = 0; k < nZ; k++) {
		for (int i= -nX + 1; i < nX; i++) {
			for (int j = -nY+1; j < nY; j ++) {
				ChQuaternion<> myRot = ChQuaternion<>(MyRand(), MyRand(), MyRand(), MyRand());
				myRot.Normalize();
				ChVector<> myPos = ChVector<>(i * maxDim, j * maxDim , k * maxDim);
//				ChVector<> myPos = ChVector<>(0, 0, bucket_interior_halfDim.z + (i%3) * sLenghWithTol.z)
//						+ ChVector<>(i * sLenghWithTol.x, j * sLenghWithTol.z , k * sLenghWithTol.y);

			  Smarticle * smarticle0 = new Smarticle(&mphysicalSystem, smarticleCount + 3 /* 1 and 2 are the first two objects */,
					  rho_smarticle, mat_g, l_smarticle, w_smarticle, t_smarticle, t2_smarticle,
					  myPos,
					  myRot);
			  smarticleCount++;
			//  smarticle0 = new Smarticle(&mphysicalSystem, 1, 1000, mat_g,
			//		  1, 1, .2, .05, S_BOX, ChVector<>(1,1,0), Q_from_AngAxis(CH_C_PI / 3, VECT_Y) * Q_from_AngAxis(CH_C_PI / 3, VECT_X));
			  smarticle0->Create();
			  mySmarticlesVec.push_back(smarticle0);
			}
		}
	}
	/////////////////
	// test body
	/////////////////

////  //////
//	double r = .1;
//	double l = 1;
//	double len = l;
//	double w = 3;
//  ChVector<> posRel = ChVector<>(-w/2 + r, l/2 - r, 1);
//  ChSharedPtr<ChBody> m_arm = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
//	m_arm->SetMaterialSurface(mat_g);
//
//	m_arm->SetPos(posRel);
////	m_arm->SetRot(QUNIT);
//  m_arm->SetCollide(true);
//  m_arm->SetBodyFixed(false);
//
//	double vol = utils::CalcCylinderVolume(r, len);
//	ChVector<> gyr = utils::CalcCylinderGyration(r, len).Get_Diag();
////	r = 1;
////	double vol = utils::CalcSphereVolume(r);
////	ChVector<> gyr = utils::CalcSphereGyration(r).Get_Diag();
//
//	double mass = 1000 * vol;
//
//	// create body
//    m_arm->SetMass(mass);
//    m_arm->SetInertiaXX(mass * gyr);
//
//  m_arm->GetCollisionModel()->ClearModel();
//	utils::AddCylinderGeometry(m_arm.get_ptr(), r, len, posRel, QUNIT);
////	utils::AddSphereGeometry(m_arm.get_ptr(), r);
//
//    m_arm->GetCollisionModel()->BuildModel();
//    mphysicalSystem.AddBody(m_arm);
//
//
//
//
//
//
//
//
//
//    ChSharedPtr<ChLinkLockRevolute> shipGroundPrismatic(new ChLinkLockRevolute);
//    shipGroundPrismatic->Initialize(
//    		smarticle0.GetArm(0), ground, true, ChCoordsys<>(ChVector<>(0, l/2, 0)), ChCoordsys<>(posRel + ChVector<>(0, l/2, 0)));
//    shipGroundPrismatic->SetName("ship_ground_prismatic");
//    mphysicalSystem.AddLink(shipGroundPrismatic);


//  ChSharedPtr<ChLinkLockRevolute> shipGroundPrismatic(new ChLinkLockRevolute);
//  shipGroundPrismatic->Initialize(
//  		smarticle0.GetArm(0), smarticle0.GetArm(1), true, ChCoordsys<>(ChVector<>(0, l/2, 0)), ChCoordsys<>(posRel + ChVector<>(0, l/2, 0)));
//  shipGroundPrismatic->SetName("ship_ground_prismatic");
//  mphysicalSystem.AddLink(shipGroundPrismatic);

}
// =============================================================================

void SavePovFilesMBD(ChSystemParallelDVI& mphysicalSystem,
                     int tStep) {
  int out_steps = std::ceil((1.0 / dT) / out_fps);
  printf("tStep %d , outstep %d \n", tStep, out_steps);

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
double Find_Max_Z(ChSystemParallelDVI& mphysicalSystem) {

	double zMax = -999999999;
	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
	for (int i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
		ChBody* bodyPtr = *(myIter + i);
		if ( strcmp(bodyPtr->GetName(), "smarticle_arm") == 0 ) {
			if (zMax < bodyPtr->GetPos().z) {
				zMax = bodyPtr->GetPos().z;
			}
		}
	}
	return zMax;
}
// =============================================================================
bool IsIn(ChVector<> pt, ChVector<> min, ChVector<> max) {
	if ((pt > max) || (pt < min)) {
		return false;
	}
	return true;
}
// =============================================================================
void PrintFractions(ChSystemParallelDVI& mphysicalSystem, int tStep) {
	const std::string vol_frac = out_dir + "/volumeFraction.txt";
	int stepSave = 10;
	if (tStep % stepSave != 0) return;

	std::ofstream vol_frac_of;
	if (tStep == 0) {
	  vol_frac_of.open(vol_frac);
	} else {
	  vol_frac_of.open(vol_frac, std::ios::app);
	}

	double zMax = Find_Max_Z(mphysicalSystem);
	zMax = std::min(zMax, bucket_ctr.z + bucket_interior_halfDim.z);

	int countInside = 0;
	double totalVolume = 0;
	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
	for (int i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
		ChBody* bodyPtr = *(myIter + i);
		if ( strcmp(bodyPtr->GetName(), "smarticle_arm") == 0 ) {
			if ( IsIn(bodyPtr->GetPos(), bucket_ctr - bucket_interior_halfDim, bucket_ctr + bucket_interior_halfDim) ) {
				countInside ++;
				totalVolume += bodyPtr->GetMass() / bodyPtr->GetDensity();
			}
		}
	}

	double volumeFraction = totalVolume / (8 * bucket_interior_halfDim.x * bucket_interior_halfDim.y * zMax);

	vol_frac_of << mphysicalSystem.GetChTime() << ", " << volumeFraction << std::endl;





	  vol_frac_of.close();
}
// =============================================================================

int main(int argc, char* argv[]) {
	  time_t rawtime;
	  struct tm* timeinfo;
	  time(&rawtime);
	  timeinfo = localtime(&rawtime);
	  const std::string simulationParams = out_dir + "/simulation_specific_parameters.txt";
	  simParams.open(simulationParams);
	  simParams << " Job was submitted at date/time: " << asctime(timeinfo) << std::endl;
	  ChTimerParallel step_timer;

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

  // Create a ChronoENGINE physical system
  ChSystemParallelDVI mphysicalSystem;
  InitializeMbdPhysicalSystem(mphysicalSystem, argc, argv);

  std::vector<Smarticle*> mySmarticlesVec;
  CreateMbdPhysicalSystemObjects(mphysicalSystem, mySmarticlesVec);



#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();

//	ChVector<> CameraLocation = ChVector<>(0, -10, 4);
//	ChVector<> CameraLookAt = ChVector<>(0, 0, -1);
	ChVector<> CameraLocation = sizeScale * ChVector<>(-.1, -.06, .06);
	ChVector<> CameraLookAt = sizeScale * ChVector<>(0, 0, -.01);
	gl_window.Initialize(1280, 720, "Smarticles", &mphysicalSystem);
	gl_window.SetCamera(CameraLocation, CameraLookAt, ChVector<>(0, 0, 1)); //camera
	gl_window.SetRenderMode(opengl::WIREFRAME);

// Uncomment the following two lines for the OpenGL manager to automatically
// run the simulation in an infinite loop.
// gl_window.StartDrawLoop(time_step);
// return 0;
#endif


  int stepEnd = int(tFinal / dT);  // 1.0e6;//2.4e6;//600000;//2.4e6 * (.02 * paramsH.sizeScale) /
  // ***************************** Simulation loop ********************************************

  ChSharedPtr<ChFunction> fun1 = ChSharedPtr<ChFunction>(new ChFunction_Const(0));
  ChSharedPtr<ChFunction> fun2 = ChSharedPtr<ChFunction>(new ChFunction_Ramp(0,1));
  ChSharedPtr<ChFunction> fun3 = ChSharedPtr<ChFunction>(new ChFunction_Ramp(0,-1));

  ChSharedPtr<ChFunction> fun4 = ChSharedPtr<ChFunction>(new ChFunction_Const(CH_C_PI / 2));
  ChSharedPtr<ChFunction> fun5 = ChSharedPtr<ChFunction>(new ChFunction_Const(-CH_C_PI / 2));

//  for (int i = 0; i < mySmarticlesVec.size(); i++) {
//	  mySmarticlesVec[i]->SetActuatorFunction(0, fun4);
//	  mySmarticlesVec[i]->SetActuatorFunction(1, fun4);
//
//  }

  double omega_bucket = 2 * CH_C_PI * 30;  // 30 Hz vibration similar to Gravish 2012, PRL
  double vibration_amp = sizeScale * 0.0002;

  for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
	  double t = mphysicalSystem.GetChTime();

	  // move bucket
	  double x_bucket = vibration_amp * sin(omega_bucket * t);
	  double xDot_bucket = omega_bucket * vibration_amp * cos(omega_bucket * t);
	  bucket->SetPos(ChVector<>(x_bucket, 0, 0));
	  bucket->SetPos_dt(ChVector<>(xDot_bucket, 0, 0));

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

#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(dT);
      gl_window.Render();
    }
#else
    mphysicalSystem.DoStepDynamics(dT);
#endif

	  step_timer.stop("step time");
	  printf("step: %d\n", tStep);

	  PrintFractions(mphysicalSystem, tStep);
  }
  for (int i = 0; i < mySmarticlesVec.size(); i++) {
	  delete mySmarticlesVec[i];

  }
  simParams.close();
  return 0;
}
