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
#include <stdlib.h>  // system
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

std::ofstream simParams;

double gravity = -9.81;
double dT = .001;
double contact_recovery_speed = 0.1;
double tFinal = 2000;

bool povray_output = true;
int out_fps = 120;
const std::string out_dir = "PostProcess";
const std::string pov_dir_mbd = out_dir + "/povFilesSmarticles";

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
  bool thread_tuning = false;

  //	uint max_iteration = 20;//10000;
  int max_iteration_normal = 0;
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

void CreateMbdPhysicalSystemObjects(ChSystemParallelDVI& mphysicalSystem) {
	  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
		mat_g->SetFriction(0.1);
		mat_g->SetCohesion(0);
		mat_g->SetCompliance(0.0);
		mat_g->SetComplianceT(0.0);
		mat_g->SetDampingF(0.2);

	/////////////////
	// Ground body
	/////////////////

	ChVector<> boxDim(10, 10, .2);
	ChVector<> boxLoc(0, 0, -3);

  ChSharedPtr<ChBody> ground = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
//  ground->SetIdentifier(-1);
  ground->SetBodyFixed(true);
  ground->SetCollide(true);
  ground->SetMaterialSurface(mat_g);
  ground->GetCollisionModel()->ClearModel();

  utils::AddBoxGeometry(ground.get_ptr(), boxDim, boxLoc, ChQuaternion<>(1, 0, 0, 0));

  ground->GetCollisionModel()->BuildModel();
  mphysicalSystem.AddBody(ground);

	/////////////////
	// Smarticle body
	/////////////////

//  Smarticle smarticle0(&mphysicalSystem, 1, 1000, mat_g,
//		  1, 1, .1, .05, S_CYLINDER);
//  smarticle0.Create();

  //////
	double r = .1;
	double l = 1;
	double len = l;
	double w = 1;
  ChVector<> posRel = ChVector<>(-w/2 + r, l/2 - r, 0);
  ChSharedBodyPtr arm;
  arm = ChSharedBodyPtr(new ChBody(new collision::ChCollisionModelParallel));
	arm->SetPos(posRel);
	arm->SetRot(QUNIT);
  arm->SetCollide(true);
  arm->SetBodyFixed(false);
	arm->SetMaterialSurface(mat_g);

	arm->GetCollisionModel()->ClearModel();

	double vol = utils::CalcCylinderVolume(r, len);
	ChVector<> gyr = utils::CalcCylinderGyration(r, len).Get_Diag();
	utils::AddCylinderGeometry(arm.get_ptr(), r, len, posRel, QUNIT);

	double mass = 1000 * vol;

	// create body
    arm->SetMass(mass);
    arm->SetInertiaXX(mass * gyr);
    arm->GetCollisionModel()->BuildModel();
    mphysicalSystem.AddBody(arm);



}
// =============================================================================

void SavePovFilesMBD(ChSystemParallelDVI& mphysicalSystem,
                     int tStep) {
  int out_steps = std::ceil((1.0 / dT) / out_fps);

  static int out_frame = 0;

  // If enabled, output data for PovRay postprocessing.
  if (povray_output && tStep % out_steps == 0) {
		if (tStep / out_steps == 0) {
			const std::string rmCmd = std::string("rm ") + pov_dir_mbd + std::string("/*.dat");
			system(rmCmd.c_str());
		}

    char filename[100];
    sprintf(filename, "%s/data_%03d.dat", pov_dir_mbd.c_str(), out_frame + 1);
    utils::WriteShapesPovray(&mphysicalSystem, filename);

    out_frame++;
  }
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

  // Create a ChronoENGINE physical system
  ChSystemParallelDVI mphysicalSystem;
  InitializeMbdPhysicalSystem(mphysicalSystem, argc, argv);
  CreateMbdPhysicalSystemObjects(mphysicalSystem);


#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();

//	ChVector<> CameraLocation = ChVector<>(0, -10, 4);
//	ChVector<> CameraLookAt = ChVector<>(0, 0, -1);
	ChVector<> CameraLocation = ChVector<>(0, -5, 20);
	ChVector<> CameraLookAt = ChVector<>(0, 0, -1);
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

  for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
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


  }
  simParams.close();
  return 0;
}
