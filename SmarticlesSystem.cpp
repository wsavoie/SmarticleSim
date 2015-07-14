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
#include <omp.h>
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"     //Arman: why is this
#include "chrono_utils/ChUtilsInputOutput.h"  //Arman: Why is this
#include "chrono_utils/ChUtilsGenerators.h"

#include <ctime>
#include <stdlib.h>  // system, rand, srand, RAND_MAX
#include "core/ChFileutils.h" // for MakeDirectory
#include "Smarticle.h"
#include "SmarticleU.h"


//#include <fstream> // Arman: I don't know why this is not required
//#include <cstring> // Arman: I don't know why this is not required
//#include <stdio.h> // Arman: I don't know why this is not required



//#ifdef CHRONO_PARALLEL_HAS_OPENGL
//#undef CHRONO_PARALLEL_HAS_OPENGL
//#endif

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#if USE_PARALLEL
#define CH_SYSTEM ChSystemParallelDVI
#else
#define CH_SYSTEM ChSystem
#undef CHRONO_PARALLEL_HAS_OPENGL
#endif
//***********************************
// Use the namespace of Chrono

using namespace chrono;
//using namespace chrono::collision;
//using namespace std;

enum SmarticleType {SMART_ARMS , SMART_U};
// =============================================================================
std::ofstream simParams;
ChSharedPtr<ChBody> bucket;




	double sizeScale = 1;
	double gravity = -9.81 * sizeScale;
	double vibration_freq = 10;
	double dT = std::min(0.0001, 1.0 / vibration_freq / 200);;//std::min(0.0005, 1.0 / vibration_freq / 200);
	double contact_recovery_speed = 0.5 * sizeScale;
	double tFinal = 100;
	double rho_smarticle = 7850 / (sizeScale * sizeScale * sizeScale);
	int numLayers = 400;

	SmarticleType smarticleType = SMART_U;

	bool povray_output = true;
	int out_fps = 120;
	const std::string out_dir = "PostProcess";
	const std::string pov_dir_mbd = out_dir + "/povFilesSmarticles";

	ChVector<> bucket_ctr = ChVector<>(0,0,0);
	//ChVector<> Cbucket_interior_halfDim = sizeScale * ChVector<>(.05, .05, .025);
	ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(.05, .05, .025);
	//ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(.1, .1, .05);
	double bucket_thick = sizeScale * .005;

	// smarticle geometry
	double w_smarticle 	= sizeScale * 0.0117;
	double l_smarticle 	= 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	double t_smarticle 	= sizeScale * .00127;
	double t2_smarticle	= sizeScale * .0005;


// =============================================================================
void MySeed(double s = time(NULL)) { srand(s); }
double MyRand() { return float(rand()) / RAND_MAX; }
// =============================================================================
void SetArgumentsForMbdFromInput(int argc, char* argv[], int& threads, int& max_iteration_sliding, int& max_iteration_bilateral) {
  if (argc > 1) {
	const char* text = argv[1];
	double mult_l = atof(text);
	l_smarticle = mult_l * w_smarticle;
  }
  if (argc > 2) {
    const char* text = argv[2];
    threads = atoi(text);
  }
  if (argc > 3) {
    const char* text = argv[3];
    max_iteration_sliding = atoi(text);
  }
  if (argc > 4) {
    const char* text = argv[4];
    max_iteration_bilateral = atoi(text);
  }
}
// =============================================================================
void InitializeMbdPhysicalSystem_NonParallel(ChSystem& mphysicalSystem, int argc, char* argv[]) {
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
  int max_iteration_bilateral = 1000;

  // ----------------------
  // Set params from input
  // ----------------------

  SetArgumentsForMbdFromInput(argc, argv, threads, max_iteration_sliding, max_iteration_bilateral);

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
		  " l_smarticle mult for w (w = mult x l): " << l_smarticle /  w_smarticle << std::endl <<
		  " dT: " << dT << std::endl << std::endl;

  // ---------------------
  // Edit mphysicalSystem settings.
  // ---------------------

  // Modify some setting of the physical system for the simulation, if you want
  mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR); // LCP_ITERATIVE_SOR_MULTITHREAD , LCP_ITERATIVE_SOR
  mphysicalSystem.SetIterLCPmaxItersSpeed(50);
  mphysicalSystem.SetIterLCPmaxItersStab(5);   // unuseful for Anitescu, only Tasora uses this
  mphysicalSystem.SetParallelThreadNumber(1);
  mphysicalSystem.SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  mphysicalSystem.SetIterLCPwarmStarting(true);
  mphysicalSystem.SetUseSleeping(false);
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
		  " max_iteration_bilateral: " << max_iteration_bilateral << std::endl <<
		  " l_smarticle: " << l_smarticle << std::endl <<
		  " l_smarticle mult for w (w = mult x l): " << l_smarticle /  w_smarticle << std::endl <<
		  " dT: " << dT << std::endl << std::endl;

  // ---------------------
  // Edit mphysicalSystem settings.
  // ---------------------

  double tolerance = 0.001;  // 1e-3;  // Arman, move it to paramsH
  double collisionEnvelop = .4 * t2_smarticle;
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

  mphysicalSystem.GetSettings()->collision.collision_envelope = collisionEnvelop;
  mphysicalSystem.GetSettings()->collision.bins_per_axis = _make_int3(40, 40, 40);  // Arman check
}

// =============================================================================
void AddParticlesLayer(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec) {

	ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
	mat_g->SetFriction(0.5);
	/////////////////
	// Smarticle body
	/////////////////
	ChVector<> smarticleLengths(l_smarticle, w_smarticle, t_smarticle); // l, w, t
	ChVector<> sLenghWithTol = 1.3 * ChVector<>(smarticleLengths.x, smarticleLengths.y, 2 * smarticleLengths.z);

//	int nX = bucket_interior_halfDim.x / sLenghWithTol.x;
//	int nY = bucket_interior_halfDim.y / sLenghWithTol.z;
//	int nZ = 1;
	double maxDim = 1.3 * std::max(sLenghWithTol.x, sLenghWithTol.y);
	int nX = bucket_interior_halfDim.x / maxDim;
	int nY = bucket_interior_halfDim.y / maxDim;

	int smarticleCount = mySmarticlesVec.size();
	for (int i= -nX + 1; i < nX; i++) {
		for (int j = -nY+1; j < nY; j ++) {
			ChQuaternion<> myRot = ChQuaternion<>(MyRand(), MyRand(), MyRand(), MyRand());
			myRot.Normalize();
			ChVector<> myPos = ChVector<>(i * maxDim, j * maxDim , bucket_ctr.z + 5 * bucket_interior_halfDim.z); //1.5 * bucket_interior_halfDim.z to make sure it is above the pile
//				ChVector<> myPos = ChVector<>(0, 0, bucket_interior_halfDim.z + (i%3) * sLenghWithTol.z)
//						+ ChVector<>(i * sLenghWithTol.x, j * sLenghWithTol.z , k * sLenghWithTol.y);

			if (smarticleType == SMART_ARMS) {
				Smarticle * smarticle0  = new Smarticle(&mphysicalSystem);
				smarticle0->Properties(smarticleCount + 3 /* 1 and 2 are the first two objects, i.e. ground and bucket */,
								  rho_smarticle, mat_g, l_smarticle, w_smarticle, t_smarticle, t2_smarticle,
								  myPos,
								  myRot);
				smarticle0->Create();
				mySmarticlesVec.push_back((Smarticle*)smarticle0);
			} else if (smarticleType == SMART_U) {
				SmarticleU * smarticle0  = new SmarticleU(&mphysicalSystem);
				smarticle0->Properties(smarticleCount + 3 /* 1 and 2 are the first two objects */,
								  rho_smarticle, mat_g, l_smarticle, w_smarticle, t_smarticle, t2_smarticle,
								  myPos,
								  myRot);
				smarticle0->Create();
				mySmarticlesVec.push_back(smarticle0);
			} else {
				std::cout << "Error! Smarticle type is not set correctly" << std::endl;
			}
			smarticleCount++;
		}
	}
}

// =============================================================================
void CreateMbdPhysicalSystemObjects(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec) {
	  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
		mat_g->SetFriction(0.5);
//		mat_g->SetCohesion(0);
//		mat_g->SetCompliance(0.0);
//		mat_g->SetComplianceT(0.0);
//		mat_g->SetDampingF(0.2);

	/////////////////
	// Ground body
	/////////////////

	// ground
	ChVector<> boxDim = sizeScale * ChVector<>(0.1, 0.1, .002);
	ChVector<> boxLoc = sizeScale * ChVector<>(0, 0, -bucket_interior_halfDim.z - boxDim.z*5); // 1.1 to add 10% clearance between bucket and ground
	ChSharedPtr<ChBody> ground;
	if (USE_PARALLEL) {
		ground = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	} else {
		ground = ChSharedPtr<ChBody>(new ChBody);
	}
	ground->SetMaterialSurface(mat_g);
	ground->SetPos(boxLoc);

	// ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(true);

	ground->GetCollisionModel()->ClearModel();
	utils::AddCylinderGeometry(ground.get_ptr(), boxDim.x, boxDim.z, ChVector<>(0,0,0), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	ground->GetCollisionModel()->SetFamily(1);
	ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	ground->GetCollisionModel()->BuildModel();
	mphysicalSystem.AddBody(ground);

	// bucket
	if (USE_PARALLEL) {
		bucket = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	} else {
		bucket = ChSharedPtr<ChBody>(new ChBody);
	}

	// 1: create bucket
	bucket = utils::CreateBoxContainer(&mphysicalSystem, 1, mat_g, bucket_interior_halfDim, bucket_thick, bucket_ctr, QUNIT, true, false, true, false);
	bucket->GetCollisionModel()->SetFamily(1);
	bucket->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	bucket->SetBodyFixed(false);

	// 2: create plate
//	bucket->SetMaterialSurface(mat_g);
//	bucket->SetBodyFixed(true);
//	bucket->SetCollide(true);
//	bucket->GetCollisionModel()->ClearModel();
//	utils::AddBoxGeometry(bucket.get_ptr(), ChVector<>(bucket_interior_halfDim.x, bucket_interior_halfDim.y, bucket_thick), bucket_ctr - ChVector<>(0, 0, bucket_interior_halfDim.z));
//	bucket->GetCollisionModel()->BuildModel();
//	mphysicalSystem.AddBody(bucket);


	//	/////////////////
//	// Smarticle body
//	/////////////////

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
//	ChVector<> myPos = ChVector<>(nX / 2 * maxDim, nY / 2 * maxDim , 30 / 2 * maxDim + maxDim);
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






	/////////////////
	// test body
	/////////////////

////  //////
//	double r = .1;
//	double l = 1;
//	double len = l;
//	double w = 3;
//  ChVector<> posRel = ChVector<>(-w/2 + r, l/2 - r, 1);
	ChSharedPtr<ChBody> m_arm;
	if (USE_PARALLEL) {
		m_arm = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	} else {
		m_arm = ChSharedPtr<ChBody>(new ChBody);
	}
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
  printf("tStep %d , outstep %d, num bodies %d \n", tStep, out_steps, mphysicalSystem.Get_bodylist()->size());

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
		//smarticleTypeName = "smarticle_arm";
		smarticleTypeName = "smarticle_arm";
	} else if (smarticleType == SMART_U) {
		//smarticleTypeName = "smarticle_u";s
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
				zMax = bodyPtr->GetPos().z;
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
	ChVector<> bucketMin = bucket->GetPos() - bucket_interior_halfDim;

	ChVector<> mBuckCtr = bucket->GetPos();
//	std::cout << "bucket stuff " << mBuckCtr.x << ", " << mBuckCtr.y << ", " << mBuckCtr.z << std::endl;

	zMax = std::min(zMax, mBuckCtr.z + bucket_interior_halfDim.z);

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
//			if ( IsIn(bodyPtr->GetPos(), mBuckCtr - bucket_interior_halfDim, mBuckCtr + bucket_interior_halfDim) ) {
//				countInside ++;
//				totalVolume1 += bodyPtr->GetMass() / bodyPtr->GetDensity();
//			}
//		}
//	}

	double totalVolume2 = 0;
	int countInside2 = 0;
	for (int i = 0; i < mySmarticlesVec.size(); i ++) {
		Smarticle* sPtr = mySmarticlesVec[i];
		if ( IsIn(sPtr->Get_cm(), mBuckCtr - bucket_interior_halfDim, mBuckCtr + bucket_interior_halfDim) ) {
			countInside2 ++;
			totalVolume2 += sPtr->GetVolume();
		}
	}

	double volumeFraction = totalVolume2 / (4 * bucket_interior_halfDim.x * bucket_interior_halfDim.y * (zMax - (mBuckCtr.z - bucket_interior_halfDim.z)));

	vol_frac_of << mphysicalSystem.GetChTime() << ", " << countInside2  << ", " << volumeFraction << ", " << zMax << std::endl;

	vol_frac_of.close();
}
// =============================================================================

int main(int argc, char* argv[]) {
	  time_t rawtime;
	  struct tm* timeinfo;
	  time(&rawtime);
	  timeinfo = localtime(&rawtime);
	  ChTimerParallel step_timer;

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

  // Create a ChronoENGINE physical system
  CH_SYSTEM mphysicalSystem;
#if (USE_PARALLEL)
	  InitializeMbdPhysicalSystem_Parallel(mphysicalSystem, argc, argv);
#else
	  InitializeMbdPhysicalSystem_NonParallel(mphysicalSystem, argc, argv);
#endif

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
	gl_window.viewer->render_camera.camera_scale = 2.0/(1000.0)*sizeScale;
	gl_window.viewer->render_camera.near_clip = .001;
	gl_window.SetRenderMode(opengl::WIREFRAME);

// Uncomment the following two lines for the OpenGL manager to automatically
// run the simulation in an infinite loop.
// gl_window.StartDrawLoop(time_step);
// return 0;
#endif


  int stepEnd = int(tFinal / dT);  // 1.0e6;//2.4e6;//600000;//2.4e6 * (.02 * paramsH.sizeScale) /
  // ***************************** Simulation loop ********************************************

  ChSharedPtr<ChFunction_Const> fun1 = ChSharedPtr<ChFunction_Const>(new ChFunction_Const(0));
  ChSharedPtr<ChFunction_Ramp> fun2 = ChSharedPtr<ChFunction_Ramp>(new ChFunction_Ramp(0,1));
  ChSharedPtr<ChFunction_Ramp> fun3 = ChSharedPtr<ChFunction_Ramp>(new ChFunction_Ramp(0,-1));

  ChSharedPtr<ChFunction_Const> fun4 = ChSharedPtr<ChFunction_Const>(new ChFunction_Const(CH_C_PI / 2));
  ChSharedPtr<ChFunction_Const> fun5 = ChSharedPtr<ChFunction_Const>(new ChFunction_Const(-CH_C_PI / 2));

//  for (int i = 0; i < mySmarticlesVec.size(); i++) {
//	  mySmarticlesVec[i]->SetActuatorFunction(0, fun4);
//	  mySmarticlesVec[i]->SetActuatorFunction(1, fun4);
//
//  }

  double omega_bucket = 2 * CH_C_PI * vibration_freq;  // 30 Hz vibration similar to Gravish 2012, PRL
  double vibration_amp = sizeScale * 0.001;

  double timeForVerticalDisplcement = 1.5 * sqrt(2 * w_smarticle / mphysicalSystem.Get_G_acc().Length()); // 1.5 for safety proximity
  int numGeneratedLayers = 0;

//  for (int tStep = 0; tStep < 1; tStep++) {
  for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
	  double t = mphysicalSystem.GetChTime();

	  if (  (fmod(mphysicalSystem.GetChTime(), timeForVerticalDisplcement) < dT)  &&
			  (numGeneratedLayers < numLayers) ){
		  AddParticlesLayer(mphysicalSystem, mySmarticlesVec);
		  numGeneratedLayers ++;
	  }

	   printf("\n");
	  // move bucket
	  double x_bucket = vibration_amp * sin(omega_bucket * t);
	  double xDot_bucket = omega_bucket * vibration_amp * cos(omega_bucket * t);
	  bucket->SetPos(ChVector<>(x_bucket, 0, 0));
	  bucket->SetPos_dt(ChVector<>(xDot_bucket, 0, 0));
	  bucket->SetRot(QUNIT);
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

	  time(&rawtimeCurrent);
	  double timeDiff = difftime(rawtimeCurrent, rawtime);

	  step_timer.stop("step time");
	  std::cout << "step time: " << step_timer.GetTime("step time") << ", time passed: " << int(timeDiff)/3600 <<":"<< (int(timeDiff) % 3600) / 60 << ":" << (int(timeDiff) % 60) <<std::endl;

	  std::cout << "particle velocity: " << ((SmarticleU*)mySmarticlesVec[0])->GetSmarticleBodyPointer()->GetPos_dt().Length() << std::endl;
	  FixBodies(mphysicalSystem, tStep);
	  PrintFractions(mphysicalSystem, tStep, mySmarticlesVec);

	  std::cout.flush();

  }
  for (int i = 0; i < mySmarticlesVec.size(); i++) {
	  delete mySmarticlesVec[i];

  }
  mySmarticlesVec.clear();
  simParams.close();
  return 0;
}
