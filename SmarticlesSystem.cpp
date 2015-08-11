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
//   Multibody dynamics engine
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
#include "CheckPointSmarticles.h"
//#include "D:/ChronoCode/libs/png++-0.2.7/png.hpp";
//#include "D:\ChronoCode\libs\pngwriter-release-0.5.5\src\pngwriter.h"

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif
#define GLEW_STATIC
#include <GL/glew.h>
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")

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
using namespace chrono;
enum SmarticleType {SMART_ARMS , SMART_U};
enum BucketType { CYLINDER, BOX };
SmarticleType smarticleType = SMART_U;
BucketType bucketType = CYLINDER;
// =============================================================================
double Find_Max_Z(CH_SYSTEM& mphysicalSystem);
std::ofstream simParams;
ChSharedPtr<ChBody> bucket;

	double sizeScale = 1;
	double gravity = -9.81 * sizeScale;
	
	double vibration_freq = 30;
	bool read_from_file = true;
	double omega_bucket = 2 * CH_C_PI * vibration_freq;  // 30 Hz vibration similar to Gravish 2012, PRL
	//double vibration_amp = sizeScale * 0.00055;
	double mGamma = 1.23 * gravity;
	double vibration_amp = mGamma / (omega_bucket*omega_bucket);
	


	//double dT = std::min(0.001, 1.0 / vibration_freq / 200);;//std::min(0.0005, 1.0 / vibration_freq / 200);
	double dT = 0.001;//std::min(0.0005, 1.0 / vibration_freq / 200);
	double contact_recovery_speed = 0.2 * sizeScale;
	double tFinal = 10.0;
	double vibrateStart= tFinal-5.0;
	double collapseStart = .2;
	double removeWallStart = .05;

	double rho_smarticle = 7850.0 / (sizeScale * sizeScale * sizeScale);
	double rho_cylinder = 1180.0 / (sizeScale * sizeScale * sizeScale);
	ChSharedPtr<ChMaterialSurface> mat_g;
	int numLayers = 100;
	double armAngle1 = 90;
	double armAngle2 = 90;
	

	bool povray_output = true;
	int out_fps = 120;
	const std::string out_dir = "PostProcess";
	const std::string pov_dir_mbd = out_dir + "/povFilesSmarticles";

	ChVector<> bucket_ctr = ChVector<>(0,0,0);
	//ChVector<> Cbucket_interior_halfDim = sizeScale * ChVector<>(.05, .05, .025);
	double bucket_rad = sizeScale*0.022;
	ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, .010);

	
	//ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(.1, .1, .05);
	double bucket_half_thick = sizeScale * .005;
	double h = bucket_interior_halfDim.z * 2 + 2*bucket_half_thick; //from entangled paper height of available volume for smarticles in bucket
	double d = bucket_interior_halfDim.y * 2; //from entangled paper width of availble volume for smarticles in bucket

	// smarticle geometry//
	double w_smarticle 	= sizeScale * 0.0117;
	double l_smarticle 	= 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
	double t_smarticle 	= sizeScale * .00127;
	double t2_smarticle	= sizeScale * .0005;
	double m_smarticle = (t_smarticle)* (t2_smarticle)* (w_smarticle + 2 * l_smarticle);
	double collisionEnvelope = .4 * t2_smarticle;

// =============================================================================
void MySeed(double s = time(NULL)) { srand(s); }
double MyRand() { return float(rand()) / RAND_MAX; }
// =============================================================================
void SetArgumentsForMbdFromInput(int argc, char* argv[], int& threads, int& max_iteration_sliding, int& max_iteration_bilateral, double& dt, int& num_layers, double& mangle1,double& mangle2,double& gamma, bool& readFile) {
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
		mangle1 = atof(text);
	}
	if (argc > 5){
		const char* text = argv[5];
		mangle2 = atof(text);
	}
	if (argc > 6){
		const char* text = argv[6];
		gamma = atof(text)*gravity;
		vibration_amp= gamma / (omega_bucket*omega_bucket);
	}
	if (argc > 7){
		const char* text = argv[7];
		readFile = atoi(text);
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


  // ---------------------
  // Print the rest of parameters
  // ---------------------

	int dummyNumber0;
	int dummyNumber1;
	int dummyNumber2;
  SetArgumentsForMbdFromInput(argc, argv, dummyNumber0, dummyNumber1, dummyNumber2, dT,numLayers, armAngle1,armAngle2,mGamma,read_from_file);
	simParams << std::endl <<
		" l_smarticle: " << l_smarticle << std::endl <<
		" l_smarticle mult for w (w = mult x l): " << l_smarticle / w_smarticle << std::endl <<
		" dT: " << dT << std::endl <<
		" tFinal: " << tFinal << std::endl <<
		" vibrate start: " << vibrateStart << std::endl <<
		" Gamma: " << mGamma << std::endl <<
		" numlayers: " << numLayers << std::endl <<
		" read from file: " << read_from_file << std::endl <<
		" arm angle 1 and 2: " << armAngle1 << " " << armAngle2 << std::endl <<
		"	non-parallel" << std::endl;

	if (argc > 7)
	{
		simParams << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6] << " " << argv[7] << std::endl;
	}

  // ---------------------
  // Edit mphysicalSystem settings.
  // ---------------------

  // Modify some setting of the physical system for the simulation, if you want
  mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR); // LCP_ITERATIVE_SOR_MULTITHREAD , LCP_ITERATIVE_SOR  (LCP_ITERATIVE_SOR_MULTITHREAD does not work)
  mphysicalSystem.SetIterLCPmaxItersSpeed(50);
  mphysicalSystem.SetIterLCPmaxItersStab(5);   // unuseful for Anitescu, only Tasora uses this
  mphysicalSystem.SetParallelThreadNumber(1);
  mphysicalSystem.SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
  mphysicalSystem.SetIterLCPwarmStarting(true);
  mphysicalSystem.SetUseSleeping(false);
  mphysicalSystem.Set_G_acc(ChVector<>(0, 0, gravity));
}
// =============================================================================
void InitializeMbdPhysicalSystem_Parallel(ChSystemParallelDVI& mphysicalSystem, int argc, char* argv[]) {
	// initializd random seeder
	MySeed();
  // Desired number of OpenMP threads (will be clamped to maximum available)
  int threads = 1;
  // Perform dynamic tuning of number of threads?
  bool thread_tuning = true;

  //	uint max_iteration = 20;//10000;
  int max_iteration_normal = 250;
  int max_iteration_sliding = 250;
  int max_iteration_spinning = 0;
  int max_iteration_bilateral = 50;

  // ----------------------
  // Set params from input
  // ----------------------
  SetArgumentsForMbdFromInput(argc, argv, threads, max_iteration_sliding, max_iteration_bilateral, dT,numLayers, armAngle1,armAngle2,mGamma,read_from_file);

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
	//mphysicalSystem.GetSettings()->max_threads = std::min(max_threads, int(3.0 * threads / 2));
	mphysicalSystem.GetSettings()->max_threads =max_threads;
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
		" Gamma: " << mGamma << std::endl <<
		" numlayers: " << numLayers << std::endl <<
		" read from file: " << read_from_file << std::endl << 
		" arm angle 1 and 2: " << armAngle1<< " " << armAngle2 << std::endl << 
		" parallel" << std::endl << std::endl;

	if (argc > 10)
	{
		simParams << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6] << " " << argv[7] << " " << argv[8] << " " << argv[9] << " " << argv[10] << std::endl;
	}

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

// =============================================================================
void AddParticlesLayer(CH_SYSTEM& mphysicalSystem, std::vector<Smarticle*> & mySmarticlesVec) {
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
	int numPerLayer = 9;
	//this filling method works better for smaller diameter where diameter < 3*width of staple	
	if (w_smarticle * 6 > bucket_interior_halfDim.x * 2)
	{
		for (int i = 0; i < numPerLayer; i++)
		{
			ChQuaternion<> myRot = ChQuaternion<>(MyRand(), MyRand(), MyRand(), MyRand());
			myRot.Normalize();

			ChVector<> myPos = ChVector<>(bucket_ctr.x + MyRand()*bucket_interior_halfDim.x - MyRand()*bucket_interior_halfDim.x / 2.0,
				bucket_ctr.y + MyRand()*bucket_interior_halfDim.y - MyRand()*bucket_interior_halfDim.y / 2.0,
				std::min(9*bucket_interior_halfDim.z ,z)+i*w_smarticle/4);

			if (smarticleType == SMART_ARMS) {
				Smarticle * smarticle0 = new Smarticle(&mphysicalSystem);
				smarticle0->Properties(smarticleCount,
					rho_smarticle, mat_g,
					collisionEnvelope,
					l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
					myPos,
					myRot);
				smarticle0->Create();
				mySmarticlesVec.push_back((Smarticle*)smarticle0);
			}
			else if (smarticleType == SMART_U) {
				SmarticleU * smarticle0 = new SmarticleU(&mphysicalSystem);
				smarticle0->Properties(smarticleCount,
					rho_smarticle, mat_g,
					collisionEnvelope,
					l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
					myPos,
					myRot);
					//smarticle0->SetAngle(armAngle, true);
				smarticle0->SetAngle(armAngle1,armAngle2, true);
				smarticle0->Create();
				mySmarticlesVec.push_back(smarticle0);
			}
			else {
				std::cout << "Error! Smarticle type is not set correctly" << std::endl;
			}
			smarticleCount++;
		}
	}
	else
	{
		for (int i = -nX + 1; i < nX; i++) {
			for (int j = -nY + 1; j < nY; j++) {
				ChQuaternion<> myRot = ChQuaternion<>(MyRand(), MyRand(), MyRand(), MyRand());
				myRot.Normalize();


				ChVector<> myPos = ChVector<>(i * maxDim + bucket_ctr.x + 2 * MyRand()*w_smarticle - w_smarticle
					, j * maxDim + bucket_ctr.y + 2 * MyRand()*w_smarticle - w_smarticle
					, z + 0.5 * maxDim);

				//ChVector<> myPos = ChVector<>(i * maxDim + bucket_ctr.x, j * maxDim + bucket_ctr.y, bucket_ctr.z + 6.0 * bucket_interior_halfDim.z + 2 * bucket_half_thick);


				//ChVector<> myPos = ChVector<>(i * maxDim, j * maxDim, bucket_ctr.z + 6.0 * bucket_interior_halfDim.z + 2 * bucket_half_thick);
				// ***  added 2*bucket_half_thick to make sure stuff are initialized above bucket. Remember, bucket is inclusive, i.e. the sizes are extende 2*t from each side

				if (smarticleType == SMART_ARMS) {
					Smarticle * smarticle0 = new Smarticle(&mphysicalSystem);
					smarticle0->Properties(smarticleCount,
						rho_smarticle, mat_g,
						collisionEnvelope,
						l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
						myPos,
						myRot);
					smarticle0->Create();
					mySmarticlesVec.push_back((Smarticle*)smarticle0);
				}
				else if (smarticleType == SMART_U) {
					SmarticleU * smarticle0 = new SmarticleU(&mphysicalSystem);
					smarticle0->Properties(smarticleCount,
						rho_smarticle, mat_g,
						collisionEnvelope,
						l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
						myPos,
						myRot);
					//smarticle0->SetAngle(armAngle, true);
					smarticle0->SetAngle(armAngle1,armAngle2, true);
					smarticle0->Create();
					mySmarticlesVec.push_back(smarticle0);
				}
				else {
					std::cout << "Error! Smarticle type is not set correctly" << std::endl;
				}
				smarticleCount++;
			}
		}
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
	double half_height = bucket_interior_halfDim.z;
	double box_side = bucket_rad * 2.0 * tan(CH_C_PI / num_boxes);//side length of cyl
	double o_lap = 0;
	if (overlap){ o_lap = t * 2; }
	double ang = 2.0 * CH_C_PI / num_boxes;
	ChVector<> box_size = (0, 0, 0); //size of plates
	ChVector<> pPos = (0, 0, 0);  //position of each plate
	ChQuaternion<> quat = QUNIT; //rotation of each plate
	cyl_container->GetCollisionModel()->ClearModel();
	cyl_container->SetMaterialSurface(wallMat);
	for (int i = 0; i < num_boxes; i++)
	{

		box_size = ChVector<>((box_side + t) / 2.0,
			t,
			2 * half_height + o_lap);

		pPos = bucket_ctr + ChVector<>(sin(ang * i) * (t + bucket_rad),
			cos(ang*i)*(t + bucket_rad),
			2 * half_height);

		quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));

		//this is here to make half the cylinder invisible.
		bool m_visualization = false;
		if (ang*i < CH_C_PI || ang*i > 3.0 * CH_C_PI / 2.0)
		{
			m_visualization = true;
		}
		utils::AddBoxGeometry(cyl_container.get_ptr(), box_size, pPos, quat, m_visualization);
	}

	//Add ground piece
	//
	if (!read_from_file)	{
		//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad + 2 * t, t, ChVector<>(0, 0, -t), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
		utils::AddBoxGeometry(cyl_container.get_ptr(), Vector(bucket_rad+t, bucket_rad + t, t), Vector(0, 0, -t), QUNIT, true);
	}
	else{
		cyl_container->SetPos(cyl_container->GetPos() + Vector(0, 0, vibration_amp*sin(CH_C_PI / 2.0)));//to place box in way which is maximum downward so upon creation it has no chance of starting particles inside
	}

	//add up volume of bucket and multiply by rho to get mass;
	double cyl_volume = CH_C_PI*(2 * box_size.z - 2 * t)*(2 * box_size.z - 2 * t)*((2 * bucket_rad + 2 * t)*(2 * bucket_rad + 2 * t) - bucket_rad*bucket_rad) + (CH_C_PI)*(bucket_rad + 2 * t)*(bucket_rad + 2 * t) * 2 * t;
	cyl_container->SetMass(rho_cylinder*cyl_volume);
	//utils::AddBoxGeometry(cyl_container.get_ptr(), Vector(bucket_rad, bucket_rad + t, t), Vector(0, 0, -t), QUNIT, true);

	//checks top,bottom, and middle location
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos() + Vector(0,0,2 * bucket_interior_halfDim.z + 2 * bucket_half_thick), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos(), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos() + Vector(0, 0, bucket_interior_halfDim.z), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	
	//ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, bucket_interior_halfDim.z);
	
	cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
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
	ground->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	ground->GetCollisionModel()->BuildModel();
	mphysicalSystem.AddBody(ground);
	// bucket
	if (USE_PARALLEL) {
		bucket = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	} else {
		bucket = ChSharedPtr<ChBody>(new ChBody);
	}


	// 1: create bucket
		mat_g->SetFriction(0.4); //steel- plexiglass   (plexiglass was outer cylinder material)
	if (bucketType == BOX){
		bucket = utils::CreateBoxContainer(&mphysicalSystem, 1, mat_g, bucket_interior_halfDim, bucket_half_thick, bucket_ctr, QUNIT, true, false, true, false);
	}
	if (bucketType == CYLINDER){
		//http://www.engineeringtoolbox.com/friction-coefficients-d_778.html to get coefficients
		bucket = create_cylinder_from_blocks(25, 1, true, &mphysicalSystem, mat_g);
//		utils::CreateCylindricalContainerFromBoxes(&mphysicalSystem, 1, mat_g, ChVector<>(bucket_interior_halfDim.x, bucket_interior_halfDim.y, 2 * bucket_interior_halfDim.z), bucket_half_thick, 25, rho_cylinder, collisionEnvelope, bucket_ctr, QUNIT, true, true, false, true);
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
	//  bucket->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
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
	//					  rho_smarticle, mat_g,
	//  					collisionEnvelope,
	//  					l_smarticle, w_smarticle, t_smarticle, t2_smarticle,

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
				zMax = bodyPtr->GetPos().z + bucket->GetPos().z;
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
//rad is (radius from center, , max z)
bool IsInRadial(ChVector<> pt, ChVector<> centralPt, ChVector<> rad)
{
	ChVector<> dist = pt - centralPt;
	double xydist = (std::sqrt(dist.x * dist.x + dist.y + dist.y));
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
		if (bodyPtr->GetPos().z < -4.0 * bucket_interior_halfDim.z) {
			bodyPtr->SetBodyFixed(true);
			continue;
		}
		if (bodyPtr->GetPos().x > 5 * bucket_rad || bodyPtr->GetPos().x < -5 * bucket_rad) {
			bodyPtr->SetBodyFixed(true);
			continue;
		}

		if (bodyPtr->GetPos().y > 5 * bucket_rad || bodyPtr->GetPos().y < -5 * bucket_rad) {
			bodyPtr->SetBodyFixed(true);
			continue;
		}
		if (bodyPtr->GetRot_dt().GetVector().Length2() > 10000)
		{
			bodyPtr->SetRot_dt(QUNIT);
		}
	}
}
void FixBodies(CH_SYSTEM& mphysicalSystem, int tStep, std::vector<Smarticle*> mySmarticlesVec) {
	for (int i = 0; i < mySmarticlesVec.size(); i++) {
		Smarticle* sPtr = mySmarticlesVec[i];
		if (sPtr->Get_cm().z < -4.0 * bucket_interior_halfDim.z) {
			sPtr->SetBodyFixed(true); //Could/Should we write the destructor to remove these from system and then remove them from the vector too?
				continue;
			}
		if (sPtr->Get_cm().x > 5 * bucket_rad || sPtr->Get_cm().x < -5 * bucket_rad) {
			sPtr->SetBodyFixed(true);
			continue;
		}

		if (sPtr->Get_cm().y > 5 * bucket_rad || sPtr->Get_cm().y < -5 * bucket_rad) {
			sPtr->SetBodyFixed(true);
			continue;
		}
	}
}

// =============================================================================
void PrintFractionsAndCOM(CH_SYSTEM& mphysicalSystem, int tStep, std::vector<Smarticle*> mySmarticlesVec,ChVector<> bucketMin) 
{
	double zCom = 0;
	const std::string vol_frac = out_dir + "/volumeFraction.txt";
	int stepSave = 10;
	if (tStep % stepSave != 0) return;

	std::ofstream vol_frac_of;
	if (tStep == 0) {
		vol_frac_of.open(vol_frac.c_str());
	}
	else {
		vol_frac_of.open(vol_frac.c_str(), std::ios::app);
	}

	double zMax = Find_Max_Z(mphysicalSystem);
//	ChVector<> bucketMin = bucket->GetPos();

	zMax = std::min(zMax, bucketMin.z + 2 * bucket_interior_halfDim.z);

	ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, bucket_interior_halfDim.z);

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
				zCom += sPtr->Get_cm().z*m_smarticle;

			}
		}
		
		volumeFraction = totalVolume2 / (4 * bucket_interior_halfDim.x * bucket_interior_halfDim.y * (zMax - bucketMin.z));
	}
	if (bucketType == CYLINDER)
	{
		for (int i = 0; i < mySmarticlesVec.size(); i++) {
			Smarticle* sPtr = mySmarticlesVec[i];
			
			double rad = bucket_rad; 
			
			if (read_from_file){
				rad = 6 * bucket_rad+2*bucket_half_thick;//allows for larger region to test for com in when column collapses
			}

			//isinradial rad parameter is Vector(bucketrad,zmin,zmax)
			if (IsInRadial(sPtr->Get_cm(), bucketCtr, ChVector<>(rad, bucketMin.z, bucketMin.z + 2 * bucket_interior_halfDim.z))) {
				countInside2++;
				totalVolume2 += sPtr->GetVolume();
				zCom += sPtr->Get_cm().z*m_smarticle;
			}
		}

		//volumeFraction = totalVolume2 / (CH_C_PI * bucket_rad * bucket_rad * 2 * bucket_interior_halfDim.z);
		volumeFraction = totalVolume2 / (CH_C_PI * bucket_rad * bucket_rad * zMax);
		
	}
	zCom = zCom / (countInside2*m_smarticle);
	vol_frac_of << mphysicalSystem.GetChTime() << ", " << countInside2 << ", " << volumeFraction << ", " << zMax <<", "<< zCom<< std::endl;

	vol_frac_of.close();
}
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

	//zMax = std::min(zMax, bucketMin.z + 2 * bucket_interior_halfDim.z + 2*bucket_half_thick);
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
			if (IsInRadial(sPtr->Get_cm(), bucketCtr, ChVector<>(bucket_rad, bucketMin.z, bucketMin.z + 2 * bucket_interior_halfDim.z))) {
			//if (IsInRadial(sPtr->Get_cm(), bucketCtr, ChVector<>(bucket_rad, bucketMin.z, bucketMin.z+2*bucket_interior_halfDim.z+2*bucket_half_thick))) {
				countInside2++;
				totalVolume2 += sPtr->GetVolume();
			}
		}

		//volumeFraction = totalVolume2 / (CH_C_PI * bucket_rad * bucket_rad * 2 * bucket_interior_halfDim.z);
		volumeFraction = totalVolume2 / (CH_C_PI * bucket_rad * bucket_rad * zMax);
	}

	vol_frac_of << mphysicalSystem.GetChTime() << ", " << countInside2  << ", " << volumeFraction << ", " << zMax << std::endl;

	vol_frac_of.close();
}
// =============================================================================
// move bucket
void vibrate_bucket(double t,ChSharedPtr<chrono::ChBody> body) {
	//double x_bucket = vibration_amp*sin(omega_bucket * t);
	//this allows for that at t=vibration start time, the bucket wont jump some amount
	double phase = -omega_bucket*vibrateStart;
	double x_bucket = vibration_amp*sin(omega_bucket * t+phase);
	double xDot_bucket = vibration_amp*omega_bucket*cos(omega_bucket * t+phase);
	double xDDot_bucket = vibration_amp*omega_bucket*omega_bucket*-1 * sin(omega_bucket * t+phase);
	body->SetPos(ChVector<>(0, 0, x_bucket));
	body->SetPos_dt(ChVector<>(0, 0, xDot_bucket));
	body->SetPos_dtdt(ChVector<>(0, 0, xDDot_bucket));
	body->SetRot(QUNIT);
}
//bool screenshot2(char *fileName){
//png::image< png::rgb_pixel > image(128, 128);
//for (size_t y = 0; y < image.get_height(); ++y)
//{
//	for (size_t x = 0; x < image.get_width(); ++x)
//	{
//		image[y][x] = png::rgb_pixel(x, y, x + y);
//		// non-checking equivalent of image.set_pixel(x, y, ...);
//	}
//}
//image.write("rgb.png");
//}
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
		
  std::vector<Smarticle*> mySmarticlesVec;
  CreateMbdPhysicalSystemObjects(mphysicalSystem, mySmarticlesVec);
	
#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	
//	ChVector<> CameraLocation = ChVector<>(0, -10, 4);
//	ChVector<> CameraLookAt = ChVector<>(0, 0, -1);
	ChVector<> CameraLocation = sizeScale * ChVector<>(-.1, -.06, .1);
	ChVector<> CameraLookAt = sizeScale * ChVector<>(0, 0, -.01);
	char appTitle[240];
	sprintf(appTitle,"Smarticle: lw: %g, angs: %g,%g, numlayers: %d, dT: %g, gamma: %g", l_smarticle / w_smarticle, armAngle1,armAngle2, numLayers, dT,mGamma);
	gl_window.Initialize(1280, 720, appTitle, &mphysicalSystem);
	gl_window.viewer->contact_renderer.SetPointSize(.001);
	gl_window.viewer->cloud.SetPointSize(0.001);

	
	gl_window.SetCamera(CameraLocation, CameraLookAt, ChVector<>(0, 0, 1)); //camera
	gl_window.viewer->render_camera.camera_scale = 2.0/(1000.0)*sizeScale;
	gl_window.viewer->render_camera.near_clip = .001;
	gl_window.SetRenderMode(opengl::SOLID);


// Uncomment the following two lines for the OpenGL manager to automatically
// run the simulation in an infinite loop.
// gl_window.StartDrawLoop(time_step);
// return 0;
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

	double timeForVerticalDisplcement = 0.05; // 1.5 for safety proximity
 
	int numGeneratedLayers = 0;
	
	

//  CheckPointSmarticles_Read(mphysicalSystem, mySmarticlesVec);

  printf("************** size sys %d \n", mySmarticlesVec.size());

	// bucket
	ChSharedPtr<ChBody> bucket_bott;
	if (USE_PARALLEL) {
		bucket_bott = ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));
	}
	else {
		bucket_bott = ChSharedPtr<ChBody>(new ChBody);
	}
	

	if (read_from_file)
	{ 
		CheckPointSmarticles_Read(mphysicalSystem, mySmarticlesVec);

		bucket_bott->GetCollisionModel()->ClearModel();

		utils::AddBoxGeometry(bucket_bott.get_ptr(), Vector(bucket_rad+2*bucket_half_thick, bucket_rad + 2 * bucket_half_thick, bucket_half_thick), Vector(0, 0, -bucket_half_thick), QUNIT, true);
		//utils::AddCylinderGeometry(bucket_bott.get_ptr(), bucket_rad*6 +2*bucket_half_thick, bucket_half_thick, bucket->GetPos() + ChVector<>(0, 0, -bucket_half_thick), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
		bucket_bott->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
		bucket_bott->GetCollisionModel()->BuildModel();

		bucket_bott->SetCollide(true);
		bucket_bott->SetBodyFixed(true);
		bucket_bott->GetCollisionModel()->SetFamily(1);
		bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

		mphysicalSystem.AddBody(bucket_bott);
	}
	bool bucket_exist = true;
//  for (int tStep = 0; tStep < 1; tStep++) {
	
  for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
		double t = mphysicalSystem.GetChTime();

		int sSize1 = mySmarticlesVec.size();
		if (read_from_file)
		{	
			if (t >= removeWallStart){
			
					//bucket_bott->SetBodyFixed(false);
					//vibrate_bucket(t,bucket_bott);
					if (bucket_exist){
						bucket->SetPos(ChVector<>(100, 0, 0));
						bucket_exist = false;
					}
					if (t >collapseStart)
					{
						bucket_bott->SetBodyFixed(false);
						vibrate_bucket(t,bucket_bott);
					}
				}
				else{
					bucket->SetBodyFixed(true);
				}
		}
		else{

			if ((fmod(mphysicalSystem.GetChTime(), timeForVerticalDisplcement) < dT) &&
				(numGeneratedLayers < numLayers)){
				AddParticlesLayer(mphysicalSystem, mySmarticlesVec);
				numGeneratedLayers++;
			}

			if (smarticleType == SMART_ARMS)
			{

				for (int i = 0; i < mySmarticlesVec.size(); i++) {
					double omega = 10;
					if (tStep < 500) {
						mySmarticlesVec[i]->SetActuatorFunction(0, -omega, dT);
					}
					else 
					{
						mySmarticlesVec[i]->SetActuatorFunction(0, omega, dT);
						printf("\n");

						
					}

				}

			}
			if (t > vibrateStart){
				bucket->SetBodyFixed(false);
				vibrate_bucket(t,bucket);
			}
			else{
				bucket->SetBodyFixed(true);
			}
		}
	  
		//if (numGeneratedLayers == numLayers)
		//{
		//	//start shaking
		//}
		printf("\n");

		


		 


		
	



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
	  step_timer.Reset();
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
		char filename[100];
		sprintf(filename, "screenshot%d.bmp", tStep);
		//screenshot(filename);
	  step_timer.stop("step time");
	  std::cout << "step time: " << step_timer.GetTime("step time") << ", time passed: " << int(timeDiff)/3600 <<":"<< (int(timeDiff) % 3600) / 60 << ":" << (int(timeDiff) % 60) <<std::endl;
	  printf("num contacts %d, time %f\n", mphysicalSystem.GetNcontacts(), mphysicalSystem.GetChTime());


		FixBodies(mphysicalSystem, tStep);
		if (read_from_file)
		{
			if (bucket_exist)
				PrintFractionsAndCOM(mphysicalSystem, tStep, mySmarticlesVec, bucket->GetPos());
			else
				PrintFractionsAndCOM(mphysicalSystem, tStep, mySmarticlesVec, bucket_bott->GetPos());
				
		}
		else{
			PrintFractions(mphysicalSystem, tStep, mySmarticlesVec);
		}

	  std::cout.flush();


		CheckPointSmarticles_Write(mySmarticlesVec,
			tStep,
			mat_g,
			l_smarticle,
			w_smarticle,
			t_smarticle,
			t2_smarticle,
			collisionEnvelope,
			rho_smarticle,
			armAngle1,
			armAngle2);

  }
  for (int i = 0; i < mySmarticlesVec.size(); i++) {
	  delete mySmarticlesVec[i];

  }
  mySmarticlesVec.clear();
	simParams << "completed"<<std::endl;
  simParams.close();
  return 0;
}


