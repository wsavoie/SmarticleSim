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

#include "chrono_utils/ChUtilsCreators.h"     //Arman: why is this
#include "chrono_utils/ChUtilsInputOutput.h"  //Arman: Why is this
#include "chrono_utils/ChUtilsGenerators.h"

#include <ctime>
#include <stdlib.h>  // system, rand, srand, RAND_MAX
#include "core/ChFileutils.h" // for MakeDirectory
#include "Smarticle.h"
#include "SmarticleU.h"
#include "CheckPointSmarticles.h"

#include <memory>

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

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
SmarticleType smarticleType = SMART_ARMS;//SMART_U;
BucketType bucketType = BOX;
// =============================================================================
double Find_Max_Z(CH_SYSTEM& mphysicalSystem);
std::ofstream simParams;
ChSharedPtr<ChBody> bucket;



	double sizeScale = 1;
	double gravity = -9.81 * sizeScale;
	
	double vibration_freq = 30;
	
	double omega_bucket = 2 * CH_C_PI * vibration_freq;  // 30 Hz vibration similar to Gravish 2012, PRL
	//double vibration_amp = sizeScale * 0.00055;
	double mGamma = 2.0 * gravity;
	double vibration_amp = mGamma / (omega_bucket*omega_bucket);
	


	//double dT = std::min(0.001, 1.0 / vibration_freq / 200);;//std::min(0.0005, 1.0 / vibration_freq / 200);
	double dT = 0.001;//std::min(0.0005, 1.0 / vibration_freq / 200);
	double contact_recovery_speed = 0.2 * sizeScale;
	double tFinal = 1000;
	double vibrateStart= tFinal-5.0;

	double rho_smarticle = 7850.0 / (sizeScale * sizeScale * sizeScale);
	double rho_cylinder = 1180.0 / (sizeScale * sizeScale * sizeScale);
	ChSharedPtr<ChMaterialSurface> mat_g;
	int numLayers = 100;
	double armAngle = 90;
	double sOmega = .1;  // smarticle omega
	

	bool povray_output = true;
	int out_fps = 120;
	const std::string out_dir = "PostProcess";
	const std::string pov_dir_mbd = out_dir + "/povFilesSmarticles";

	ChVector<> bucket_ctr = ChVector<>(0,0,0);
	//ChVector<> Cbucket_interior_halfDim = sizeScale * ChVector<>(.05, .05, .025);
	double bucket_rad = sizeScale*0.044;
	ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, .010);

	
	//ChVector<> bucket_interior_halfDim = sizeScale * ChVector<>(.1, .1, .05);
	double bucket_half_thick = sizeScale * .005;
	double h = bucket_interior_halfDim.z * 2 + 2*bucket_half_thick; //from entangled paper height of available volume for smarticles in bucket
	double d = bucket_interior_halfDim.y * 2; //from entangled paper width of availble volume for smarticles in bucket

	// smarticle geometry
	double w_smarticle 	= sizeScale * 0.0117;
	double l_smarticle 	= 1 * w_smarticle; // [0.02, 1.125] * w_smarticle;
//	double t_smarticle 	= sizeScale * .00127;
//	double t2_smarticle	= sizeScale * .0005;
	double t_smarticle 	= sizeScale * .00254;
	double t2_smarticle	= sizeScale * .001;

	double collisionEnvelope = .4 * t2_smarticle;




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
	MySeed();


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
  mphysicalSystem.SetParallelThreadNumber(1);
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
		for (int i = -nX + 1; i < nX; i++) {
			for (int j = -nY + 1; j < nY; j++) {
				ChQuaternion<> myRot = ChQuaternion<>(MyRand(), MyRand(), MyRand(), MyRand());
				myRot.Normalize();


				ChVector<> myPos = ChVector<>(i * maxDim + bucket_ctr.x + MyRand()*w_smarticle - 0.5 * w_smarticle
					, j * maxDim + bucket_ctr.y + MyRand() * w_smarticle - 0.5 * w_smarticle
					, z + maxDim);

				//ChVector<> myPos = ChVector<>(i * maxDim + bucket_ctr.x, j * maxDim + bucket_ctr.y, bucket_ctr.z + 6.0 * bucket_interior_halfDim.z + 2 * bucket_half_thick);


				//ChVector<> myPos = ChVector<>(i * maxDim, j * maxDim, bucket_ctr.z + 6.0 * bucket_interior_halfDim.z + 2 * bucket_half_thick);
				// ***  added 2*bucket_half_thick to make sure stuff are initialized above bucket. Remember, bucket is inclusive, i.e. the sizes are extende 2*t from each side



//				ChSharedPtr<ChBody>(new ChBody(new collision::ChCollisionModelParallel));



				ChSharedPtr<SmarticleMotionPiece> myMotionDefault(new SmarticleMotionPiece);
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
				myMotion->SetMotionType(SQUARE_G);


				if (smarticleType == SMART_ARMS) {
					Smarticle * smarticle0 = new Smarticle(&mphysicalSystem);
					smarticle0->Properties(smarticleCount,
						rho_smarticle, mat_g,
						collisionEnvelope,
						l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
						sOmega,
						myPos,
						myRot);


					smarticle0->Create();
					smarticle0->AddMotion(myMotionDefault);
					smarticle0->AddMotion(myMotion);
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
					smarticle0->Create();
					mySmarticlesVec.push_back(smarticle0);
					if (!USE_PARALLEL) {
						smarticle0->GetSmarticleBodyPointer()->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
					}
				}
				else {
					std::cout << "Error! Smarticle type is not set correctly" << std::endl;
				}
				smarticleCount++;
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
	ChVector<> box_size = (0,0,0); //size of plates
	ChVector<> pPos = (0,0,0);  //position of each plate
	ChQuaternion<> quat=QUNIT; //rotation of each plate
	
	cyl_container->GetCollisionModel()->ClearModel();
	cyl_container->SetMaterialSurface(wallMat);
	for (int i = 0; i < num_boxes; i++)
	{

		box_size= ChVector<>((box_side + t) / 2.0,
			t,
			4*half_height + o_lap);

		pPos = bucket_ctr + ChVector<>(sin(ang * i) * (t + bucket_rad),
			cos(ang*i)*(t + bucket_rad),
			half_height);

		quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));

		//this is here to make half the cylinder invisible.
		bool m_visualization = false;
		if (ang*i < CH_C_PI  || ang*i > 3.0 * CH_C_PI / 2.0) 
		{
			m_visualization = true;
		}
		utils::AddBoxGeometry(cyl_container.get_ptr(), box_size, pPos, quat, m_visualization);

	}
	//Add ground piece
	//

	
	utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad + 2 * t, t, ChVector<>(0, 0, -t), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//add up volume of bucket and multiply by rho to get mass;
	double cyl_volume = CH_C_PI*(2 * box_size.z - 2 * t)*(2 * box_size.z - 2 * t)*((2 * bucket_rad + 2 * t)*(2 * bucket_rad + 2 * t) - bucket_rad*bucket_rad) + (CH_C_PI)*(bucket_rad + 2 * t)*(bucket_rad + 2 * t) * 2 * t;
	cyl_container->SetMass(rho_cylinder*cyl_volume);
	//utils::AddBoxGeometry(cyl_container.get_ptr(), Vector(bucket_rad, bucket_rad + t, t), Vector(0, 0, -t), QUNIT, true);

	//checks top,bottom, and middle location
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos() + Vector(0,0,2 * bucket_interior_halfDim.z + 2 * bucket_half_thick), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos(), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos() + Vector(0, 0, bucket_interior_halfDim.z), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	
	//ChVector<> bucketCtr = bucketMin + ChVector<>(0, 0, bucket_interior_halfDim.z);
	
	
	
	/*bool IsInRadial(ChVector<> pt, ChVector<> centralPt, ChVector<> rad)
	{
		bucketCtr, ChVector<>(bucket_rad, bucketMin.z, bucketMin.z + 2 * bucket_interior_halfDim.z + 2 * bucket_half_thick))
			ChVector<> dist = pt - centralPt;
		double xydist = (std::sqrt(dist.x * dist.x + dist.y + dist.y));

		if ((xydist < rad.x) && ((pt.z > rad.y) && (pt.z < rad.z))) {
			return true;
		}*/
	
	
	
	//utils::AddCylinderGeometry(cyl_container.get_ptr(), bucket_rad, 0, cyl_container->GetPos(), Q_from_AngAxis(CH_C_PI / 2, VECT_X));
		
		//, bucketMin.z + 2 * bucket_interior_halfDim.z + 2 * bucket_half_thick
		//z, bucketMin.z + 2 * bucket_interior_halfDim.z + 2 * bucket_half_thick

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

	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(hw1 + ht, 0, h1 + hh2), QUNIT, true); // uppper part, max_x plate
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(-hw1 - ht, 0, h1 + hh2), QUNIT, true); // uppper part, min_x plate
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, hw2 + ht, h1 + hh2), QUNIT, true); // uppper part, min_x plate
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, -hw2 - ht, h1 + hh2), QUNIT, true); // uppper part, min_x plate

	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, -hw2 - ht, hh1), QUNIT, true); // uppper part, min_x plate
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, hw2 + ht, hh1), QUNIT, true); // uppper part, min_x plate
	double mtheta = atan((hw1 - hw3) / h1);
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(hw3 + hh1 * tan(mtheta) + ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(mtheta, VECT_Y), true); // uppper part, min_x plate
	utils::AddBoxGeometry(cyl_container.get_ptr(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(-hw3 - hh1 * tan(mtheta) - ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(-mtheta, VECT_Y), true); // uppper part, min_x plate

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
		mat_g->SetFriction(0.4); //steel- plexiglass   (plexiglass was outer cylinder material)
	if (bucketType == BOX){
		bucket = utils::CreateBoxContainer(&mphysicalSystem, 1, mat_g, bucket_interior_halfDim, bucket_half_thick, bucket_ctr, QUNIT, true, false, true, false);
//		bucket = Create_hopper(&mphysicalSystem, mat_g, bucket_interior_halfDim.x, bucket_interior_halfDim.y, 0.5 * bucket_interior_halfDim.x, bucket_interior_halfDim.z, 2 * bucket_interior_halfDim.z,  true);

	}
	if (bucketType == CYLINDER){
		//http://www.engineeringtoolbox.com/friction-coefficients-d_778.html to get coefficients

		bucket = create_cylinder_from_blocks(25, 1, true, &mphysicalSystem, mat_g);
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
void SetEnvelopeForSystemObjects(ChSystem& mphysicalSystem) {
	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin();
	for (int i = 0; i < mphysicalSystem.Get_bodylist()->size(); i++) {
		(*myIter)->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
		myIter++;
	}

}
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
		mySmarticlesVec[i]->MoveLoop();
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
#if(!USE_PARALLEL)
  SetEnvelopeForSystemObjects(mphysicalSystem);
#endif

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();

//	ChVector<> CameraLocation = ChVector<>(0, -10, 4);
//	ChVector<> CameraLookAt = ChVector<>(0, 0, -1);
	ChVector<> CameraLocation = sizeScale * ChVector<>(-.1, -.06, .1);
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

	  int sSize1 = mySmarticlesVec.size();
	  if (  (fmod(mphysicalSystem.GetChTime(), timeForVerticalDisplcement) < dT)  &&
			  (numGeneratedLayers < numLayers) ){
		  AddParticlesLayer(mphysicalSystem, mySmarticlesVec);
		  numGeneratedLayers ++;
	  }

//  CheckPointSmarticles_Read(mphysicalSystem, mySmarticlesVec);

  printf("************** size sys %d \n", mySmarticlesVec.size());
//  for (int tStep = 0; tStep < 1; tStep++) {
	
  for (int tStep = 0; tStep < stepEnd + 1; tStep++) {
	  double t = mphysicalSystem.GetChTime();

//	  int sSize1 = mySmarticlesVec.size();
//	  if (  (fmod(mphysicalSystem.GetChTime(), timeForVerticalDisplcement) < dT)  &&
//			  (numGeneratedLayers < numLayers) ){
//		  AddParticlesLayer(mphysicalSystem, mySmarticlesVec);
//		  numGeneratedLayers ++;
//	  }
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

    UpdateSmarticles(mphysicalSystem, mySmarticlesVec);
	  time(&rawtimeCurrent);
	  double timeDiff = difftime(rawtimeCurrent, rawtime);

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

