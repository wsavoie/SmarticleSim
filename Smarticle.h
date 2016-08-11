/*
 * Smarticle.h
 *
 *  Created on: May 22, 2015
 *      Author: Arman Pazouki
 */

#ifndef SMARTICLE_H_
#define SMARTICLE_H_

//#include "core/ChVector.h"
#include <iostream>
#include "assets/ChTexture.h"
#include "chrono_irrlicht/ChIrrApp.h"//changed path from unit to chrono to reflect changes in updated chrono
#include "chrono_irrlicht/ChIrrTools.h"
#include <irrlicht.h>
//#include "physics/ChSystem.h"  // Arman: take care of this later
//#include "chrono_parallel/physics/ChSystemParallel.h"
#include <memory>
#include <deque>
#include "common.h"
#include "Controller.h"

namespace chrono {


	enum MoveType { GLOBAL = 0, GUI1 = 1, GUI2 = 2, GUI3 = 3, VIB = 4, EXTRA1 = 5, EXTRA2 = 6, MIDT = 7, OT = 8 }; //IF ADDING MORE ALWAYS KEEP OT AS LAST INDEX!!!!

	// structs to attach motion to smarticles
	class Smarticle {
	public:
		int steps=0; //number of steps this smarticle has existed, updated in smarticleUpdate
		static std::shared_ptr<ChTexture> mtextureOT;
		static std::shared_ptr<ChTexture> mtextureMid;
		static std::shared_ptr<ChTexture> mtextureArm;
		bool active=false;
		static double pctActive;
		bool armBroken;
		double timeSinceLastGait = 0;
		// Construct a smarticle and add it to ChSystem.
		Smarticle(
				ChSystem* otherSystem
				//			ChSystemParallelDVI* otherSystem
				);

		// Destructor. Nothing happens
		~Smarticle();

		virtual void Properties(
				int sID,
				double other_density,
				std::shared_ptr<ChMaterialSurface> surfaceMaterial,
				double other_envelop,
				double other_l,
				double other_w,
				double other_r,
				double other_r2,
				ChVector<> pos = ChVector<>(0, 0, 0),
				ChQuaternion<> rot = QUNIT,
				double angle1 = PI_2,
				double angle2 = PI_2);

		virtual void Properties(
				int sID,
				double other_density,
				std::shared_ptr<ChMaterialSurface> surfaceMaterial,
				double other_envelop,
				double other_l,
				double other_w,
				double other_r,
				double other_r2,
				double other_omega=0,
				ChVector<> pos = ChVector<>(0, 0, 0),
				ChQuaternion<> rot = QUNIT,
				double angle1 = PI_2,
				double angle2 = PI_2);

		virtual void Properties(
			int sID,
			int mdumID,
			double other_density,
			std::shared_ptr<ChMaterialSurface> surfaceMaterial,
			double other_envelop,
			double other_l,
			double other_w,
			double other_r,
			double other_r2,
			double other_omega = 0,
			bool willVersion = true,
			ChVector<> pos = ChVector<>(0, 0, 0),
			ChQuaternion<> rot = QUNIT,
			double angle1 = PI_2,
			double angle2 = PI_2,
			double other_torThresh2=1,
			double other_angLow=0,
			double other_angHigh=120);

		double initialAng0;
		double initialAng1;
		static const size_t numEngs = 2;
		static const size_t numSegs = 3;
		//Controller armsControl (new Controller());

		virtual ChVector<> GetReactTorqueVector(int id);

		virtual double GetArmTorque(int index); 
		virtual double GetTotalTorque();

		virtual bool ChangeArmColor(double torque01, double torque12, bool LA, bool MA, bool OA);
		bool MoveOverStress(double torque01, double torque12);
		bool MoveMidStress(double torque01, double torque12);
		bool MoveLowStress(double torque01, double torque12,double timeSinceChanged);

		virtual void SetDefaultOmega(double omega);
		
		virtual void SetOmega(int idx, double momega, bool angularFreq=true);
		virtual void SetOmega(double momega, bool angularFreq = true);
		virtual void SetOmega1(double momega1, bool angularFreq = true);
		virtual void SetOmega2(double momega2, bool angularFreq = true);
		double GetActuatorOmega(int id);
		virtual double GetOmega(int index, bool angularFreq = true);
		virtual double GetOmega1(bool angularFreq = true);
		virtual double GetOmega2(bool angularFreq = true);
		double GetNextOmega(int id);
		double GetMotTorque(int id);
		virtual std::shared_ptr<ChBody> GetSmarticleBodyPointer();

		// create the smarticle by creating arms, adding joint between them, and functions
		virtual void Create();

		// get arm shared pointer

		virtual std::shared_ptr<ChBody> GetArm(int armID);

		// get joint shared pointer
		// jointID belongs to {0, 1}, i.e. the joint between 0 and 1, or between 1 and 2
		virtual std::shared_ptr<ChLinkLockRevolute> GetRevoluteJoint(int jointID);

		// get actuator function
		// actuatorID belongs to {0, 1}, i.e. the actuatorID between 0 and 1, or between 1 and 2
		virtual std::shared_ptr<ChFunction> GetActuatorFunction(int actuatorID);
		virtual void SetActuatorFunction(int actuatorID, std::shared_ptr<ChFunction> actuatorFunction);
		Controller* armsController;
		// Smarticle volume
		virtual double GetVolume();
		virtual ChVector<> Get_cm();
		virtual ChVector<> Get_InitPos();

		virtual double GetMass();
		virtual double GetDensity() {return density;};
		virtual void TransportSmarticle(ChVector<>);
		virtual void RotateSmarticle(ChQuaternion<>);
		virtual void RotateSmarticleBy(ChQuaternion<>);
		virtual void SetEdges();
		virtual void SetSpeed(ChVector<> newSpeed);


		//smarticle arm angle
		virtual void SetInitialAngles();
		virtual double GetInitialAngle(int id);
		virtual void SetAngles(double mangle1, double mangle2, bool degrees = false);
		virtual void SetAngle(std::pair<double, double> mangles, bool degrees = false);
		virtual void SetAngle(double mangle, bool degrees = false);
		virtual void SetAngle(int id, double mangle, bool degrees = false);
		virtual void SetAngle1(double mangle1, bool degrees = false);
		virtual void SetAngle2(double mangle2, bool degrees = false);

		std::vector<ChBody*> body_list;

		
		virtual int GetID();
		virtual double GetAngle(int id, bool degrees = false);
		virtual double GetAngle1(bool degrees = false);
		virtual double GetAngle2(bool degrees = false);

		//body fixing
		virtual void SetBodyFixed(bool mev);

		void MoveToAngle(double, double);
		virtual bool GetArm0OT();
		virtual bool GetArm2OT();
		void MoveLoop();
		////////////Will smarticle implementation////////////
		void AssignState(int guiState);
		static std::vector<std::pair<double, double>> global;
		static std::vector<std::pair<double, double>> gui1;//gui option 1
		static std::vector<std::pair<double, double>> gui2;//gui option 2
		static std::vector<std::pair<double, double>> gui3;//gui option 3
		static std::vector<std::pair<double, double>> midTorque;//midtorque1
		static std::vector<std::pair<double, double>> extra1; //extraGait1
		static std::vector<std::pair<double, double>> extra2; //extraGait1
		bool visualize=false;
		bool successfulMotion = false;
		bool prevSuccessful = false;
		std::vector<std::pair<double, double>> ot; //over torque

		std::vector<std::pair<double, double>> vib; //vibrate this HAS to be particle specific so cannot be static?
		
		std::vector<int> moveTypeIdxs;//this vector keeps the current values of the move types
		MoveType moveType;
		MoveType prevMoveType;
		double torqueLimit; //torqueThres^2 to avoid using sqrts
		double angLow;
		int specialState = -1;
		bool lowStressChange = true;
		double gaitLengthChangeTime=.25;

		double angHigh;
		static double distThresh;
		static unsigned int global_GUI_value;
		void CheckLTTimer(double torque1, double torque2);
		void CheckOTTimer();

		//////////////

		double activateStress = 0;;
		double LTThresh; //Low  torque threshold
		double MTThresh; //Mid  torque threshold
		double OTThresh; //Over torque threshold
		void UpdateState();
		//arm0
		ChVector<> arm0OuterEdge;
		ChVector<> arm0Edge;
		ChVector<> arm0Front;
		ChVector<> arm0Back;

		//arm1
		ChVector<> arm10Shared;
		ChVector<> arm12Shared;
		ChVector<> arm1Front;
		ChVector<> arm1Back;


		//arm2
		ChVector<> arm2OuterEdge;
		ChVector<> arm2Front;
		ChVector<> arm2Back;
		ChVector<> arm2Edge;

		//armAxes[arms][edges]
		ChVector<>armAxes[3][4];
		//armVerts[arms][points]
		ChVector<>armVerts[3][4];
		//{ { VNULL, VNULL, VNULL, VNULL },
		//{ VNULL, VNULL, VNULL, VNULL },
		//{ VNULL, VNULL, VNULL, VNULL } };
		////////////////

		std::vector<std::pair<double, double>> *mv;
		std::deque<std::tuple<double,double,double,double>> torques;
		std::deque<double> torque1;
		std::deque<double> torque2;
		std::tuple<double,double,double,double> torqueAvg;
		void updateTorqueDeque(); 
		void updateTorqueAvg(std::tuple <double, double,double,double > oldT);
		///////////////////////////////////////////////////////////

		void SetNextAngle(int id, double ang);
		void GenerateVib(double angle1, double angle2);
		double GetNextAngle(int id);
		double GetCurrAngle(int id);
		double GetExpAngle(int id);
		bool NotAtDesiredPos(int id, double ang, double exp);
		void addInterpolatedPathToVector(double arm0i, double arm2i, double arm0f,double arm2f);
		std::vector<double> linspace(double a, double b, int n);
		std::pair<double, double> populateMoveVector();
		//populateMoveVector(std::vector<std::pair<double, double>> &mglobal, std::vector<std::pair<double, double>> &mOT, std::vector<std::pair<double, double>> &mGUI1);
		bool MoveToAngle2(std::vector<std::pair<double, double>> *v, double momega1,double momega2, MoveType mtype);

		double ChooseOmegaAmount(double momega, double currAng, double destAng);
		virtual void setCurrentMoveType(MoveType newMoveType);
		void ControllerMove(int guiState, double torque01, double torque12);
		std::shared_ptr<ChLinkEngine> getLinkActuator(int id);
		double defaultOmega;
		double omegaLim;
		//////////////////////////////////////////////////////
	private:
		// create smarticle arm, set collision, surface, and mass property.
		// armID = 0 (left arm), 1 (middle arm), 2 (right arm)
		void CreateArm(
				int armID, 			// 0: left arm, 1: middle arm, 2: right arm
				double len, 			// arm length
				ChVector<> posRel, 	// relative initPosition of the arm wrt the smarticle initPos, which is the center of the center arm
									// Y-axis is parallel to the arms. Z-axis is perpendicular to smarticle plane.
				ChQuaternion<> armRelativeRot = QUNIT	// relative rotation of the arm wrt smarticle
				);

		void CreateArm2(
			int armID, 			// 0: left arm, 1: middle arm, 2: right arm
			double len, 			// arm length
			double mr,
			double mr2,
			ChVector<> posRel, 	// relative initPosition of the arm wrt the smarticle initPos, which is the center of the center arm
			// Y-axis is parallel to the arms. Z-axis is perpendicular to smarticle plane.
			ChQuaternion<> armRelativeRot = QUNIT	// relative rotation of the arm wrt smarticle
			);

		//void CreateArm1(int armID, double len, ChVector<> posRel, ChQuaternion<> armRelativeRot=QUNIT);
		void CreateJoints();
		//void CreateJoints1(ChQuaternion<>, ChQuaternion<>);
		void CreateActuators();
		//void CreateActuators1(ChQuaternion<>, ChQuaternion<>);
		

		


	public:
		// location and orientation (location of the center of the middle arm)
		ChVector<> initPos;
		ChQuaternion<> rotation;

		// length
		double l;  // arm length including the thickness
		double w;  // mid-segment length including thickness
		double r;  // in-plane half-thickness of arm
		double r2;  // off-plane  half-thickness of arm, i.e. prependicular to the object
		double angle1; //angle between center segment and outer segments radians
		double angle2; //angle between center segment and outer segments
		double volume;
		double mass;
		double collisionEnvelope;
		bool arm0OT;
		bool arm2OT;

		//double LTThresh; //Low  torque threshold
		//double MTThresh; //Mid  torque threshold
		//double OTThresh; //Over torque threshold
		// material property
		double density;
		std::shared_ptr<ChMaterialSurface> mat_g;




		
		///< pointer to the Chrono system
		ChSystem* m_system;  // Arman : take care of this later

	 private:
		double OTTimer;//timer which keeps current value of time in the OT phase
		double OTMaxTime; //time for smarticle to be in the OT phase
		bool OTRunning; //if OT is running

	
		std::vector<int> OTVal;//vector containing OT moves to switch between
		int OTValIdx;					//current index of OTVal

		double LTTimer;//timer which keeps current value of time in the LT phase
		double LTMaxTime; //time for smarticle to be in the LT phase
		bool LTRunning; //if LT is running
		std::vector<int> LTVal;//vector containing LT moves to switch between
		int LTValIdx;					//current index of LTVal

		int smarticleID;			// smarticleID is not bodyID. smarticle is composed of 3 bodies.
		int dumID;
		double jointClearance; // space at joint
		double offPlaneoffset;
		double omega1;
		double omega2;
		std::vector <double> nextOmega;
		std::vector <double> nextAngle;
		std::vector <double> currTorque;
		// bodies
	 std::shared_ptr<ChBody> arm0;	// left arm
	 std::shared_ptr<ChBody> arm1;	// middle arm
	 std::shared_ptr<ChBody> arm2;	// right arm
	 std::shared_ptr<ChBody> smarticle;
		// joints
		std::shared_ptr<ChLinkLockRevolute> link_revolute01; 	// revolute joint between arms 0 and 1
		std::shared_ptr<ChLinkLockRevolute> link_revolute12; 	// revolute joint between arms 0 and 1

		// Actuators
		std::shared_ptr<ChLinkEngine> link_actuators[numEngs];

		// joints functions
		std::shared_ptr<ChFunction> EngFunctions[numEngs];

		// assets
		std::shared_ptr<ChTexture> arm0_textureAsset;
		std::shared_ptr<ChTexture> arm1_textureAsset;
		std::shared_ptr<ChTexture> arm2_textureAsset;

	};

}

#endif /* SMARTICLE_H_ */
