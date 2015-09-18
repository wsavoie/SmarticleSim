/*
 * Smarticle.h
 *
 *  Created on: May 22, 2015
 *      Author: Arman Pazouki
 */

#ifndef SMARTICLE_H_
#define SMARTICLE_H_

#include "core/ChVector.h"
#include "assets/ChTexture.h"
#include "chrono_irrlicht/ChIrrApp.h"//changed path from unit to chrono to reflect changes in updated chrono
#include "chrono_irrlicht/ChIrrTools.h"
#include <irrlicht.h>
//#include "physics/ChSystem.h"  // Arman: take care of this later
#include "chrono_parallel/physics/ChSystemParallel.h"
#include <memory>

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#define USE_PARALLEL false
#define irrlichtVisualization true
namespace chrono {

	enum MotionType {SQUARE_G, CIRCLE_G, RELEASE_G, LOOP_G};

	enum MoveType { GLOBAL=0, GUI1=1, GUI2=2, GUI3=3, VIB=4, OT=5}; //IF ADDING MORE ALWAYS KEEP OT AS LAST INDEX!!!!
	// structs to attach motion to smarticles
	class JointMotion : public ChShared {
	public:
		double theta1;			// lower limit of the motion
		double theta2;			// upper limit of the motion
		double omega;			// joint angular velocity

		JointMotion() {}
		~JointMotion() {}
	};

	class SmarticleMotionPiece : public ChShared{

	public:
		JointMotion joint_01;	// joint 1 motion,
		JointMotion joint_12;	// joint 1 motion
		double timeInterval;	// time of action
		double startTime;		// start time of the motion
		double dT;
		SmarticleMotionPiece() {}
		~SmarticleMotionPiece() {}

		virtual void SetMotionType(MotionType myMotion) {motionType = myMotion;}
		virtual MotionType GetMotionType() {return motionType;}

		int motionSubSegment;
	private:
		int motionSegment;
		MotionType motionType;
	};



	class Smarticle {
	public:

		static ChSharedPtr<ChTexture> mtextureOT;
		static ChSharedPtr<ChTexture> mtextureMid;
		static ChSharedPtr<ChTexture> mtextureArm;
		bool active=false;
		static double pctActive;
		bool armBroken;
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
				ChSharedPtr<ChMaterialSurface> surfaceMaterial,
				double other_envelop,
				double other_l,
				double other_w,
				double other_r,
				double other_r2,
				ChVector<> pos = ChVector<>(0, 0, 0),
				ChQuaternion<> rot = QUNIT,
				double angle1= CH_C_PI/2,
				double angle2= CH_C_PI/2);

		virtual void Properties(
				int sID,
				double other_density,
				ChSharedPtr<ChMaterialSurface> surfaceMaterial,
				double other_envelop,
				double other_l,
				double other_w,
				double other_r,
				double other_r2,
				double other_omega=0,
				ChVector<> pos = ChVector<>(0, 0, 0),
				ChQuaternion<> rot = QUNIT,
				double angle1= CH_C_PI/2,
				double angle2= CH_C_PI/2);

		virtual void Properties(
			int sID,
			int mdumID,
			double other_density,
			ChSharedPtr<ChMaterialSurface> surfaceMaterial,
			double other_envelop,
			double other_l,
			double other_w,
			double other_r,
			double other_r2,
			double other_omega = 0,
			bool willVersion = true,
			ChVector<> pos = ChVector<>(0, 0, 0),
			ChQuaternion<> rot = QUNIT,
			double angle1 = CH_C_PI / 2,
			double angle2 = CH_C_PI / 2,
			double other_torThresh2=1,
			double other_angLow=0,
			double other_angHigh=120);


		virtual ChVector<> GetReactTorqueVectors01();
		virtual ChVector<> GetReactTorqueVectors12();
		virtual double GetReactTorqueLen01();
		virtual double GetReactTorqueLen12();

		virtual void SetDefaultOmega(double omega);
		virtual void SetOmega(double momega1, double momega2, bool angularFreq = true);
		virtual void SetOmega(double momega, bool angularFreq = true);
		virtual void SetOmega1(double momega1, bool angularFreq = true);
		virtual void SetOmega2(double momega2, bool angularFreq = true);
		virtual double GetOmega1(bool angularFreq = true);
		virtual double GetOmega2(bool angularFreq = true);
		virtual ChSharedBodyPtr GetSmarticleBodyPointer();
		// create the smarticle by creating arms, adding joint between them, and functions
		virtual void Create();

		// get arm shared pointer
		virtual ChSharedBodyPtr GetArm(int armID);

		// get joint shared pointer
		// jointID belongs to {0, 1}, i.e. the joint between 0 and 1, or between 1 and 2
		virtual ChSharedPtr<ChLinkLockRevolute> GetRevoluteJoint(int jointID);

		// get actuator function
		// actuatorID belongs to {0, 1}, i.e. the actuatorID between 0 and 1, or between 1 and 2
		virtual ChSharedPtr<ChFunction> GetActuatorFunction(int actuatorID);
		virtual void SetActuatorFunction(int actuatorID, ChSharedPtr<ChFunction> actuatorFunction);
		virtual void SetActuatorFunction(int actuatorID, double omega, double dT);
		virtual void SetActuatorFunction(int actuatorID, double omega);


		// Smarticle volume
		virtual double GetVolume();
		virtual ChVector<> Get_cm();
		virtual ChVector<> Get_InitPos();

		virtual double GetMass();
		virtual double GetDensity() {return density;};
		virtual void AddMotion(ChSharedPtr<SmarticleMotionPiece> s_motionPiece);
		//	virtual void SetCurrentMotion(ChSharedPtr<SmarticleMotionPiece> s_motionPiece); // to be implemented
		//	virtual ChSharedPtr<SmarticleMotionPiece> s_motionPiece GetCurrentMotion(); // to be implemented
		virtual void TransportSmarticle(ChVector<>);
		virtual void SetSpeed(ChVector<> newSpeed);
		virtual void UpdateMySmarticleMotion();


		//smarticle arm angle
		virtual void SetAngle(double mangle1, double mangle2, bool degrees = false);
		virtual void SetAngle(std::pair<double, double> mangles, bool degrees = false);
		virtual void SetAngle(double mangle, bool degrees = false);
		virtual void SetAngle1(double mangle1, bool degrees = false);
		virtual void SetAngle2(double mangle2, bool degrees = false);

		virtual ChSharedPtr<SmarticleMotionPiece> Get_Current_Motion();
		virtual int GetID();

		virtual double GetAngle1(bool degrees = false);
		virtual double GetAngle2(bool degrees = false);
		//body fixing
		virtual void SetBodyFixed(bool mev);

		void MoveToAngle(double, double);
		virtual bool GetArm0OT();
		virtual bool GetArm2OT();
		void MoveLoop();

		bool MoveToRange();
		void MoveSquare();
		void MoveCircle();
		void MoveRelease();
		////////////Will smarticle implementation////////////
		////vectors containing move instructions, these will be read in from a csv file, tried to make these static but wont compile...
		static std::vector<std::pair<double, double>> global;
		static std::vector<std::pair<double, double>> gui1;//gui option 1
		static std::vector<std::pair<double, double>> gui2;//gui option 2
		static std::vector<std::pair<double, double>> gui3;//gui option 3
		std::vector<std::pair<double, double>> ot; //over torque
		std::vector<std::pair<double, double>> oa; //over angle
		std::vector<std::pair<double, double>> vib; //vibrate this HAS to be particle specific so cannot be static?

		std::vector<int> moveTypeIdxs;//this vector keeps the current values of the move types
		MoveType moveType;
		MoveType prevMoveType;
		double torqueThresh2; //torqueThres^2 to avoid using sqrts
		double angLow;
		double angHigh;
		static double distThresh;
		static unsigned int global_GUI_value;
		///////////////////////////////////////////////////////////

		std::pair<double, double> populateMoveVector();
		//populateMoveVector(std::vector<std::pair<double, double>> &mglobal, std::vector<std::pair<double, double>> &mOT, std::vector<std::pair<double, double>> &mGUI1);
		bool MoveToAngle2(std::vector<std::pair<double, double>> *v, double momega1,double momega2, MoveType mtype);

		double ChooseOmegaAmount(double momega, double currAng, double destAng);
		virtual void setCurrentMoveType(MoveType newMoveType);
		void MoveLoop2(int guiState);

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

		//void CreateArm1(int armID, double len, ChVector<> posRel, ChQuaternion<> armRelativeRot=QUNIT);
		void CreateJoints();
		//void CreateJoints1(ChQuaternion<>, ChQuaternion<>);
		void CreateActuators();
		//void CreateActuators1(ChQuaternion<>, ChQuaternion<>);


	protected:
		// location and orientation (location of the center of the middle arm)
		ChVector<> initPos;
		ChQuaternion<> rotation;

		// length
		double l;  // arm length including the thickness
		double w;  // mid-segment length including thickness
		double r;  // in-plane half-thickness of arm
		double r2;  // off-plane  half-thickness of arm, i.e. prependicular to the object
		double angle1; //angle between center segment and outer segments
		double angle2; //angle between center segment and outer segments
		double volume;
		double mass;
		double collisionEnvelop;
		bool arm0OT;
		bool arm2OT;

		// material property
		double density;
		ChSharedPtr<ChMaterialSurface> mat_g;





		///< pointer to the Chrono system
		ChSystem* m_system;  // Arman : take care of this later

	 private:
		// ID
		int smarticleID;			// smarticleID is not bodyID. smarticle is composed of 3 bodies.
		int dumID;
		double jointClearance; // space at joint
		double defaultOmega;
		double omega1;
		double omega2;


		// bodies
	 ChSharedBodyPtr arm0;	// left arm
	 ChSharedBodyPtr arm1;	// middle arm
	 ChSharedBodyPtr arm2;	// right arm
	 ChSharedBodyPtr smarticle;
		// joints
		ChSharedPtr<ChLinkLockRevolute> link_revolute01; 	// revolute joint between arms 0 and 1
		ChSharedPtr<ChLinkLockRevolute> link_revolute12; 	// revolute joint between arms 0 and 1

		// Actuators
		ChSharedPtr<ChLinkEngine> link_actuator01;	// actuator joint between arms 0 and 1
		ChSharedPtr<ChLinkEngine> link_actuator12;	// actuator joint between arms 0 and 1

		// joints functions
		ChSharedPtr<ChFunction> function01;
		ChSharedPtr<ChFunction> function12;

		std::vector<ChSharedPtr<SmarticleMotionPiece>> motion_vector;
		ChSharedPtr<SmarticleMotionPiece> current_motion;





	};

}

#endif /* SMARTICLE_H_ */
