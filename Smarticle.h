/*
 * Smarticle.h
 *
 *  Created on: May 22, 2015
 *      Author: Arman Pazouki
 */

#ifndef SMARTICLE_H_
#define SMARTICLE_H_

#include "core/ChVector.h"
//#include "physics/ChSystem.h"  // Arman: take care of this later
#include "chrono_parallel/physics/ChSystemParallel.h"
#include <memory>

#ifndef true 
#define true 1
#endif

#ifndef false 
#define false 0
#endif

#define USE_PARALLEL true

namespace chrono {

	enum MotionType {SQUARE_G, CIRCLE_G, RELEASE_G, LOOP_G};

	enum MoveType { GLOBAL=0, GUI1=1, GUI2=2, GUI3=3,OT=4}; //IF ADDING MORE ALWAYS KEEP OT AS LAST INDEX!!!!
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

		virtual void SetDefaultOmega(double omega);


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
		virtual double GetDensity() {return density;};
		virtual void AddMotion(ChSharedPtr<SmarticleMotionPiece> s_motionPiece);
		//	virtual void SetCurrentMotion(ChSharedPtr<SmarticleMotionPiece> s_motionPiece); // to be implemented
		//	virtual ChSharedPtr<SmarticleMotionPiece> s_motionPiece GetCurrentMotion(); // to be implemented

		virtual void UpdateMySmarticleMotion();


		//smarticle arm angle
		virtual void SetAngle(double mangle1, double mangle2, bool degrees);
		virtual void SetAngle(double mangle, bool degrees);
		virtual void SetAngle1(double mangle1, bool degrees);
		virtual void SetAngle2(double mangle2, bool degrees);

		virtual ChSharedPtr<SmarticleMotionPiece> Get_Current_Motion();

		virtual double GetAngle1(bool degrees);
		virtual double GetAngle2(bool degrees);
		//body fixing
		virtual void SetBodyFixed(bool mev);

		void MoveToAngle(double, double);

		void MoveLoop();

		bool MoveToRange();
		void MoveSquare();
		void MoveCircle();
		void MoveRelease();
		////////////Will smarticle implementation////////////
		////vectors containing move instructions, these will be read in from a csv file, tried to make these static but wont compile...
		//TODO figure out how to make these static vars, all smarticles dont need their own copies if all share the same one
		std::vector<std::pair<double, double>> global;
		std::vector<std::pair<double, double>> gui1;//gui option 1
		std::vector<std::pair<double, double>> gui2;//gui option 2
		std::vector<std::pair<double, double>> gui3;//gui option 3
		std::vector<std::pair<double, double>> ot; //over torque
		std::vector<std::pair<double, double>> oa; //over angle

		
		std::vector<int> moveTypeIdxs;//this vector keeps the current values of the move types
		MoveType moveType;
		MoveType prevMoveType;
		double torqueThresh2; //torqueThres^2 to avoid using sqrts
		double angLow;
		double angHigh;
		double distThresh;
		///////////////////////////////////////////////////////////

		void populateMoveVector(std::vector<std::pair<double, double>> *moveVector, MoveType moveVectorType);
		bool MoveToAngle2(std::vector<std::pair<double, double>> *v, double omega, MoveType mtype);
		double Smarticle::ChooseOmegaAmount(double momega, double currAng, double destAng);
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
		void CreateJoints();
		void CreateActuators();


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

		// material property
		double density;
		ChSharedPtr<ChMaterialSurface> mat_g;




		///< pointer to the Chrono system
		ChSystem* m_system;  // Arman : take care of this later

	 private:
		// ID
		int smarticleID;			// smarticleID is not bodyID. smarticle is composed of 3 bodies.
		double jointClearance; // space at joint
		double defaultOmega;


		// bodies
	 ChSharedBodyPtr arm0;	// left arm
	 ChSharedBodyPtr arm1;	// middle arm
	 ChSharedBodyPtr arm2;	// right arm

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
