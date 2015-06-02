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

namespace chrono {

class Smarticle {
public:
  // Construct a smarticle and add it to ChSystem.
  Smarticle(
//		  	ChSystem* otherSystem,
			ChSystemParallelDVI* otherSystem,
			int sID,
			double other_density,
			ChSharedPtr<ChMaterialSurface> surfaceMaterial,
			double other_l,
			double other_w,
			double other_r,
			double other_r2 = 0,
			ChVector<> pos = ChVector<>(0, 0, 0),
			ChQuaternion<> rot = QUNIT);

  // Destructor. Nothing happens
  ~Smarticle();

  // create the smarticle by creating arms, adding joint between them, and functions
  void Create();

  // get arm shared pointer
  virtual ChSharedBodyPtr GetArm(int armID);

  // get joint shared pointer
  // jointID belongs to {0, 1}, i.e. the joint between 0 and 1, or between 1 and 2
  virtual ChSharedPtr<ChLinkLockRevolute> GetRevoluteJoint(int jointID);

  // get actuator function
  // actuatorID belongs to {0, 1}, i.e. the actuatorID between 0 and 1, or between 1 and 2
  virtual ChSharedPtr<ChFunction> GetActuatorFunction(int actuatorID);
  virtual void SetActuatorFunction(int actuatorID, ChSharedPtr<ChFunction> actuatorFunction);

  // Smarticle volume
  double GetVolume();
  ChVector<> Get_cm();
  double GetDensity() {return density;};

private:
  // create smarticle arm, set collision, surface, and mass property.
  // armID = 0 (left arm), 1 (middle arm), 2 (right arm)
  void CreateArm(
		  int armID, 			// 0: left arm, 1: middle arm, 2: right arm
		  double len, 			// arm length
		  ChVector<> posRel, 	// relative position of the arm wrt the smarticle position, which is the center of the center arm
								// Y-axis is parallel to the arms. Z-axis is perpendicular to smarticle plane.
		  ChQuaternion<> armRelativeRot = QUNIT	// relative rotation of the arm wrt smarticle
		  );
  void CreateJoints();
  void CreateActuators();

protected:
  // location and orientation (location of the center of the middle arm)
  ChVector<> position;
  ChQuaternion<> rotation;

  // ID
  int smarticleID;

  // length
  double l;  // arm length including the thickness
  double w;  // mid-segment length including thickness
  double r;  // in-plane half-thickness of arm
  double r2;  // off-plane  half-thickness of arm, i.e. prependicular to the object

  double volume;
  double mass;

  // material property
  double density;
  ChSharedPtr<ChMaterialSurface> mat_g;


  ///< pointer to the Chrono system
//  ChSystem* m_system;  // Arman : take care of this later
  ChSystemParallelDVI*  m_system;  // Arman : use shared ptr. also go through chrono and modify

 private:
  double jointClearance; // space at joint


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

};
}

#endif /* SMARTICLE_H_ */
