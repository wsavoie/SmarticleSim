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


enum ArmType { S_CYLINDER, S_CAPSULE, S_BOX };

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
			ArmType aType = S_CYLINDER,
			ChVector<> pos = ChVector<>(0, 0, 0),
			ChQuaternion<> rot = QUNIT);

  // Destructor. Nothing happens
  ~Smarticle();

  // create the smarticle by creating arms, adding joint between them, and functions
  void Create();

  // get arm shared pointer
  ChSharedBodyPtr GetArm(int armID);

  // get joint shared pointer
  // jointID belongs to {0, 1}, i.e. the joint between 0 and 1, or between 1 and 2
  ChSharedPtr<ChLinkLockRevolute> GetRevoluteJoint(int jointID);

  // get actuator function
  // actuatorID belongs to {0, 1}, i.e. the actuatorID between 0 and 1, or between 1 and 2
  ChSharedPtr<ChFunction> GetActuatorFunction(int actuatorID);
  void SetActuatorFunction(int actuatorID, ChSharedPtr<ChFunction> actuatorFunction);

private:
  // create smarticle arm, set collision, surface, and mass property.
  // armID = 0 (left arm), 1 (middle arm), 2 (right arm)
  void CreateArm(int armID);
  void CreateJoints();
  void CreateActuators();

 private:
  int smarticleID;

  // location and orientation (location of the center of the middle arm)
  ChVector<> position;
  ChQuaternion<> rotation;

  // length
  double l;  // arm length including the thickness
  double w;  // mid-segment length including thickness
  double r;  // radius of the cross-section, if ArmType is S_CYLINDER or S_CAPSULE
          // in-plane thickness if ArmType is S_BOX
  double r2;  // off-plane  thickness, i.e. prependicular to the object, if
              // ArmType is S_BOX

  // material property
  double density;
  ChSharedPtr<ChMaterialSurface> mat_g;

  // geometry
  ArmType armType;

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

  ///< pointer to the Chrono system
//  ChSystem* m_system;  // Arman : take care of this later
  ChSystemParallelDVI*  m_system;  // Arman : use shared ptr. also go through chrono and modify
};
}

#endif /* SMARTICLE_H_ */
