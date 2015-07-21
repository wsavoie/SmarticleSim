/*
 * SmarticleU.h
 *
 *  Created on: Jun 2, 2015
 *      Author: arman
 */

#ifndef SMARTICLEU_H_
#define SMARTICLEU_H_

#include "core/ChVector.h"
#include "Smarticle.h"

namespace chrono {

class SmarticleU : public Smarticle {

public:
	// Construct a smarticleU and add it to ChSystem.
	SmarticleU(
			  ChSystem* otherSystem
			  ) : Smarticle(otherSystem) {
		volume = GetVolume();
	}

	// Destructor. Nothing happens
	~SmarticleU() {}

	// create the smarticleU by creating a u-shaped particle without joint and actuators
	virtual void Create();

	// get center of mass
	virtual ChVector<> Get_cm();
	virtual double GetVolume();
	virtual ChSharedBodyPtr GetSmarticleBodyPointer() {
		return smarticleU;
	}
	virtual void SetAngle(double mangle, bool degrees);
	virtual double GetAngle(bool degrees);

	//////////////////////
	// no support zone
	//////////////////////

	virtual ChSharedBodyPtr GetArm(int armID) {
		printf("Warning!! SmarticleU does not have independent arms\n");
		ChSharedBodyPtr tmpPtr;
		return tmpPtr;
	}

	virtual ChSharedPtr<ChLinkLockRevolute> GetRevoluteJoint(int jointID) {
		printf("Warning!! SmarticleU does not have joint\n");
		ChSharedPtr<ChLinkLockRevolute> tmpPtr;
		return tmpPtr;
	}

	virtual ChSharedPtr<ChFunction> GetActuatorFunction(int actuatorID) {
		printf("Warning!! SmarticleU does not have actuator\n");
		ChSharedPtr<ChFunction> tmpPtr;
		return tmpPtr;
	}

	virtual void SetActuatorFunction(int actuatorID, ChSharedPtr<ChFunction> actuatorFunction) {
		printf("Warning!! SmarticleU does not have actuator. Actuator cannot be set.\n");
	}
	


private:
	  // bodies
	ChSharedBodyPtr smarticleU;
};
}



#endif /* SMARTICLEU_H_ */
