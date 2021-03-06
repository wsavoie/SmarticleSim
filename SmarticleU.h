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
			std::shared_ptr<CH_SYSTEM> otherSystem
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
		virtual std::shared_ptr<ChBody> GetSmarticleBodyPointer() {
			return smarticleU;
		}



		//////////////////////
		// no support zone
		//////////////////////

		virtual std::shared_ptr<ChBody> GetArm(int armID) {
			printf("Warning!! SmarticleU does not have independent arms\n");
			std::shared_ptr<ChBody> tmpPtr;
			return tmpPtr;
		}

		virtual std::shared_ptr<ChLinkLockRevolute> GetRevoluteJoint(int jointID) {
			printf("Warning!! SmarticleU does not have joint\n");
			std::shared_ptr<ChLinkLockRevolute> tmpPtr;
			return tmpPtr;
		}

		virtual std::shared_ptr<ChFunction> GetActuatorFunction(int actuatorID) {
			printf("Warning!! SmarticleU does not have actuator\n");
			std::shared_ptr<ChFunction> tmpPtr;
			return tmpPtr;
		}

		virtual void SetActuatorFunction(int actuatorID, std::shared_ptr<ChFunction> actuatorFunction) {
			printf("Warning!! SmarticleU does not have actuator. Actuator cannot be set.\n");
		}



	private:
		// bodies
		std::shared_ptr<ChBody> smarticleU;
	};
}



#endif /* SMARTICLEU_H_ */
