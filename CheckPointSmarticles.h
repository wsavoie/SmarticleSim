/*
 * CheckPointSmarticles.h
 *
 *  Created on: Jul 14, 2015
 *      Author: Arman Pazouki
 */

#ifndef CHECKPOINTSMARTICLES_H_
#define CHECKPOINTSMARTICLES_H_

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "Smarticle.h"
#include "SmarticleU.h"

namespace chrono {


void CheckPointSmarticles_Write(
		std::vector<Smarticle*> & mySmarticlesVec,
		int tStep,
		ChSharedPtr<ChMaterialSurface> mat_g,
		double l_smarticle,
		double w_smarticle,
		double t_smarticle,
		double t2_smarticle,
		double collisionEnvelop,
		double rho_smarticle);

void CheckPointSmarticles_Read(
		ChSystemParallelDVI & mphysicalSystem,
		std::vector<Smarticle*> & mySmarticlesVec);

}
#endif /* CHECKPOINTSMARTICLES_H_ */
