/*
 * CheckPointSmarticles.h
 *
 *  Created on: Jul 14, 2015
 *      Author: Arman Pazouki
 */

#ifndef CHECKPOINTSMARTICLES_H_
#define CHECKPOINTSMARTICLES_H_

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "physics/ChSystem.h"
#include "Smarticle.h" //do we need this if smarticleU imports smarticle?
#include "SmarticleU.h"

#if USE_PARALLEL
#define CH_SYSTEM ChSystemParallelDVI
#else
#define CH_SYSTEM ChSystem
#endif

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
		double rho_smarticle,
		double angle1,
		double angle2);
void CheckPointSmarticles_Read(
	CH_SYSTEM& mphysicalSystem,
	std::vector<Smarticle*> & mySmarticlesVec);
//void CheckPointSmarticles_Read(
//		#include "SmarticleU.h",
//		std::vector<Smarticle*> & mySmarticlesVec);

}
#endif /* CHECKPOINTSMARTICLES_H_ */
