/*
 * CheckPointSmarticles.cpp
 *
 *  Created on: Jul 14, 2015
 *      Author: Arman Pazouki
 */
#include "CheckPointSmarticles.h"
#include <fstream>
#include <string>
#include "SmarticleU.h"

namespace chrono {

//====================================================================================
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
		double angle2) {


	//*******************************************************************
	int tStepsCheckPoint = 1000;

	if (tStep % tStepsCheckPoint != 0) {
		return;
	}
#ifdef _WIN32
		std::system("mkdir checkPointFiles");
#else
		std::system("mkdir -p checkPointFiles");
#endif
		if (tStep / tStepsCheckPoint == 0) {
#ifdef _WIN32
		std::system("rm checkPointFiles/*.csv");
#else
		std::system("rm checkPointFiles/*.csv");
#endif
		}

	char fileCounter[5];
	int dumNumChar = sprintf(fileCounter, "%d", int(tStep / tStepsCheckPoint) );

	char nameCheckPoint[255];
	sprintf(nameCheckPoint, "checkPointFiles/smarticles");
	strcat(nameCheckPoint, fileCounter);
	strcat(nameCheckPoint, ".csv");

	std::cout << "******************************************************************************************** " <<std::endl;
	std::ofstream outSmarticles;
	outSmarticles.open(nameCheckPoint);

	outSmarticles <<
		l_smarticle << std::endl <<
		w_smarticle << std::endl <<
		t_smarticle << std::endl <<
		t2_smarticle << std::endl <<
		collisionEnvelop << std::endl <<
		rho_smarticle << std::endl <<
		mat_g->GetKfriction() << std::endl <<
		angle1 << std::endl <<
		angle2 << std::endl<<
			'#' << std::endl;

	for (int i = 0; i < mySmarticlesVec.size(); i ++) {
		SmarticleU* mSmart = (SmarticleU*)mySmarticlesVec[i];
		outSmarticles <<
				mSmart->GetSmarticleBodyPointer()->GetPos().x << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetPos().y << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetPos().z << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetRot().e0 << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetRot().e1 << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetRot().e2 << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetRot().e3 << ", " << std::endl;
	}

	outSmarticles.close();

}

//====================================================================================
// note : for now, this is only for smarticleU and only for ChSystemParallelDVI
void CheckPointSmarticles_Read(
		CH_SYSTEM& mphysicalSystem,
		std::vector<Smarticle*> & mySmarticlesVec) {
	std::ifstream inSmarticles;
	inSmarticles.open("smarticles.csv");
	double l_smarticle, w_smarticle, t_smarticle, t2_smarticle, collisionEnvelop, friction,angle1,angle2;
	double rho_smarticle;
	ChSharedPtr<ChMaterialSurface> mat_g = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	char ddCh;
	inSmarticles >>
	l_smarticle >>
	w_smarticle >>
	t_smarticle >>
	t2_smarticle >>
	collisionEnvelop >>
	rho_smarticle >>
	friction>>
	angle1>>
	angle2;
	printf("l_smarticle %f w_smarticle %f t_smarticle %f t2_smarticle %f collisionEnvelop %f rho_smarticle %f friction %f angle1 %f angle2 %f",
			l_smarticle, w_smarticle, t_smarticle, t2_smarticle, collisionEnvelop, rho_smarticle, friction, angle1, angle2);
	mat_g->SetFriction(friction);

	ddCh = '!';
	while (ddCh != '#') {
		inSmarticles >> ddCh;
	}
	std::string ddSt;
	getline(inSmarticles, ddSt);

	int smarticleCount = 0;
	ChVector<> p3;
	ChQuaternion<> q4;
	inSmarticles >> p3.x >> ddCh >> p3.y >> ddCh >> p3.z >> ddCh >>
	q4.e0 >> ddCh >> q4.e1 >> ddCh >> q4.e2 >> ddCh >>  q4.e3 >> ddCh;
	while (inSmarticles.good()) {
		SmarticleU * smarticle0  = new SmarticleU(&mphysicalSystem);
		smarticle0->Properties(smarticleCount,
						  rho_smarticle, mat_g,
						  collisionEnvelop,
						  l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
						  p3,
						  q4);
		smarticle0->SetAngle(angle1, angle2, true);
		smarticle0->Create();
		mySmarticlesVec.push_back(smarticle0);

		smarticleCount ++;
		inSmarticles >> p3.x >> ddCh >> p3.y >> ddCh >> p3.z >> ddCh >>
		q4.e0 >> ddCh >> q4.e1 >> ddCh >> q4.e2 >> ddCh >>  q4.e3 >> ddCh;
	}

	printf("num smarticles: %d\n", mySmarticlesVec.size());

}


}



