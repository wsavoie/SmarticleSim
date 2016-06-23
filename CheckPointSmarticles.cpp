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
		std::shared_ptr<ChMaterialSurface> mat_g,
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
	#if defined(_WIN64) || defined(_WIN32)
			std::system("mkdir checkPointFiles");
	#else
			std::system("mkdir -p checkPointFiles");
	#endif
	if (tStep / tStepsCheckPoint == 0) 
	{
		#if defined(_WIN64) || defined(_WIN32)
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

	for (size_t i = 0; i < mySmarticlesVec.size(); i ++) {
		SmarticleU* mSmart = (SmarticleU*)mySmarticlesVec[i];
		outSmarticles <<
				mSmart->GetSmarticleBodyPointer()->GetPos().x << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetPos().y << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetPos().z << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetRot().e0 << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetRot().e1 << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetRot().e2 << ", " <<
				mSmart->GetSmarticleBodyPointer()->GetRot().e3 << ", " <<
				mSmart->GetAngle1(true) << ", " <<
				mSmart->GetAngle2(true) << ", " <<
				std::endl;
	}

	outSmarticles.close();

}
void CheckPointSmarticlesDynamic_Write(
	std::vector<Smarticle*> & mySmarticlesVec,
	int tStep,
	std::shared_ptr<ChMaterialSurface> mat_g,
	double l_smarticle,
	double w_smarticle,
	double t_smarticle,
	double t2_smarticle,
	double collisionEnvelop,
	double rho_smarticle) {
	//*******************************************************************
	int tStepsCheckPoint = 1000;

	if (tStep % tStepsCheckPoint != 0) {
		return;
	}
#if defined(_WIN64) || defined(_WIN32)
	std::system("mkdir checkPointFiles");
#else
	std::system("mkdir -p checkPointFiles");
#endif
	if (tStep / tStepsCheckPoint == 0)
	{
#if defined(_WIN64) || defined(_WIN32)
		std::system("rm checkPointFiles/*.csv");
#else
		std::system("rm checkPointFiles/*.csv");
#endif
	}

	char fileCounter[5];
	int dumNumChar = sprintf(fileCounter, "%d", int(tStep / tStepsCheckPoint));

	char nameCheckPoint[255];
	sprintf(nameCheckPoint, "checkPointFiles/smarticles");
	strcat(nameCheckPoint, fileCounter);
	strcat(nameCheckPoint, ".csv");

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
		Smarticle::global_GUI_value << std::endl <<
		'#' << std::endl;
		
	//inSmarticles >> p3.x >> ddCh >> p3.y >> ddCh >> p3.z >> ddCh >>
	//	q4.e0 >> ddCh >> q4.e1 >> ddCh >> q4.e2 >> ddCh >> q4.e3 >> ddCh >>
	//	angle1 >> ddCh >> angle2 >> ddCh >>
	//	dumId >> globalidx >> ddCh >> gui1idx >> ddCh >> gui2idx >> ddCh >>
	//	gui3idx >> ddCh >> prevMoveType >> ddCh >> currMoveType >> ddCh;


	for (size_t i = 0; i < mySmarticlesVec.size(); i++) {
		Smarticle* mSmart = (Smarticle*)mySmarticlesVec[i];
		double yold0 = 0;
		double yold1 = 0;
		if (mSmart->armsController)
		{
			yold0 = mSmart->armsController->yold[0];
			yold1 = mSmart->armsController->yold[1];
		}
		outSmarticles <<
			mSmart->GetSmarticleBodyPointer()->GetPos().x << ", " <<
			mSmart->GetSmarticleBodyPointer()->GetPos().y << ", " <<
			mSmart->GetSmarticleBodyPointer()->GetPos().z << ", " <<
			mSmart->GetSmarticleBodyPointer()->GetRot().e0 << ", " <<
			mSmart->GetSmarticleBodyPointer()->GetRot().e1 << ", " <<
			mSmart->GetSmarticleBodyPointer()->GetRot().e2 << ", " <<
			mSmart->GetSmarticleBodyPointer()->GetRot().e3 << ", " <<
			mSmart->GetAngle1() << ", " <<
			mSmart->GetAngle2() << ", " <<
			mSmart->GetID() << ", " <<
			mSmart->moveTypeIdxs.at(MoveType::GLOBAL) << ", " <<
			mSmart->moveTypeIdxs.at(MoveType::GUI1) << ", " <<
			mSmart->moveTypeIdxs.at(MoveType::GUI2) << ", " <<
			mSmart->moveTypeIdxs.at(MoveType::GUI3) << ", " <<
			mSmart->moveTypeIdxs.at(MoveType::VIB) << ", " <<
			mSmart->moveTypeIdxs.at(MoveType::EXTRA1) << ", " <<
			mSmart->moveTypeIdxs.at(MoveType::EXTRA2) << ", " <<
			mSmart->moveTypeIdxs.at(MoveType::MIDT) << ", " <<
			mSmart->moveTypeIdxs.at(MoveType::OT) << ", " <<
			mSmart->prevMoveType << ", " <<
			mSmart->moveType << ", " <<
			mSmart->GetOmega1() << ", "<<
			yold0 << ", " <<
			yold1 << ", " <<
			mSmart->active << ", " <<
			std::endl;
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
	auto mat_g = std::make_shared<ChMaterialSurface>();
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

	ChVector<> p3;
	ChQuaternion<> q4;
	inSmarticles >> p3.x >> ddCh >> p3.y >> ddCh >> p3.z >> ddCh >>
	q4.e0 >> ddCh >> q4.e1 >> ddCh >> q4.e2 >> ddCh >>  q4.e3 >> ddCh;
	while (inSmarticles.good()) {
		SmarticleU * smarticle0  = new SmarticleU(&mphysicalSystem);
		smarticle0->Properties(mySmarticlesVec.size(),
						  rho_smarticle, mat_g,
						  collisionEnvelop,
						  l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
						  p3,
						  q4);
		smarticle0->SetAngles(angle1, angle2, false);
		smarticle0->Create();
		mySmarticlesVec.emplace_back(smarticle0);

		inSmarticles >> p3.x >> ddCh >> p3.y >> ddCh >> p3.z >> ddCh >>
		q4.e0 >> ddCh >> q4.e1 >> ddCh >> q4.e2 >> ddCh >>  q4.e3 >> ddCh;
	}

	printf("num smarticles: %d\n", mySmarticlesVec.size());

}

void CheckPointSmarticlesDynamic_Read(
	CH_SYSTEM& mphysicalSystem,
	std::vector<Smarticle*> & mySmarticlesVec, ChIrrApp& application) {
	std::ifstream inSmarticles;
	std::pair<double, double> angPair;
	inSmarticles.open("smarticles.csv");
	double l_smarticle, w_smarticle, t_smarticle, t2_smarticle, collisionEnvelop, friction, angle1, angle2, yold0, yold1;
	int globalidx, gui1idx, gui2idx, gui3idx, dumId, vibidx, extra1idx, extra2idx, midtidx, otidx;
	bool activeSmart;
	unsigned int currMoveType, prevMoveType, gui_value;
	double rho_smarticle;
	auto mat_g = std::make_shared<ChMaterialSurface>();
	char ddCh;
	double omega;
	inSmarticles >>
		l_smarticle >>
		w_smarticle >>
		t_smarticle >>
		t2_smarticle >>
		collisionEnvelop >>
		rho_smarticle >>
		friction>>
		gui_value;
	//TODO initialize angle1 and angle2
	angle1 = 0;
	angle2 = 0;
	printf("l_smarticle %f w_smarticle %f t_smarticle %f t2_smarticle %f collisionEnvelop %f rho_smarticle %f friction %f angle1 %f angle2 %f",
		l_smarticle, w_smarticle, t_smarticle, t2_smarticle, collisionEnvelop, rho_smarticle, friction, angle1, angle2);
	mat_g->SetFriction((float) friction);
	Smarticle::global_GUI_value = gui_value;
	
	ddCh = '!';
	while (ddCh != '#') {
		inSmarticles >> ddCh;
	}
	std::string ddSt;
	getline(inSmarticles, ddSt);
	ChVector<> p3;
	ChQuaternion<> q4;
	//TODO torque threshold
	inSmarticles >> p3.x >> ddCh >> p3.y >> ddCh >> p3.z >> ddCh >>
		q4.e0 >> ddCh >> q4.e1 >> ddCh >> q4.e2 >> ddCh >> q4.e3 >> ddCh >>
		angle1 >> ddCh >> angle2 >> ddCh >> dumId >> ddCh >>
		globalidx >> ddCh >> gui1idx >> ddCh >> gui2idx >> ddCh >>
		gui3idx >> ddCh >> vibidx >> ddCh >> extra1idx >> ddCh >>
		extra2idx >> ddCh >> midtidx >> ddCh >> otidx >> ddCh >>
		prevMoveType >> ddCh >> currMoveType >> ddCh >> omega >>
		ddCh >> yold0 >> ddCh >> yold1 >> ddCh >> activeSmart >> ddCh;
	//TODO add 
	while (inSmarticles.good()) {
		
		Smarticle * smarticle0 = new Smarticle(&mphysicalSystem);
		smarticle0->active = activeSmart;
		smarticle0->Properties(mySmarticlesVec.size(), dumId,
			rho_smarticle, mat_g,
			collisionEnvelop,
			l_smarticle, w_smarticle, 0.5 * t_smarticle, 0.5 * t2_smarticle,
			omega,true,p3,
			q4,angle1*D2R,angle2*D2R);
		smarticle0->populateMoveVector();
		smarticle0->SetAngles(angle1, angle2, false);
		smarticle0->visualize = true;
	
		smarticle0->Create();

		smarticle0->vib.emplace_back(angle1*D2R, angle2*D2R);
		smarticle0->AssignState(VIB);
		smarticle0->GenerateVib(angle1*D2R, angle2*D2R);

		smarticle0->moveTypeIdxs.at(MoveType::GLOBAL) = globalidx;
		smarticle0->moveTypeIdxs.at(MoveType::GUI1)		= gui1idx;
		smarticle0->moveTypeIdxs.at(MoveType::GUI2)		= gui2idx;
		smarticle0->moveTypeIdxs.at(MoveType::GUI3)		= gui3idx;
		smarticle0->vib.emplace_back(angle1*D2R, angle2*D2R);

		if (smarticle0->active)
		{
			smarticle0->armsController->yold[0] = yold0;
			smarticle0->armsController->yold[1] = yold1;
		}
		smarticle0->setCurrentMoveType(VIB);
		smarticle0->activateStress = percentToChangeStressState;

		mySmarticlesVec.emplace_back(smarticle0);
		application.DrawAll();

		/*application.GetVideoDriver()->endScene();
		application.GetVideoDriver()->beginScene(true, true,
			video::SColor(255, 140, 161, 192));*/
#if irrlichtVisualization
		application.AssetBindAll();
		application.AssetUpdateAll();
		//application.AssetUpdate(smarticle0->GetArm(0));
		//application.AssetUpdate(smarticle0->GetArm(1));
		//application.AssetUpdate(smarticle0->GetArm(2));

#endif

		inSmarticles >> p3.x >> ddCh >> p3.y >> ddCh >> p3.z >> ddCh >>
			q4.e0 >> ddCh >> q4.e1 >> ddCh >> q4.e2 >> ddCh >> q4.e3 >> ddCh >>
			angle1 >> ddCh >> angle2 >> ddCh >> dumId >> ddCh >>
			globalidx >> ddCh >> gui1idx >> ddCh >> gui2idx >> ddCh >>
			gui3idx >> ddCh >> vibidx >> ddCh >> extra1idx >> ddCh >> 
			extra2idx >> ddCh >> midtidx >> ddCh >> otidx >> ddCh >> 
			prevMoveType >> ddCh >> currMoveType >> ddCh >> omega >> 
			ddCh >> yold0 >> ddCh >> yold1 >> ddCh >> activeSmart >> ddCh;
	}

	GetLog()<< "num smarticles:"<<  mySmarticlesVec.size();

}

}



