/*
 * Smarticle.cpp
 *
 *  Created on: May 22, 2015
 *      Author: Arman Pazouki
 */

#include "Smarticle.h"
#include "chrono_utils/ChUtilsGeometry.h"
#include "chrono_utils/ChUtilsCreators.h"


using namespace chrono;

Smarticle::Smarticle(
//		ChSystem* otherSystem,
		ChSystemParallelDVI* otherSystem,
		int sID,
		double other_density,
		ChSharedPtr<ChMaterialSurface> surfaceMaterial,
		double other_l,
		double other_w,
		double other_r,
		double other_r2,
		ArmType aType,
		ChVector<> pos,
		ChQuaternion<> rot) :
				m_system(otherSystem),
				smarticleID(sID),
				density(other_density),
				mat_g(surfaceMaterial),
				l(other_l),
				w(other_w),
				r(other_r),
				r2(other_r2),
				armType(aType),
				position(pos),
				rotation(rot) {


	// create 3 bodies to be connected



	// add joints

}

void Smarticle::CreateArm(int armID) { // 0: left arm, 1: middle arm, 2: right arm
	ChVector<> posRel; 	// 	relative postion of the arm wrt the smarticle position.
						//	Y-axis is parallel to the arms. Z-axis is perpendicular to smarticle plane.
	ChVector<> gyr;  	// components gyration
	double vol;			// components volume
	double jClearance = 1.1 * r;
	if (armType == S_BOX) {
		jClearance = 1.1 * r2;
	}

	ChQuaternion<> armRelativeRot = QUNIT;

	double len;
	chrono::ChSharedBodyPtr arm;
	switch (armID) {
	case 0: {
		posRel = ChVector<>(-w/2 - l/2 - jClearance, 0, 0); // Arman!! : the arms overlap with this measure. But right now we are considering ideal case
		len = l;
		arm = arm0;
	} break;
	case 1: {
		posRel = ChVector<>(0, 0, 0); // Arman!! : the arms overlap with this measure. But right now we are considering ideal case
		len = w;
		arm = arm1;
	} break;
	case 2: {
		posRel = ChVector<>(w/2 + l/2 + jClearance, 0, 0); // Arman!! : the arms overlap with this measure. But right now we are considering ideal case
		len = l;
		arm = arm2;
	} break;
	default:
		std::cout << "Error! smarticle can only have 3 arms with ids from {0, 1, 2}" << std::endl;
		break;
	}

	// NOTE: you can not combine this switch case with the next one. For some reason, it seams mass property need to be added before clearing the collision model
	switch (armType) {
	case S_CYLINDER: {
		vol = utils::CalcCylinderVolume(r, len/2.0);
		gyr = utils::CalcCylinderGyration(r, len/2.0).Get_Diag();
	} break;
	case S_CAPSULE: {
		vol = utils::CalcCapsuleVolume(r, len/2.0);
		gyr = utils::CalcCapsuleGyration(r, len/2.0).Get_Diag();
	} break;
	case S_BOX: {
		vol = utils::CalcBoxVolume(ChVector<>(len/2.0, r, r2));
		gyr = utils::CalcBoxGyration(ChVector<>(len/2.0, r, r2)).Get_Diag();
	} break;
	default:
		std::cout << "Error! smarticle shape not supported. Please choose from {cylinder, capsule, box}" << std::endl;
		break;
	}

	// create body, set position and rotation, add surface property, and clear/make collision model
	arm = chrono::ChSharedBodyPtr(new chrono::ChBody(new chrono::collision::ChCollisionModelParallel));
	ChVector<> posArm = rotation.Rotate(posRel) + position;

	arm->SetPos(posArm);
	arm->SetRot(rotation*armRelativeRot);
    arm->SetCollide(true);
    arm->SetBodyFixed(false);
	arm->SetMaterialSurface(mat_g);

	double mass = density * vol;

	// create body
    arm->SetMass(mass);
    arm->SetInertiaXX(mass * gyr);

	arm->GetCollisionModel()->ClearModel();

	// calc mass properties
	switch (armType) {
	case S_CYLINDER: {
		utils::AddCylinderGeometry(arm.get_ptr(), r, len/2.0, ChVector<>(0, 0, 0));
	} break;
	case S_CAPSULE: {
		utils::AddCapsuleGeometry(arm.get_ptr(), r, len/2.0, ChVector<>(0, 0, 0));
	} break;
	case S_BOX: {
		utils::AddBoxGeometry(arm.get_ptr(), ChVector<>(len/2.0, r, r2), ChVector<>(0, 0, 0));
	} break;
	default:
		std::cout << "Error! smarticle shape not supported. Please choose from {cylinder, capsule, box}" << std::endl;
		break;
	}

    // finalize collision
//    arm->GetCollisionModel()->SetFamily(smarticleID);
//    arm->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(smarticleID);
    arm->GetCollisionModel()->BuildModel();
    m_system->AddBody(arm);
}

void Smarticle::Create() {
	// Create Arms
	CreateArm(0);
	CreateArm(1);
	CreateArm(2);

//	link_revolute01 = ChSharedPtr<ChLinkLockRevolute>(new ChLinkLockRevolute);
//	ChVector<> pR01(-w/2 + r, 0, 0);
//	link_revolute01->Initialize(arm0, arm1,
//        ChCoordsys<>(rotation.Rotate(pR01) + position, rotation * QUNIT)); //Q_from_AngAxis(CH_C_PI / 2, VECT_Y))
//	m_system->AddLink(link_revolute01);
//
//	link_revolute12 = ChSharedPtr<ChLinkLockRevolute>(new ChLinkLockRevolute);
//	ChVector<> pR12(w/2 - r, 0, 0);
//	link_revolute12->Initialize(arm1, arm2,
//        ChCoordsys<>(rotation.Rotate(pR12) + position, rotation * QUNIT)); //Q_from_AngAxis(CH_C_PI / 2, VECT_Y))
//	m_system->AddLink(link_revolute12);
}

Smarticle::~Smarticle() {}

