/*
 * SmarticleU.cpp
 *
 *  Created on: Jun 2, 2015
 *      Author: arman
 */

#include "SmarticleU.h"
#include "chrono_utils/ChUtilsGeometry.h"
#include "chrono_utils/ChUtilsCreators.h"

using namespace chrono;


void SmarticleU::Create() {
	ChVector<> box1_dim = ChVector<>(w / 2.0, r, r2);
	ChVector<> box2_dim = ChVector<>(r2, r, l / 2.0);
	ChVector<> box3_dim = ChVector<>(r2, r, l / 2.0);
	ChVector<> box1_loc = ChVector<>(0, 0, 0);					;
	ChVector<> box2_loc = ChVector<>(-w / 2.0 + r2, 0, l / 2.0 + r2);
	ChVector<> box3_loc = ChVector<>(w / 2.0 - r2, 0, l / 2.0 + r2) ;

	ChVector<> gyr1 = utils::CalcBoxGyration(box1_dim).Get_Diag();
	ChVector<> gyr2 = utils::CalcBoxGyration(box2_dim).Get_Diag();
	ChVector<> gyr3 = utils::CalcBoxGyration(box3_dim).Get_Diag();

	double m1 = density * utils::CalcBoxVolume(box1_dim);
	double m2 = density * utils::CalcBoxVolume(box2_dim);
	double m3 = density * utils::CalcBoxVolume(box3_dim);

	mass = m1 + m2 + m3;
	cm = (m1 * box1_loc + m2 * box2_loc + m3 * box3_loc) / mass + position;

	ChVector<> rel_loc1 = box1_loc - (cm - position);
	ChVector<> rel_loc2 = box2_loc - (cm - position);
	ChVector<> rel_loc3 = box3_loc - (cm - position);

	ChVector<> mInertia;
	mInertia.x =
			m1 * (gyr1.x + ChVector<>(0, rel_loc1.y, rel_loc1.z).Length2()) +
			m2 * (gyr2.x + ChVector<>(0, rel_loc2.y, rel_loc2.z).Length2()) +
			m3 * (gyr3.x + ChVector<>(0, rel_loc3.y, rel_loc3.z).Length2()) ;

	mInertia.y =
			m1 * (gyr1.y + ChVector<>(rel_loc1.x, 0, rel_loc1.z).Length2()) +
			m2 * (gyr2.y + ChVector<>(rel_loc2.x, 0, rel_loc2.z).Length2()) +
			m3 * (gyr3.y + ChVector<>(rel_loc3.x, 0, rel_loc3.z).Length2()) ;

	mInertia.z =
			m1 * (gyr1.z + ChVector<>(rel_loc1.x, rel_loc1.y, 0).Length2()) +
			m2 * (gyr2.z + ChVector<>(rel_loc2.x, rel_loc2.y, 0).Length2()) +
			m3 * (gyr3.z + ChVector<>(rel_loc3.x, rel_loc3.y, 0).Length2()) ;

	// create body, set position and rotation, add surface property, and clear/make collision model
	smarticleU = ChSharedBodyPtr(new ChBody(new collision::ChCollisionModelParallel));

	smarticleU->SetName("smarticle_u");
	smarticleU->SetPos(position);
	smarticleU->SetRot(rotation);
    smarticleU->SetCollide(true);
    smarticleU->SetBodyFixed(false);
	smarticleU->SetMaterialSurface(mat_g);

	// create body
    smarticleU->SetMass(mass);
    smarticleU->SetInertiaXX(mInertia);

	smarticleU->GetCollisionModel()->ClearModel();
	utils::AddBoxGeometry(smarticleU.get_ptr(), box1_dim, box1_loc);
	utils::AddBoxGeometry(smarticleU.get_ptr(), box2_dim, box2_loc);
	utils::AddBoxGeometry(smarticleU.get_ptr(), box3_dim, box3_loc);
	smarticleU->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2

    smarticleU->GetCollisionModel()->BuildModel();
    m_system->AddBody(smarticleU);
}

ChVector<> SmarticleU::Get_cm() {
	return cm;
}

double SmarticleU::GetVolume() {
	return mass / density;
}