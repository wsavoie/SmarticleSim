/*
 * SmarticleU.cpp
 *
 *  Created on: Jun 2, 2015
 *      Author: arman
 */

#include "SmarticleU.h"
#include "utils/ChUtilsGeometry.h"
#include "utils/ChUtilsCreators.h"
//#include "physics/ChSystem.h"  // Arman: take care of this later



using namespace chrono;


void SmarticleU::Create() {
	ChVector<> box1_dim = ChVector<>(w / 2.0, r, r2);
	ChVector<> box2_dim = ChVector<>(r2, r, l / 2.0);
	ChVector<> box3_dim = ChVector<>(r2, r, l / 2.0);

	// smarticle initPos is the location of the center of the center segment
	ChVector<> box1_loc = ChVector<>(0, 0, 0);
	
	
	ChVector<> box2_loc = ChVector<>(
		(-w / 2.0 + r2) - (l / 2.0+r2)*cos(angle1),
		0,
		(l / 2.0 + r2)*sin(angle1));
	
	
	ChVector<> box3_loc = ChVector<>(
		(w / 2.0 - r2) + (l / 2.0+r2)*cos(angle2),
		0,
		(l / 2.0 + r2)*sin(angle2));

	// relative location of the boxes wrt smarticle initPos,
	ChQuaternion<> quat2 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, -angle1 +PI_2, 0));
	ChQuaternion<> quat3 = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, angle2 - PI_2, 0));


	ChVector<> gyr1 = utils::CalcBoxGyration(box1_dim,box1_loc).Get_Diag();
	ChVector<> gyr2 = utils::CalcBoxGyration(box2_dim,box2_loc,quat2).Get_Diag();
	ChVector<> gyr3 = utils::CalcBoxGyration(box3_dim,box3_loc,quat3).Get_Diag();

	double m1 = density * utils::CalcBoxVolume(box1_dim);
	double m2 = density * utils::CalcBoxVolume(box2_dim);
	double m3 = density * utils::CalcBoxVolume(box3_dim);

	mass = m1 + m2 + m3;

	ChVector<> cmRel;						// cm in reletive reference frame at smarticle initPos
	cmRel = (m1 * box1_loc + m2 * box2_loc + m3 * box3_loc) / mass;

	ChVector<> rel_loc1 = box1_loc - cmRel; // relative location wrt CM, needed for parallel axis theorem
	ChVector<> rel_loc2 = box2_loc - cmRel;
	ChVector<> rel_loc3 = box3_loc - cmRel;

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
			
	// create body, set initPos and rotation, add surface property, and clear/make collision model

		smarticleU = ChSharedBodyPtr(new ChBody);
	


//	ChVector<> cm = initPos;		// cm in abs reference frame
	smarticleU->SetName("smarticle_u");
	smarticleU->SetPos(initPos);
	smarticleU->SetRot(rotation);
    smarticleU->SetCollide(true);
    smarticleU->SetBodyFixed(false);
	smarticleU->SetMaterialSurface(mat_g);

	smarticleU->GetCollisionModel()->ClearModel();
	//smarticleU->GetCollisionModel()->SetDefaultSuggestedEnvelope(.4*r2);
	// initialize collision geometry wrt cm
	
	utils::AddBoxGeometry(smarticleU.get_ptr(), box1_dim, rel_loc1);
	utils::AddBoxGeometry(smarticleU.get_ptr(), box2_dim, rel_loc2,quat2);
	utils::AddBoxGeometry(smarticleU.get_ptr(), box3_dim, rel_loc3,quat3);
	smarticleU->GetCollisionModel()->SetFamily(2); // just decided that smarticle family is going to be 2

	smarticleU->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelop);
    smarticleU->GetCollisionModel()->BuildModel();  // this function overwrites the intertia

    // change mass and inertia property
    smarticleU->SetMass(mass);
    smarticleU->SetInertiaXX(mInertia);
    smarticleU->SetDensity(density);

    m_system->AddBody(smarticleU);
}

ChVector<> SmarticleU::Get_cm() {
	return smarticleU->GetPos();
}

double SmarticleU::GetVolume() {
	return (2 * r) * (2 * r2 )* (w + 2 * l);
}
void SmarticleU::SetAngle(double mangle1, double mangle2, bool degrees = false)// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& make sure on use I supply rads without final arg
{
	if (degrees)
	{
		angle1 = mangle1*D2R;
		angle2 = mangle2*D2R;
	}
	else
	{
		angle1 = mangle1;
		angle2 = mangle2;
	}
}
void SmarticleU::SetAngle(double mangle, bool degrees = false)
{
	if (degrees)
	{
		angle1 = mangle*D2R;
		angle2 = mangle*D2R;
	}
	else
	{
		angle1 = mangle;
		angle2 = mangle;
	}
}
void SmarticleU::SetAngle1(double mangle1, bool degrees = false)
{
	if (degrees) { angle1 = mangle1*D2R; }
	else{ angle1 = mangle1; }
}
void SmarticleU::SetAngle2(double mangle2, bool degrees = false)
{
	if (degrees) { angle2 = mangle2*D2R; }
	else{ angle2 = mangle2; }
}

double SmarticleU::GetAngle1(bool degrees = false)
{
	if (degrees)
		return angle1*R2D;
	else
		return angle1;
}
double SmarticleU::GetAngle2(bool degrees = false)
{
	if (degrees)
		return angle2*R2D;
	else
		return angle2;
}
