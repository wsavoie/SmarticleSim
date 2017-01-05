#include "SystemGeometry.h"

using namespace chrono;

//SystemGeometry::SystemGeometry()
//{
//	bucketID= 1;
//	mat_wall = std::make_shared<ChMaterialSurface>();
//
//	bucketTexture = std::make_shared<ChTexture>();
//	sphereTexture = std::make_shared<ChTexture>();
//	groundTexture = std::make_shared<ChTexture>();
//	floorTexture = std::make_shared<ChTexture>();
//	collisionEnvelope=1;
//	bucket_rad=1;
//	bucket_half_thick=1;
//	bucket_interior_halfDim = ChVector<>(0, 0, 0);
//	w_smarticle = 1;
//	l_smarticle = 1;
//	t_smarticle = 1;
//	t2_smarticle = 1;
//	rho_cylinder = 1;
//	wall_fric = 1;
//
//	
////##TODO make a method to redefine some params based off smartSizes or put them in common
////##TODO put rotateBucket in systemGeom?
////##TODO put vibrateBucket in systemGeom?
////##TODO put setupBucketActuator in systemGeom
////##TODO put remove set ID and set family from the individual methods and put it in regular methods
////##TODO in empty thing, difference between wall container and cyl_container
//}

class ChFunctionCustom : public ChFunction{
public:
	ChFunctionCustom(){ y = 0; y_dx = 0; y_dxdx = 0; }
	ChFunctionCustom(double m_x, double m_x_dx, double m_x_dxdx)
		: y(m_x), y_dx(m_x_dx), y_dxdx(m_x_dxdx){};
	~ChFunctionCustom(){};

	/// "Virtual" copy constructor (covariant return type).
	//virtual ChFunction_Sine* Clone() const override { return new ChFunction_Sine(*this); }

	virtual FunctionType Get_Type() const override { return FUNCT_SINE; }
		virtual ChFunctionCustom* Clone() const override { return new ChFunctionCustom(*this); };

	void Copy(ChFunction* source) {
		Set_y(source->Get_y(0));
		Set_y_dx(source->Get_y_dx(0));
		Set_y_dxdx(source->Get_y_dxdx(0));
	}
	virtual ChFunction* new_Duplicate() {
		ChFunctionCustom* m_func;
		/*	m_func = new ChFunctionCustom;
		m_func->Copy(this);*/
		return (m_func);
	}

	virtual int Get_Type() { return 1; }
	void Set_y(double x){ y = x; }
	void Set_y_dx(double x){ y_dx = x; }
	void Set_y_dxdx(double x){ y_dxdx = x; }
	virtual double Get_y(double x) const override { return y; }
	virtual double Get_y_dx(double x) const override { return y_dx; }
	virtual double Get_y_dxdx(double x) const override { return y_dxdx; }
	//virtual double Get_y(double x) const override;
	//virtual double Get_y_dx(double x) const override;
	//virtual double Get_y_dxdx(double x) const override;

private:
	double y;
	double y_dx;
	double y_dxdx;
};


SystemGeometry::SystemGeometry(ChSystem* msys, BucketType sysType, double collisionEnv,double l_smart, double w_smart, double t_smart, double t2_smart)
{
	
	mat_wall = std::make_shared<ChMaterialSurface>();


	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));
	sphereTexture->SetTextureFilename(GetChronoDataFile("sphereTexture.png"));
	groundTexture->SetTextureFilename(GetChronoDataFile("greenwhite.png"));
	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));//custom file

	//auto mmaterial = std::make_shared<ChMaterialSurface>();
	//mmaterial->SetFriction(0.4f);
	//mat_wall = mmaterial;


	sys = msys;
	sType = sysType; 
	collisionEnvelope = collisionEnv;
	l_smarticle = l_smart;
	w_smarticle = w_smart;
	t_smarticle = t_smart;
	t2_smarticle = t2_smart;


#if stapleSize
	bucket_rad = sizeScale*w_smarticle*2;
	bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, 2 * bucket_rad / sizeScale);
	boxdim = ChVector<>(.28 / 1.5 / 2, .55245 / 2, 2 * bucket_rad / 8);
	bucket_half_thick = sizeScale * .0005; //maybe too big!
#else
	bucket_rad = sizeScale*w_smarticle * 2; //3
	bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, 2 * bucket_rad / sizeScale);
	 boxdim = ChVector<>(.28 / 1.5 *2.5, .55245, 2 * bucket_rad / 8);
	 bucket_half_thick = sizeScale * .005;
#endif
	hole_size = 1 * w_smarticle;
	rho_cylinder = 1180.0;
	wall_fric = 0.8;
	bucket_ctr = ChVector<>(0, 0, 0);
	mat_wall->SetFriction(wall_fric);
	envFamily = 1;
	bucketID = 1;
	rad = 1; 
}
std::shared_ptr<ChTexture> SystemGeometry::bucketTexture = std::make_shared<ChTexture>();
std::shared_ptr<ChTexture> SystemGeometry::sphereTexture = std::make_shared<ChTexture>();
std::shared_ptr<ChTexture> SystemGeometry::groundTexture = std::make_shared<ChTexture>();
std::shared_ptr<ChTexture> SystemGeometry::floorTexture = std::make_shared<ChTexture>();
std::shared_ptr<ChMaterialSurface> SystemGeometry::mat_wall = std::make_shared<ChMaterialSurface>();
SystemGeometry::~SystemGeometry()
{
}

std::shared_ptr<ChBody> SystemGeometry::create_Box()
{
	//blahblah
	//auto mmaterial = std::make_shared<ChMaterialSurface>();
	//mmaterial->SetFriction(0.4f);
	//mmaterial->SetCompliance(0.0000005f);
	//mmaterial->SetComplianceT(0.0000005f);
	//mmaterial->SetDampingF(0.2f);


	bucket = utils::CreateBoxContainer(sys, bucketID, mat_wall,
		boxdim, bucket_half_thick, bucket_ctr, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(box_ang, 0, 0)), true, false, true, false);

	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));

	bucket->SetCollide(true);
	bucket->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	bucket_bott->GetCollisionModel()->SetFamily(envFamily);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);


	return bucket;
}
std::shared_ptr<ChBody> SystemGeometry::create_bucketShell(int num_boxes, bool overlap)
{
	double t = bucket_half_thick;
	double r = bucket_rad;
	double h = bucket_interior_halfDim.z;

	return create_EmptyCylinder(num_boxes, overlap, true, h, t, r, bucket_ctr, true, bucketTexture);
}

std::shared_ptr<ChBody> SystemGeometry::create_EmptyCylinder(int num_boxes, bool overlap, bool createVector, double half_height, double t, double r, ChVector<> pos, bool halfVis, std::shared_ptr<ChTexture> texture)
{
	auto cyl_container = std::make_shared<ChBody>();
	cyl_container->SetIdentifier(bucketID);
	//cyl_container->SetMass(mass);
	cyl_container->SetPos(pos);
	cyl_container->SetRot(QUNIT);
	cyl_container->SetBodyFixed(false);
	cyl_container->SetCollide(true);
	//double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double wallt = t / 5; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
	//double half_height = bucket_interior_halfDim.z;
	double box_side = r * 2.0 * tan(PPI / num_boxes);//side length of cyl
	double o_lap = 0;
	if (overlap){ o_lap = t * 2; }
	double ang = 2.0 * PPI / num_boxes;
	ChVector<> box_size = (0, 0, 0); //size of plates
	ChVector<> pPos = (0, 0, 0);  //position of each plate
	ChQuaternion<> quat = QUNIT; //rotation of each plate
	cyl_container->GetCollisionModel()->ClearModel();
	cyl_container->SetMaterialSurface(mat_wall);
	
	for (int i = 0; i < num_boxes; i++)
	{

		box_size = ChVector<>((box_side + wallt) / 2.0,
			wallt,
			half_height + o_lap);

		pPos = pos + ChVector<>(sin(ang * i) * (wallt + r),
			cos(ang*i)*(wallt + r),
			half_height);

		quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));

		//this is here to make half the cylinder invisible.
		bool m_visualization = false;
		if (halfVis)
		{
			if (ang*i < 3 * PPI / 4 || ang*i > 5 * PPI / 4)
			{
				m_visualization = true;
				cyl_container->AddAsset(texture);
			}
		}
		else
		{
			m_visualization = true;
			cyl_container->AddAsset(texture);
		}
		cyl_container->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		utils::AddBoxGeometry(cyl_container.get(), box_size, pPos, quat, m_visualization);

	}

	double cyl_volume = PPI*(2 * box_size.z - 2 * t)*(2 * box_size.z - 2 * t)*((2 * r + 2 * t)*(2 * r + 2 * t) - r*r) + (PPI)*(r + 2 * t)*(r + 2 * t) * 2 * t;
	cyl_container->SetMass(rho_cylinder*cyl_volume);

	//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	cyl_container->GetCollisionModel()->BuildModel();

	if (createVector)
	{
		for (int i = 0; i < num_boxes; i++)
		{
			auto wallPiece = std::make_shared<ChBody>();
			box_size = ChVector<>((box_side + wallt) / 2.0,
				wallt,
				half_height + o_lap);

			pPos = pos + ChVector<>(sin(ang * i) * (wallt + r),
				cos(ang*i)*(wallt + r),
				half_height);

			quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));

			//this is here to make half the cylinder invisible.
			bool m_visualization = false;
			if (ang*i < 3 * PPI / 4 || ang*i > 5 * PPI / 4)
			{
				m_visualization = true;
				wallPiece->AddAsset(texture);
			}
			wallPiece->GetCollisionModel()->SetEnvelope(collisionEnvelope);
			wallPiece->SetPos(pPos);
			wallPiece->GetCollisionModel()->ClearModel();
			utils::AddBoxGeometry(wallPiece.get(), box_size, ChVector<>(0, 0, 0), quat, m_visualization);
			wallPiece->SetMass(rho_cylinder*cyl_volume);
			//wallPiece->SetPos(ChVector<>(0,0,0));


			//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
			wallPiece->GetCollisionModel()->BuildModel();
			wallPiece->GetCollisionModel()->SetFamily(envFamily); ////#############
			wallPiece->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
			sys->AddBody(wallPiece);

			wallPiece->SetRot(QUNIT);
			wallPiece->SetBodyFixed(true);
			wallPiece->SetCollide(true);
			bucket_bod_vec.emplace_back(wallPiece);

		}
	}
	return cyl_container;
}

std::shared_ptr<ChBody> SystemGeometry::create_FlatHopper(ChVector<> hdim)
{
	double wallAngle = D2R * 30;

	std::shared_ptr<ChBody> flathopper = std::make_shared<ChBody>();
	double hthick = bucket_half_thick;
	double o_lap = hthick * 2;
	double hole = hole_size / 2 + hthick;

	flathopper->SetPos(bucket_ctr);

	flathopper->SetBodyFixed(false);
	flathopper->SetCollide(true);


	flathopper->GetCollisionModel()->ClearModel();
	flathopper->SetIdentifier(bucketID);
	flathopper->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	flathopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	flathopper->GetCollisionModel()->SetFamily(envFamily);
	flathopper->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	flathopper->SetMaterialSurface(mat_wall);

	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hdim.x + o_lap, 2.5*hdim.y + o_lap, hthick), ChVector<>(0, 0, -hthick)); //bottom
	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hthick, hdim.y + o_lap, hdim.z + o_lap), //-x
		ChVector<>(-hdim.x - hthick, 0, hdim.z));
	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hthick, hdim.y + o_lap, hdim.z + o_lap), //+x
		ChVector<>(hdim.x + hthick, 0, hdim.z));
	//utils::AddBoxGeometry(flathopper.get(), ChVector<>(hdim.x + o_lap, hthick, hdim.z + o_lap), //-y
	//	ChVector<>(0, -hdim.y - hthick, hdim.z));
	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hdim.x + o_lap, hthick, hdim.z + o_lap), //+y
		ChVector<>(0, hdim.y + hthick, hdim.z));

	double hopperWallWidth = (hdim.x - hole) / 2;
	double hopperWallHeight = hopperWallWidth / sin(wallAngle);
	double hopperWallAngleHeight = hopperWallHeight*cos(wallAngle);
	///////////hopper walls////////
	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hthick, hopperWallHeight, hdim.z + o_lap), //left side
		ChVector<>(hdim.x - hopperWallWidth, -hdim.y - hopperWallHeight*cos(wallAngle), hdim.z), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, wallAngle)));

	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hthick, hopperWallHeight, hdim.z + o_lap), //right side
		ChVector<>(-(hdim.x - hopperWallWidth), -hdim.y - hopperWallAngleHeight, hdim.z), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, -wallAngle)));
	flathopper->GetCollisionModel()->BuildModel();

	flathopper->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(box_ang, 0, 0)));
	sys->AddBody(flathopper);
	ChVector<> buckBottPos = ChVector<>(bucket_ctr.x, -hdim.y - 2 * hopperWallAngleHeight + hthick, hdim.z);

	///////////hopper walls////
	///width of box = hole width + 2*xwidth of side wall


	//cyl_container->SetMass(mass);


	//ChVector<> flathopperPos(0, 0, sin(box_ang)*w - t);
	//flathopper->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	//flathopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);

	//utils::AddBoxGeometry(flathopper.get(), ChVector<>(w, w, t), flathopperPos, QUNIT, true);//bottom
	//utils::AddBoxGeometry(flathopper.get(), ChVector<>(t, w, h), flathopperPos - VECT_X*(w - t) + VECT_Z*(h - t), QUNIT, true);// -x
	//utils::AddBoxGeometry(flathopper.get(), ChVector<>(t, w, h), flathopperPos + VECT_X*(w - t) + VECT_Z*(h - t), QUNIT, true);// +x
	////utils::AddBoxGeometry(flathopper.get(), ChVector<>(w, t, h), flathopperPos - VECT_Y*(w - t) + VECT_Z*(h - t), QUNIT, true);//high side -y
	//utils::AddBoxGeometry(flathopper.get(), ChVector<>(w, t, h), flathopperPos + VECT_Y*(w - t) + VECT_Z*(h - t), QUNIT, true); //down side +y


	///bucket_bott in flathopper is the +y side!
	bucket_bott->SetBodyFixed(true);
	bucket_bott->SetCollide(true);
	bucket_bott->GetCollisionModel()->ClearModel();
	bucket_bott->SetPos(bucket_ctr);
	bucket_bott->SetMaterialSurface(mat_wall);
	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));//custom file
	bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(bucket_bott.get(), ChVector<>(hdim.x + o_lap, hthick, hdim.z + o_lap), //-y
		//ChVector<>(0, -hdim.y - hthick, hdim.z));
		buckBottPos);

	bucket_bott->AddAsset(floorTexture);
	bucket_bott->GetCollisionModel()->SetFamily(envFamily);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	bucket_bott->GetCollisionModel()->BuildModel();
	bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(box_ang, 0, 0)));
	sys->AddBody(bucket_bott);

	return flathopper;
}

std::shared_ptr<ChBody> SystemGeometry::create_Hopper(double theta, bool overlap)
{
	auto hopper = std::make_shared<ChBody>();
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double r = bucket_rad;
	double ang = theta;
	double h = bucket_interior_halfDim.z;
	double sH = (h - t) / sin(ang); //side height
	//double w = hole_size / 2 + cos(ang)*h;
	double w = cos(ang)*h - t / 2;
	hopper->SetPos(bucket_ctr);
	hopper->SetRot(QUNIT);
	hopper->SetBodyFixed(true);
	hopper->SetCollide(true);



	double o_lap = 0;
	if (overlap){ o_lap = 2 * t; }

	hopper->GetCollisionModel()->ClearModel();
	hopper->SetMaterialSurface(mat_wall);

	//bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_black_bordersBlack.png"));
	//hopper->AddAsset(bucketTexture);
	hopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);


	ChVector<> hop_pos = bucket_ctr;
	hopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	//bucket_rad + 2 * bucket_half_thick, bucket_rad + 2 * bucket_half_thick, bucket_half_thick
	//utils::AddBoxGeometry(ramp.get(), ChVector<>(w, w, t), rampPos, QUNIT, true);//bucket_bottom


	utils::AddBoxGeometry(hopper.get(), ChVector<>(t, r + t, h), hop_pos - VECT_X*(-sin(ang)*h + hole_size + cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, -ang, 0)), true);// -x negAngle
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(t, r + t, h), hop_pos - VECT_X*(sin(ang)*(h)+hole_size+cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, ang, 0)), true);// -x negAngle +theta
	utils::AddBoxGeometry(hopper.get(), ChVector<>(t, r + t, h), hop_pos + VECT_X*(-sin(ang)*h + hole_size + cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, ang, 0)), true);// +x  -theta = /


	utils::AddBoxGeometry(hopper.get(), ChVector<>(2 * -sin(ang)*h + hole_size + 2 * cos(ang)*t, t, h), hop_pos - VECT_Y*(r + o_lap) + VECT_Z*(h - o_lap), QUNIT, false);//camera pos -y
	utils::AddBoxGeometry(hopper.get(), ChVector<>(2 * -sin(ang)*h + hole_size + 2 * cos(ang)*t, t, h), hop_pos + VECT_Y*(r + o_lap) + VECT_Z*(h - o_lap), QUNIT, true); //camera target +y

	hopper->GetCollisionModel()->BuildModel();
	sys->AddBody(hopper);
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(r, t, h + o_lap), ChVector<>(-r/2, 0, h+o_lap), QUNIT, true); // front plate, max_x plate

	//utils::AddBoxGeometry(hopper.get(), ChVector<>(ht, hw2 + o_lap, hh2 + o_lap), ChVector<>(-hw1 - ht, 0, h1 + hh2), QUNIT, true); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, hw2 + ht, h1 + hh2), QUNIT, true); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(hw1 + o_lap, ht, hh2 + o_lap), ChVector<>(0, -hw2 - ht, h1 + hh2), QUNIT, false); // upper part, min_x plate

	//utils::AddBoxGeometry(hopper.get(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, -hw2 - ht, hh1), QUNIT, false); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(hw1 + o_lap, ht, hh1), ChVector<>(0, hw2 + ht, hh1), QUNIT, true); // upper part, min_x plate

	//utils::AddBoxGeometry(hopper.get(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(hw3 + hh1 * tan(mtheta) + ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(mtheta, VECT_Y), true); // upper part, min_x plate
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(ht, hw2, hh1 / cos(mtheta)), ChVector<>(-hw3 - hh1 * tan(mtheta) - ht * cos(mtheta), 0, hh1 - ht * sin(mtheta)), Q_from_AngAxis(-mtheta, VECT_Y), true); // upper part, min_x plate
	//hopper->AddAsset(bucketTexture);
	return hopper;
}

std::shared_ptr<ChBody> SystemGeometry::create_Bucket_Bott()
{
	bucket_bott->SetBodyFixed(true);
	bucket_bott->SetCollide(true);
	bucket_bott->GetCollisionModel()->ClearModel();
	bucket_bott->SetPos(bucket_ctr);
	bucket_bott->SetMaterialSurface(mat_wall);
	bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(bucket_bott.get(), Vector(bucket_rad + 2 * bucket_half_thick, bucket_rad + 2 * bucket_half_thick, bucket_half_thick), Vector(0, 0, -bucket_half_thick), QUNIT, true);
	bucket_bott->AddAsset(floorTexture);

	bucket_bott->GetCollisionModel()->SetFamily(envFamily);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	bucket_bott->GetCollisionModel()->BuildModel();
	sys->AddBody(bucket_bott);
	return bucket_bott;
}

std::shared_ptr<ChBody> SystemGeometry::create_Drum(int num_boxes, bool overlap, int ridges)
{
//essentially the same as create cyl container except made it bigger and added ridges

auto drum = std::make_shared<ChBody>();

double radMult = 1.5;
drum->SetIdentifier(bucketID);
drum->SetPos(bucket_ctr);
drum->SetRot(QUNIT);
drum->SetBodyFixed(false);
drum->SetCollide(true);
double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
double wallt = t / 5.0; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
double half_height = bucket_interior_halfDim.z / (radMult * 2);
double box_side = bucket_rad * radMult * 2 * tan(PPI / num_boxes);//side length of cyl
double o_lap = 0;
if (overlap){ o_lap = t * 2; }
double ang = 2.0 * PPI / num_boxes;
ChVector<> box_size = (0, 0, 0); //size of plates
ChVector<> ridge_size = (0, 0, 0); //size of plates
ChVector<> pPos = (0, 0, 0);  //position of each plate
ChQuaternion<> quat = QUNIT; //rotation of each plate
drum->GetCollisionModel()->ClearModel();
drum->SetMaterialSurface(mat_wall);

int ridgeNum = num_boxes / ridges;
for (int i = 0; i < num_boxes; i++)
{

	box_size = ChVector<>((box_side + wallt) / 2.0,
		wallt,
		half_height + o_lap);
	pPos = bucket_ctr + ChVector<>(sin(ang * i) * (wallt + bucket_rad*radMult),
		cos(ang*i)*(wallt + bucket_rad*radMult),
		0);

	quat = Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, ang*i));

	drum->AddAsset(bucketTexture);

	drum->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(drum.get(), box_size, pPos, quat);



	if (i%ridgeNum == 0)
	{
		ridge_size = ChVector<>((box_side) / 4.0,
			w_smarticle / 8.0,
			half_height + o_lap);
		pPos = bucket_ctr + ChVector<>(sin(ang * i) * (-wallt + bucket_rad*radMult),
			cos(ang*i)*(-wallt + bucket_rad*radMult),
			0);

		drum->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		//utils::AddBoxGeometry(drum.get(), ridge_size, pPos, quat);
	}
	drum->SetRot(Q_from_AngAxis(PI_2, VECT_X));


}
//TODO add bucketVolume as global variable and set it in each function to calculate for each shape volumefraction seamlessly
//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);

drum->GetCollisionModel()->SetEnvelope(collisionEnvelope);
//front wall made invisible so we can see inside
utils::AddBoxGeometry(drum.get(), ChVector<>(wallt + bucket_rad*radMult, wallt + bucket_rad*radMult, wallt), bucket_ctr + VECT_Z*(half_height + 2 * o_lap - 2 * t), QUNIT, false);
drum->GetCollisionModel()->BuildModel();
drum->GetCollisionModel()->SetFamily(envFamily);
drum->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
sys->AddBody(drum);


bucket_bott->SetBodyFixed(true);
bucket_bott->SetCollide(true);
bucket_bott->GetCollisionModel()->ClearModel();
bucket_bott->SetPos(bucket_ctr);
bucket_bott->SetMaterialSurface(mat_wall);
floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));//custom file
bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
bucket_bott->AddAsset(bucketTexture);
utils::AddBoxGeometry(bucket_bott.get(), ChVector<>(wallt + bucket_rad*radMult, wallt + bucket_rad*radMult, wallt), bucket_ctr - VECT_Z*(half_height + 2 * o_lap - 2 * t), QUNIT);


bucket_bott->GetCollisionModel()->BuildModel();
sys->AddBody(bucket_bott);
bucket_bott->SetRot(Q_from_AngAxis(PI_2, VECT_X));
bucket_bott->GetCollisionModel()->SetFamily(envFamily);
bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

return drum;
}

std::shared_ptr<ChBody> SystemGeometry::create_Hull(double numBoxes)
{

	auto convexShape = std::make_shared<ChBody>();
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code

	//cyl_container->SetMass(mass);
	convexShape->SetPos(bucket_ctr);
	convexShape->SetRot(QUNIT);
	convexShape->SetBodyFixed(false);
	convexShape->SetCollide(true);
 
	std::vector<ChVector<>> points;

	convexShape->GetCollisionModel()->ClearModel();

	double ang = 2 * PPI / numBoxes;
	for (size_t i = 0; i < numBoxes; i++)
	{
		convexShape->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	}

	

	ChVector<> rampSize(w_smarticle * 5, w_smarticle * 5, t);


	//utils::AddBoxGeometry(convexShape.get(), rampSize, rampPos, Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(box_ang, 0, 0)), true);



	convexShape->GetCollisionModel()->BuildModel();
	sys->AddBody(convexShape);
	return convexShape;
}

void SystemGeometry::create_Container()
{
	bucket = std::make_shared<ChBody>();
	bucket_bott = std::make_shared<ChBody>();
	bucket->SetIdentifier(bucketID);
	switch (bucketType)		//http://www.engineeringtoolbox.com/friction-coefficients-d_778.html to get coefficients
	{
		case BOX:
		{
			bucket = create_Box();
			bucket_bott = create_Bucket_Bott();
			bucket_bott->SetCollide(false);
			bucket_bott->SetPos(ChVector<>(5, 5, 5));
			bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));
			break;
		}

		case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
		{
			create_Bucket_Bott();
			bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_pinkwhite.png"));
			bucket = create_bucketShell(25, true);
			break;
		}


		case FLATHOPPER:
		{
			bucket = create_FlatHopper(boxdim);
			floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
			GetLog() << "\n\nfriction: " << bucket->GetMaterialSurface()->GetKfriction() << " " << bucket->GetMaterialSurface()->GetSfriction() << "\n";
			break;
		}
		case HOPPER:
		{
			bucket = create_Hopper(box_ang, true);
			bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_black_bordersBlack.png"));
			create_Bucket_Bott();
			break;
		}
		case HULL:
		{
			bucket = create_Hull(5);
			bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
			break;
		}
		case DRUM:
		{
			bucket = create_Drum(25, true);
			bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
			break;
		}
	}

	bucket->AddAsset(bucketTexture);
	bucket->SetBodyFixed(true);
	bucket->SetCollide(true);
	bucket->GetCollisionModel()->SetFamily(envFamily);
	bucket->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

}

void SystemGeometry::create_Ground()
{

	ChVector<> boxDim = sizeScale * ChVector<>(0.1, 0.1, .002);
	ChVector<> boxLoc = sizeScale * ChVector<>(0, 0, -5.0*this->bucket_interior_halfDim.z);

	ground = std::make_shared<ChBody>();

	ground->SetMaterialSurface(mat_wall);
	ground->SetPos(boxLoc);

	// ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(true);

	ground->GetCollisionModel()->ClearModel();
	ground->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddCylinderGeometry(ground.get(), boxDim.x, boxDim.z, ChVector<>(0, 0, 0), Q_from_AngAxis(PI_2, VECT_X));
	ground->GetCollisionModel()->SetFamily(envFamily);
	ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	ground->AddAsset(groundTexture);
	ground->GetCollisionModel()->BuildModel();
	sys->AddBody(ground);
}

BucketType SystemGeometry::getGeom()
{
	return sType;
}
void SystemGeometry::create_VibrateLink(double w, double A, double t_0, std::shared_ptr<ChBody> body)
{
	vibrate_link=std::make_shared<ChLinkLockLock>();
	vibrate_link->Initialize(truss, body, ChCoordsys<>(ChVector<>(0, 0, 0)));
	double phase = -w*t_0;
	//auto vibMotion = std::make_shared<ChFunction_Sine>();
	
	
	ChFunction_Sine* vibMotion = new ChFunction_Sine();  // phase freq ampl
	vibMotion->Set_phase(phase);
	vibMotion->Set_amp(A);
	vibMotion->Set_w(w);
	vibrate_link->SetMotion_Z(vibMotion);
	vibrate_link->SetDisabled(true);
	sys->Add(vibrate_link);

}
void SystemGeometry::create_Prismatic(std::shared_ptr<ChBody> body)
{
	////
	ChFunctionCustom* pris_motion = new ChFunctionCustom(0,1.5,0);  // phase freq ampl
	//pris_motion->Set_y(.1);
	//pris_motion->Set_y_dx(.1);
	//pris_motion->Set_y_dxdx(.1);
	////
	
	auto func = std::make_shared<ChFunctionCustom>();
	pris_link = std::make_shared<ChLinkLockLock>();
	

	//link_prismatic = std::make_shared<ChLinkLockPointLine>();
	pris_engine = std::make_shared<ChLinkLinActuator>();
	auto sinefunc = std::make_shared<ChFunction_Sine>();
	switch (bucketType)		//http://www.engineeringtoolbox.com/friction-coefficients-d_778.html to get coefficients
	{
		case STRESSSTICK: case HOOKRAISE:
		{

			pris_link->Initialize(body, truss, ChCoordsys<>(ChVector<>(0, 0, 0)));
			pris_link->SetMotion_Z(pris_motion);

			//link_prismatic->Initialize(body, truss, true, ChCoordsys<>(), ChCoordsys<>(ChVector<>(0, 0, 0), QUNIT));  // set prism as vertical (default would be aligned to z, horizontal
			//pris_engine->Initialize(body, truss, true, ChCoordsys<>(body->GetPos() + ChVector<>(0, 0, -stickLen), QUNIT), ChCoordsys<>(body->GetPos() + ChVector<>(0, 0, stickLen), QUNIT));
			//func->Set_y(0);
			//func->Set_y_dx(2.5 - .5); //the value in this is always -2.5+(value specified), dont know where -2.5 comes from....
			//pris_engine->Set_dist_funct(func);
			break;
		}

	}
	//pris_engine->SetDisabled(true);
	//sys->AddLink(link_prismatic);
	////func = std::dynamic_pointer_cast<ChFunctionCustom>(pris_engine->Get_dist_funct());
	//sys->AddLink(pris_engine);

	pris_link->SetDisabled(true);
	sys->Add(pris_link);
}

void SystemGeometry::create_Truss()
{

	truss = std::make_shared<ChBody>();
	truss->SetBodyFixed(true);
	truss->GetCollisionModel()->ClearModel();
	utils::AddCylinderGeometry(truss.get(), t2_smarticle / 2, bucket_interior_halfDim.z * 1, bucket_ctr + ChVector<>(0, 0, bucket_interior_halfDim.z), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(-PI_2, 0, 0)), false);
	truss->GetCollisionModel()->BuildModel();
	truss->AddAsset(sphereTexture);
	truss->SetCollide(false);
	sys->AddBody(truss);
}

void SystemGeometry::create_Knobs(double kpr,double rows, double length)
{
	//utils should be -PI_2?
	//utils::AddCylinderGeometry(truss.get(), t2_smarticle / 2, bucket_interior_halfDim.z * 1, bucket_ctr + ChVector<>(0, 0, bucket_interior_halfDim.z), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(PI_2, 0, 0)), true);
	auto knobstick = std::make_shared<ChBody>();

	double knobRad;
	//double rad;
	double mult = 2;
	stickLen = length;
	//double sphereStickHeight = t_smarticle*mult / 2.0 * (sphereNum + 1); //shouldnt need extra 2*rad offset because of how z is defined using i below
	if (stapleSize)
	{
		rad = t_smarticle*mult*1.5;
		knobRad = t_smarticle * 2.0;
	}
	else
	{
		rad = t_smarticle*mult / 2.0;
		knobRad = t_smarticle;

	}


	//unsigned int kpr = 4;//knobs per row
	//unsigned int rows = 15; //knob per z
	double ang = 2 * PPI / kpr;
	double hp = (stickLen - 2 * rad) / rows;//height between rows
	double pOffset = PPI / kpr; //phase offset
	for (size_t row = 0; row < rows; row++)
	{
		for (size_t col = 0; col < kpr; col++)
		{
			double theta = col*ang + row*pOffset;
			//utils::AddSphereGeometry(knobstick.get(), knobRad, bucket_ctr + ChVector<>(rad*cos(col*ang + row*pOffset), rad*sin(col*ang + row*pOffset), hp*(row + 1)), Angle_to_Quat(col*ang + row*pOffset, VECT_Y), true);
			utils::AddBoxGeometry(stick.get(), ChVector<>(knobRad*1.5, rad / 4, rad / 8), bucket_ctr + ChVector<>(rad*cos(theta), rad*sin(theta), hp*(row + 1)), Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(0, 0, theta + row % 2 * PI_2)), true);
			sphereStick.emplace_back(stick);
		}
	}

	stick->GetCollisionModel()->BuildModel();
	stick->SetCollide(true);
	//stick->GetCollisionModel()->SetFamily(envFamily);
	//stick->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	stick->GetPhysicsItem()->SetIdentifier(largeID + 1);
	sys->AddBody(stick);



	//auto link_engine = std::make_shared<ChLinkEngine>();
}
void SystemGeometry::create_CentralColumnEngine(double t_0)  //$$$$$$$$$$$$$$$$$$$$$$$NOT BEING USED
{//t_0 = vibrateStart
	double knobAmp = PI_2;
	double knobW = PPI;//// rod rotating speed knobW = PI
	double knobPhase = -knobW*t_0;
	//knobcylinderfunc->Set_amp(knobAmp);
	//knobcylinderfunc->Set_w(knobW);
	//knobcylinderfunc->Set_phase(knobPhase);

	columnEngine=std::make_shared<ChLinkEngine>();
	//ChQuaternion<double> qx= Q_from_AngAxis(PI_2, VECT_Z);
	ChQuaternion<double> qx = stick->GetRot();
	columnEngine->Initialize(stick, truss, ChCoordsys<>(ChVector<>(0, 0, 0), qx));
	//columnEngine->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK); // also works as revolute support
	columnEngine->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
	//link_engine->Set_rot_funct(knobcylinderfunc);
	sys->AddLink(columnEngine);
	auto mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(columnEngine->Get_spe_funct());
	mfun2->Set_yconst(5);
	//columnEngine->SetDisabled(true);
}
void SystemGeometry::create_CentralColumn(double length)
{

	//double rad = 0;
	double mult = 4.0;
	//std::shared_ptr<ChLinkLinActuator> pris_engine;
	//std::shared_ptr<ChLinkLockPrismatic> link_prismatic;
	//auto truss = std::make_shared<ChBody>();
	//auto link_engine = std::make_shared<ChLinkEngine>();
	//auto sinefunc = std::make_shared<ChFunction_Sine>();
	//auto func2 = std::make_shared<ChFunction_Const>();


	stick = std::make_shared<ChBody>();
	stick->SetRot(QUNIT);
	stick->SetBodyFixed(true);
	stick->SetMaterialSurface(mat_wall);
	stick->AddAsset(groundTexture);
	stick->GetCollisionModel()->ClearModel();
	stick->SetMass(4);

	if (stapleSize)
	{
		rad = t_smarticle*mult / 4;
	}
	else
	{
		rad = t_smarticle*mult / 4.0;
	}

	//bucket_interior_halfDim.z*1.5;
	stickLen = length;

	int sphereNum = stickLen / (t_smarticle / 2);

	//double sphereStickHeight = t_smarticle*mult / 2.0 * (sphereNum + 1); //shouldnt need extra 2*rad offset because of how z is defined using i below
	for (size_t i = 0; i < sphereNum; i++)
	{

		stick->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		//utils::AddSphereGeometry(stick.get(), t_smarticle / 2, sys->bucket_ctr + ChVector<>(0, 0, t_smarticle*(i + 1 / 2.0)), QUNIT, true); // upper part, min_x plate
		//utils::AddSphereGeometry(stick.get(), t_smarticle / 5, sys->bucket_ctr + ChVector<>(0, 0, t_smarticle*(i + 1 / 5.0)), QUNIT, true); // upper part, min_x plate
		//utils::AddSphereGeometry(stick.get(), t2_smarticle/2.0, sys->bucket_ctr + ChVector<>(0, 0, t2_smarticle*(i + 1 /2.0)), QUNIT, true); // upper part, min_x plate

		//if you change z height between spheres, you must change sphereStickHeight above!
		utils::AddSphereGeometry(stick.get(), rad, bucket_ctr + ChVector<>(0, 0, stickLen / sphereNum * (i)), Angle_to_Quat(ANGLESET_RXYZ, ChVector<double>(0, 0, PPI)), true);
		sphereStick.emplace_back(stick);
	}

	if (bucketType == HOOKRAISE)
	{
		int hookNum;
		hookNum = 8 - stapleSize * 2;

		for (size_t i = 0; i < hookNum; i++)
		{
			if (stapleSize)
			{

				//AddBoxGeometry
				utils::AddSphereGeometry(stick.get(), rad, bucket_ctr + ChVector<>(rad*(i + 1), 0, stickLen / sphereNum), Angle_to_Quat(ANGLESET_RXYZ, ChVector<double>(PPI, 0, 0)), true);

			}
			else
			{
				utils::AddSphereGeometry(stick.get(), rad, bucket_ctr + ChVector<>(rad*(i + 1), 0, stickLen / sphereNum), Angle_to_Quat(ANGLESET_RXYZ, ChVector<double>(PPI, 0, 0)), true);
			}
			sphereStick.emplace_back(stick);
		}
	}

	stick->GetCollisionModel()->BuildModel();
	stick->GetCollisionModel()->SetFamily(bucketID);
	stick->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	stick->GetPhysicsItem()->SetIdentifier(largeID + 1);
	stick->SetCollide(true);
	sys->AddBody(stick);
	create_Truss();
}
void SystemGeometry::vibrate_body(double t,double w, double A, double t_0, std::shared_ptr<ChBody> body)  //!!!!!!!!!!!USE CREATE_VIBRATELINK INSTEAD!!!!!!!!!!!!!!
{
	//A=amplitude
	//w=omega
	//t_0 vibration starting time
	double phase = -w*t_0;
	double x_bucket = A*sin(w * t + phase);
	double xDot_bucket = A*w*cos(w * t + phase);
	double xDDot_bucket = -1*A*w*w*sin(w * t + phase);
	
	body->SetPos(body->GetPos()+ChVector<>(0, 0, x_bucket));
	body->SetPos_dt(ChVector<>(0, 0, xDot_bucket));
	body->SetPos_dtdt(ChVector<>(0, 0, xDDot_bucket));
	body->SetRot(QUNIT);
}


void SystemGeometry::setUpBucketActuator()
{
	std::shared_ptr<ChFunction_Const> mfun2; //needs to be declared outside switch
	ChVector<> pR01(0, 0, 0);
	ChQuaternion<> qx;
	bucket_actuator = std::make_shared<ChLinkEngine>();
	switch (bucketType)
	{
	case DRUM:

		qx = bucket->GetRot();
		bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), qx));
		bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
		sys->AddLink(bucket_actuator);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_spe_funct());
		mfun2->Set_yconst(0);
		//mfun2->Set_yconst(drum_omega);
		break;
	
	case KNOBCYLINDER:

		qx = bucket->GetRot();
		bucket_actuator->Initialize(stick, truss, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), qx));
		bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
		sys->AddLink(bucket_actuator);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_spe_funct());
		mfun2->Set_yconst(0);
		//mfun2->Set_yconst(drum_omega);
		break;
	case BOX:

		qx = Q_from_AngAxis(PI_2, VECT_Y);
		bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), qx));
		bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
		//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
		sys->AddLink(bucket_actuator);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
		mfun2->Set_yconst(0);
		break;

	case FLATHOPPER:

		qx = Q_from_AngAxis(PI_2, VECT_Y);
		bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), qx));
		bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
		//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
		sys->AddLink(bucket_actuator);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
		mfun2->Set_yconst(0);
		break;
	}


}
void SystemGeometry::setUpBucketActuator(ChQuaternion<double> rot)
{
	std::shared_ptr<ChFunction_Const> mfun2; //needs to be declared outside switch
	ChVector<> pR01(0, 0, 0);
	ChQuaternion<> qx;
	bucket_actuator = std::make_shared<ChLinkEngine>();
	switch (bucketType)
	{
	case DRUM:

		//qx = bucket->GetRot()
		bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), rot));
		bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
		sys->AddLink(bucket_actuator);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_spe_funct());
		mfun2->Set_yconst(0);
		//mfun2->Set_yconst(drum_omega);
		break;

	case BOX:

		/*qx = Q_from_AngAxis(PI_2, VECT_Y);*/
		bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), rot));
		bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
		//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
		sys->AddLink(bucket_actuator);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
		mfun2->Set_yconst(0);
		break;

	case FLATHOPPER:

		//qx = Q_from_AngAxis(PI_2, VECT_Y);
		bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), rot));
		bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
		//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
		sys->AddLink(bucket_actuator);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
		mfun2->Set_yconst(0);
		break;
	}


}

void SystemGeometry::rotate_body_rot(double t, std::shared_ptr<ChBody> body, std::shared_ptr<ChLinkEngine> actuator, double ang)
{
	//static std::shared_ptr<ChFunction_Const> mfun2;
	std::shared_ptr<ChFunction_Const> mfun2;
	body->SetBodyFixed(false);
	mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
	//set rotation amount equal to box_ang-currentangle, therefore no rotation if 2 values are equal
	mfun2->Set_yconst(box_ang - ang);
	//bucket->GetRot()).x
}
void SystemGeometry::rotate_body_sp(double t, std::shared_ptr<ChBody> body, std::shared_ptr<ChLinkEngine> actuator, double w)//method is called on each iteration to rotate drum at an angular velocity of drum_omega
{
		//static std::shared_ptr<ChFunction_Const> mfun2;
		//double w = f * 2 * PPI;
		std::shared_ptr<ChFunction_Const> mfun2;
		body->SetBodyFixed(false);
		mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(actuator->Get_spe_funct());
		mfun2->Set_yconst(w);
}
