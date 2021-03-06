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

class ChFunctionCustom : public ChFunction {
public:
	ChFunctionCustom() { y = 0; y_dx = 0; y_dxdx = 0; }
	ChFunctionCustom(double m_x, double m_x_dx, double m_x_dxdx)
		: y(m_x), y_dx(m_x_dx), y_dxdx(m_x_dxdx) {};
	~ChFunctionCustom() {};

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
		ChFunctionCustom* m_func = new ChFunctionCustom();
		/*	m_func = new ChFunctionCustom;
		m_func->Copy(this);*/
		return (m_func);
	}

	virtual int Get_Type() { return 1; }
	void Set_y(double x) { y = x; }
	void Set_y_dx(double x) { y_dx = x; }
	void Set_y_dxdx(double x) { y_dxdx = x; }
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


SystemGeometry::SystemGeometry(std::shared_ptr<CH_SYSTEM> msys, BucketType sysType, double collisionEnv, double l_smart, double w_smart, double t_smart, double t2_smart)
{

	mat_wall = std::make_shared<MATSURF>();


	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));
	sphereTexture->SetTextureFilename(GetChronoDataFile("sphereTexture.png"));
	groundTexture->SetTextureFilename(GetChronoDataFile("greenwhite.png"));
	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));//custom file
	hookTexture->SetTextureFilename(GetChronoDataFile("cubetexture_black_bordersBlack.png"));//custom file
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
	hookVol = 0; //define when creating hook
	topHookID = 9998;
	bottomHookID = 9999;
	
	//initialize vars
	actuationStart = 0;
	omega_bucket = 0;
	actuation_amp = 0;
	actuationSpd = 0;
	prismaticState = -1;

#if stapleSize
	bucket_rad = sizeScale*w_smarticle * 2;
	bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, 2 * bucket_rad / sizeScale);
	boxdim = ChVector<>(.18 / 1.5 / 2, .25245 / 2, 2 * bucket_rad / 8);
	bucket_half_thick = sizeScale * .005/2; //maybe too big!
#else
	bucket_rad = sizeScale*w_smarticle * 2; //3
	bucket_interior_halfDim = sizeScale * ChVector<>(bucket_rad, bucket_rad, 2 * bucket_rad / sizeScale);
	boxdim = ChVector<>(.28 / 1.5 *2.5, .55245, 2 * bucket_rad / 8);
	bucket_half_thick = sizeScale * .005;
#endif
	hole_size = 1 * w_smarticle;
	rho_cylinder = 1180.0;
	wall_fric = 0.4;
	bucket_ctr = ChVector<>(0, 0, 0);
	mat_wall->SetFriction(wall_fric);
	envFamily = 1;
	bucketID = 133711;
	rad = 1;
}
std::shared_ptr<ChTexture> SystemGeometry::bucketTexture = std::make_shared<ChTexture>();
std::shared_ptr<ChTexture> SystemGeometry::sphereTexture = std::make_shared<ChTexture>();
std::shared_ptr<ChTexture> SystemGeometry::groundTexture = std::make_shared<ChTexture>();
std::shared_ptr<ChTexture> SystemGeometry::hookTexture = std::make_shared<ChTexture>();
std::shared_ptr<ChTexture> SystemGeometry::floorTexture = std::make_shared<ChTexture>();
std::shared_ptr<MATSURF> SystemGeometry::mat_wall = std::make_shared<MATSURF>();
SystemGeometry::~SystemGeometry()
{
}
std::shared_ptr<ChBody> SystemGeometry::create_Maze()
{
	auto cc = std::make_shared<ChBody>();
	ChVector<> sz;
	ChVector<> pos;
	double h = bucket_half_thick * 5;
	ChVector<> ip = bucket->GetPos(); //initial pos
	std::shared_ptr<ChTexture> mazeText = std::make_shared<ChTexture>();
	mazeText->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
	cc->AddAsset(mazeText);
	cc->SetPos(ip);
	cc->SetRot(QUNIT);
	cc->SetBodyFixed(true);
	cc->SetCollide(true);
	cc->SetMaterialSurface(mat_wall);
	cc->SetName("MAZE WALL");

	double LW = 0.6;
	double lFL = .6;

	//top first lane
	cc->GetCollisionModel()->ClearModel();
	sz = ChVector<>(lFL, h, h);
	pos = ip + (ChVector<>(-lFL / 2 + h, LW / 2, h*1.5));
	cc->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(cc.get(), sz, pos, QUNIT, true);

	//bott first lane
	sz = ChVector<>(lFL + LW / 2 + h, h, h);
	pos = ip + (ChVector<>(-lFL / 2 - LW / 2, -LW / 2, h*1.5));
	cc->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(cc.get(), sz, pos, QUNIT, true);

	//back first lane
	sz = ChVector<>(h, LW / 2.0, h);
	pos = ip + (ChVector<>(lFL / 2.0, 0, h*1.5));
	cc->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(cc.get(), sz, pos, QUNIT, true);


	ChVector<> L2P = ip + ChVector<>(-lFL - LW / 2 + h, LW + 3 * h, 0);

	lFL = .8;

	//left second lane
	sz = ChVector<>(h, lFL / 2, h);
	pos = L2P + (ChVector<>(0, 0, h*1.5));
	cc->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(cc.get(), sz, pos, QUNIT, true);

	//right second lane
	sz = ChVector<>(h, lFL / 2 + LW, h);
	pos = L2P + (ChVector<>(-LW - 2 * h, 0, h*1.5));
	cc->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(cc.get(), sz, pos, QUNIT, true);

	lFL = .8;

	//bott third lane

	sz = ChVector<>(lFL, h, h);
	pos = ip + (ChVector<>(-lFL / 2 + LW / 2, 2 * lFL - LW + 2 * h, h*1.5));
	cc->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(cc.get(), sz, pos, QUNIT, true);

	//top third lane
	sz = ChVector<>(lFL + LW / 2 + h, h, h);
	pos = ip + (ChVector<>(-lFL / 2 - h, 2 * lFL + 4 * h, h*1.5));
	cc->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(cc.get(), sz, pos, QUNIT, true);



	cc->GetCollisionModel()->BuildModel();
	sys->AddBody(cc);
	return cc;
}

std::shared_ptr<ChBody> SystemGeometry::create_BoxBig()
{
	//blahblah
	//auto mmaterial = std::make_shared<ChMaterialSurface>();
	//mmaterial->SetFriction(0.4f);
	//mmaterial->SetCompliance(0.0000005f);
	//mmaterial->SetComplianceT(0.0000005f);
	//mmaterial->SetDampingF(0.2f);

	boxdim = ChVector<>(.28 / 1.5 * 10, .55245 * 4, 2 * bucket_rad / 8);
	bucket = utils::CreateBoxContainer(sys.get(), bucketID, mat_wall,
		boxdim, bucket_half_thick, bucket_ctr, Angle_to_Quat(ANGLE, ChVector<>(box_ang, 0, 0)), true, false, true, false);


	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));

	bucket->SetCollide(true);
	bucket->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	bucket_bott->GetCollisionModel()->SetFamily(envFamily);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);


	return bucket;
}

std::shared_ptr<ChBody> SystemGeometry::create_Box()
{
	//blahblah
	//auto mmaterial = std::make_shared<ChMaterialSurface>();
	//mmaterial->SetFriction(0.4f);
	//mmaterial->SetCompliance(0.0000005f);
	//mmaterial->SetComplianceT(0.0000005f);
	//mmaterial->SetDampingF(0.2f);


	bucket = utils::CreateBoxContainer(sys.get(), bucketID, mat_wall,
		boxdim, bucket_half_thick, bucket_ctr, Angle_to_Quat(ANGLE, ChVector<>(box_ang, 0, 0)), true, false, true, false);


	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));

	bucket->SetCollide(true);
	bucket->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	bucket_bott->GetCollisionModel()->SetFamily(envFamily);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);


	return bucket;
}
std::shared_ptr<ChBody> SystemGeometry::create_Box2()
{
	//blahblah
	//auto mmaterial = std::make_shared<ChMaterialSurface>();
	//mmaterial->SetFriction(0.4f);
	//mmaterial->SetCompliance(0.0000005f);
	//mmaterial->SetComplianceT(0.0000005f);
	//mmaterial->SetDampingF(0.2f);

	ChVector<> box_size = (0, 0, 0); //size of plates


	boxdim= ChVector<>(bucket_rad, bucket_rad, bucket_interior_halfDim.x());
	double o_lap = bucket_half_thick * 2;
	double boxBott = boxdim.y() - bucket_half_thick;
	bucket->GetCollisionModel()->ClearModel();
	bucket->GetCollisionModel()->SetEnvelope(collisionEnvelope);

	//utils::AddBoxGeometry(bucket.get(), ChVector<>(boxdim.x() + o_lap, boxdim.y() + o_lap, bucket_half_thick), ChVector<>(0, 0, -bucket_half_thick),QUNIT,false); //bottom
	utils::AddBoxGeometry(bucket.get(), ChVector<>(bucket_half_thick, boxdim.y() + o_lap, boxdim.z() + o_lap), ChVector<>(-boxdim.x() - bucket_half_thick, 0, boxdim.z()));//right wall
	utils::AddBoxGeometry(bucket.get(), ChVector<>(bucket_half_thick, boxdim.y() + o_lap, boxdim.z() + o_lap), ChVector<>(boxdim.x() + bucket_half_thick, 0, boxdim.z())); //left wall
	utils::AddBoxGeometry(bucket.get(), ChVector<>(boxdim.x() + o_lap, bucket_half_thick, boxdim.z() + o_lap), ChVector<>(0, -boxdim.y() - bucket_half_thick, boxdim.z()), QUNIT, false); //front wall
	utils::AddBoxGeometry(bucket.get(), ChVector<>(boxdim.x() + o_lap, bucket_half_thick, boxdim.z() + o_lap), ChVector<>(0, boxdim.y() + bucket_half_thick, boxdim.z()));//back wall

	//utils::AddBoxGeometry(bucket.get(), ChVector<>(boxdim.x() + o_lap, boxdim.z() + o_lap, bucket_half_thick), ChVector<>(0, boxdim.z(), boxdim.y() + bucket_half_thick+ boxBott));
	bucket->GetCollisionModel()->BuildModel();
	//id,mat,hdim,hthick,pos,rot,collide,y_up,overlap,closed)
	//bucket = utils::CreateBoxContainer(sys.get(), bucketID, mat_wall,
	//	boxdim, bucket_half_thick, bucket_ctr, Angle_to_Quat(ANGLE, ChVector<>(box_ang, 0, 0)), true, false, true, false);


	bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));

	bucket->SetCollide(true);
	//bucket->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	bucket_bott->GetCollisionModel()->SetFamily(envFamily);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	sys->AddBody(bucket);

	return bucket;
}

void chrono::SystemGeometry::hookraise()
{
	static std::shared_ptr<ChLinkLinActuator> p_engine = std::make_shared<ChLinkLinActuator>();
	static bool setup = false;
	static bool finishedActuation = false;
	double t = sys->GetChTime();
	if (!setup)
	{
		
		p_engine->Initialize(bucket_bott, topHook, ChCoordsys<>(topHook->GetPos() + (0, 0, .2), QUNIT));
		setup = true;
	}

	

	if (t > actuationStart)
	{
		if (!finishedActuation)
		{
			prismaticState = 0;
			topHook->SetBodyFixed(false);
			//bucket_actuator->SetDisabled(false);
			p_engine->SetDisabled(false);
			link_prismatic->SetDisabled(false);
			auto func = p_engine->Get_dist_funct();

			auto pris_motion = std::make_shared<ChFunctionCustom>(0, 0, 0);
			p_engine->Set_dist_funct(pris_motion);
			GetLog() << "hello";
			//auto mfun = std::dynamic_pointer_cast<ChFunction_Const>(pris_engine->Get_dist_funct());
			//mfun->Set_yconst(actuationSpd*dT + mfun->Get_yconst());
			//
			//
			//GetLog()<<mfun->Get_y(t)<<nl;
			//GetLog()<<mfun->Get_y(t)<<","<< mfun->Get_y_dx(t) <<nl;
			//mfun->Set_y(actuationSpd * dT + mfun->Get_y(t));
			//auto mfun = std::dynamic_pointer_cast<ChFunction_Const>(pris_engine->Get_dist_funct());
			//mfun->Set_yconst(-0.0001 + mfun->Get_yconst());
			//GetLog() << pris_engine->GetDist_dt()<<nl;
			if (topHook->GetPos().z() >= (actuation_amp*.995 + topOrigHeight))
			{
				//auto pris_motion = std::make_shared<ChFunctionCustom>(0, 0 + .1, 0);  // there is an hidden -0.1 offset so 0.1 keeps it from working
				//pris_engine->Set_dist_funct(pris_motion);
				finishedActuation = true;
				prismaticState = 1;
				p_engine->SetDisabled(true);
				topHook->SetBodyFixed(true);
			}
		}
		else {
			//GetLog() <<"stopped moving";
			prismaticState = 1;
		}
	}
}
void chrono::SystemGeometry::performActuation()
{
	
	double t=sys->GetChTime();
	///add method about system actuation
	if (bucketType == DRUM)
	{
		bucket_bott->SetBodyFixed(true);
		rotate_body_sp(t, bucket, bucket_actuator, drum_omega);
	}
	if (bucketType == BOX)
	{
		bucket_bott->SetBodyFixed(true);
		rotate_body_rot(t, bucket, bucket_actuator, Quat_to_Angle(ANGLE, bucket->GetRot()).x());
	}
	//vibration movement
	if (t > actuationStart && t < actuationStart + 5)
	{
		//GetLog() << "actuation started!" << nl;

		
			switch (bucketType)
			{

				case HOOKRAISE2:
				{
					//topHook->SetPos_dt(VNULL);
					static bool finishedActuation = false;

					if (!finishedActuation)
					{
						prismaticState = 0;
						topHook->SetBodyFixed(false);
						//bucket_actuator->SetDisabled(false);
						pris_engine->SetDisabled(false);
						link_prismatic->SetDisabled(false);
						if (auto mfun = std::dynamic_pointer_cast<ChFunction_Const>(pris_engine->Get_dist_funct()))
							mfun->Set_yconst(-actuationSpd*dT + mfun->Get_yconst());

						if (topHook->GetPos().z() >= (actuation_amp*.995 + topOrigHeight))
						{

							//auto pris_motion = std::make_shared<ChFunctionCustom>(0, 0 + .1, 0);  // there is an hidden -0.1 offset so 0.1 keeps it from working
							//pris_engine->Set_dist_funct(pris_motion);
							finishedActuation = true;
							prismaticState = 1;
							//pris_engine->SetDisabled(true);
							//topHook->SetBodyFixed(true);
						}
					}
					else {
						//GetLog() <<"stopped moving";
						prismaticState = 1;
						//GetLog() << "posdt" << topHook->GetPos_dt() << "trussdt" << truss->GetPos_dt() << nl;
					}
					break;
				}
				case HOOKFRACTURE:
				{
					//topHook->SetPos_dt(VNULL);
					static bool finishedActuation = false;

					if (!finishedActuation)
					{
						prismaticState = 0;
						topHook->SetBodyFixed(false);
						//bucket_actuator->SetDisabled(false);
						pris_engine->SetDisabled(false);
						link_prismatic->SetDisabled(false);
						if (auto mfun = std::dynamic_pointer_cast<ChFunction_Const>(pris_engine->Get_dist_funct()))
							mfun->Set_yconst(-actuationSpd*dT + mfun->Get_yconst());

						if (topHook->GetPos().z() >= (actuation_amp*.995 + topOrigHeight))
						{

							//auto pris_motion = std::make_shared<ChFunctionCustom>(0, 0 + .1, 0);  // there is an hidden -0.1 offset so 0.1 keeps it from working
							//pris_engine->Set_dist_funct(pris_motion);
							finishedActuation = true;
							prismaticState = 1;
							//pris_engine->SetDisabled(true);
							//topHook->SetBodyFixed(true);
						}
					}
					else {
						//GetLog() <<"stopped moving";
						prismaticState = 1;
						//GetLog() << "posdt" << topHook->GetPos_dt() << "trussdt" << truss->GetPos_dt() << nl;
					}
					break;
				}
			case HOOKRAISE: case STRESSSTICK:
			{
				stick->SetBodyFixed(false);
				pris_link->SetDisabled(false);

				if (pris_engine->IsDisabled())
				{
					stick->SetBodyFixed(false);
					pris_engine->SetDisabled(false);

				}
				pris_engine->GetDist_dt();
				break;
			}
			case CYLINDER:
			{
				bucket_bott->SetBodyFixed(false);
				vibrate_link->SetDisabled(false);
				break;
			}
			case KNOBCYLINDER:
			{
				double rotSpeed = 2; //rads/sec
				bucket_actuator->SetDisabled(false);
				stick->SetBodyFixed(false);
				rotate_body_sp(t, stick, bucket_actuator, PPI);
				break;
			}
			case HOPPER:
			{
				bucket->SetBodyFixed(false);
				vibrate_link->SetDisabled(false);
				break;
			}
			case FLATHOPPER:
			{
				bucket_bott->SetPos(ChVector<>(1, 0, 0));
				bucket_exist = false;
				break;
			}
			default:
			{
				break;
			}
			
		}
	}
}
std::shared_ptr<ChBody> SystemGeometry::create_bucketShell(int num_boxes, bool overlap)
{
	double t = bucket_half_thick;
	double r = bucket_rad;
	double h = bucket_interior_halfDim.z();
	double o_lap = 0;
	if (overlap) { o_lap = t * 2; }
	double th = h + o_lap;//for use in cyl volume;
	double cyl_volume = PPI*(2 * th - 2 * t)*(2 * th - 2 * t)*((2 * r + 2 * t)*(2 * r + 2 * t) - r*r) + (PPI)*(r + 2 * t)*(r + 2 * t) * 2 * t;
	double m = rho_cylinder*cyl_volume;

	return create_EmptyCylinder(num_boxes, overlap, true, h, t, r, bucket_ctr, true, bucketTexture, m);
}
std::shared_ptr<ChBody> SystemGeometry::create_topHook()
{
	auto tophk = std::make_shared<ChBody>();
	sys->AddBody(tophk);
	double mass = 0.001;
	tophk->SetIdentifier(topHookID);
	tophk->SetBodyFixed(false);
	tophk->SetCollide(true);
	double height = 5;
	if (bucketType == HOOKFRACTURE)
	{
		height = 20;
	}
	
	topOrigHeight = height*t2_smarticle;
	topHookVol = 2 * t2_smarticle*t2_smarticle*w_smarticle - t2_smarticle*t2_smarticle*t2_smarticle; //width of cross volume minus the overlap

	
	tophk->SetPos(ChVector<>(0,0,topOrigHeight));
	//ChVector<> 	gyr = utils::CalcBoxGyration(ChVector<>(w_smarticle, t2_smarticle, t2_smarticle)).Get_Diag() + utils::CalcBoxGyration(ChVector<>(t2_smarticle, w_smarticle, t2_smarticle)).Get_Diag();
	ChVector<> 	gyr = utils::CalcBoxGyration(ChVector<>(w_smarticle, t2_smarticle, t2_smarticle), tophk->GetPos()).Get_Diag() + utils::CalcBoxGyration(ChVector<>(t2_smarticle, w_smarticle, t2_smarticle), tophk->GetPos()).Get_Diag();
	tophk->SetMass(mass);
	tophk->SetInertiaXX(mass * gyr);

	tophk->SetMaterialSurface(mat_wall);
	tophk->GetCollisionModel()->ClearModel();
	tophk->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	//utils::AddBoxGeometry(tophk.get(), ChVector<>(w_smarticle,t2_smarticle,t2_smarticle), bucket_bott->GetPos()+ tophk->GetPos(), QUNIT, true);
	//utils::AddBoxGeometry(tophk.get(), ChVector<>(t2_smarticle, w_smarticle, t2_smarticle), bucket_bott->GetPos() + tophk->GetPos(), QUNIT, true);
	utils::AddBoxGeometry(tophk.get(), ChVector<>(w_smarticle, t2_smarticle, t2_smarticle), tophk->GetPos(), QUNIT, true);
	utils::AddBoxGeometry(tophk.get(), ChVector<>(t2_smarticle, w_smarticle, t2_smarticle), tophk->GetPos(), QUNIT, true);

	tophk->GetCollisionModel()->BuildModel();
	

	tophk->AddAsset(hookTexture);
	
	return tophk;
}
std::shared_ptr<ChBody> SystemGeometry::create_bottomHook()
{
	auto botthk = std::make_shared<ChBody>();
	sys->AddBody(botthk);
	double mass = 0.001;
	botthk->SetIdentifier(bottomHookID);
	botthk->SetBodyFixed(true);
	botthk->SetCollide(true);
	bottomHookVol = 2 * t2_smarticle*t2_smarticle*w_smarticle - t2_smarticle*t2_smarticle*t2_smarticle; //width of cross volume minus the overlap

	double height = 10;
	double bottOrigHeight = t2_smarticle*height;
	botthk->SetPos(ChVector<>(0, 0, height*t2_smarticle));
	//ChVector<> 	gyr = utils::CalcBoxGyration(ChVector<>(w_smarticle, t2_smarticle, t2_smarticle),botthk->GetPos()).Get_Diag() + utils::CalcBoxGyration(ChVector<>(t2_smarticle, w_smarticle, t2_smarticle)).Get_Diag();
	ChVector<> 	gyr = utils::CalcBoxGyration(ChVector<>(w_smarticle, t2_smarticle, t2_smarticle), botthk->GetPos()).Get_Diag() + utils::CalcBoxGyration(ChVector<>(t2_smarticle, w_smarticle, t2_smarticle), botthk->GetPos()).Get_Diag();
	botthk->SetMass(mass);
	botthk->SetInertiaXX(mass * gyr);

	botthk->SetMaterialSurface(mat_wall);
	botthk->GetCollisionModel()->ClearModel();
	botthk->GetCollisionModel()->SetEnvelope(collisionEnvelope);

	//utils::AddBoxGeometry(botthk.get(), ChVector<>(w_smarticle, t2_smarticle, t2_smarticle), bucket_bott->GetPos() + botthk->GetPos(), QUNIT, true);
	//utils::AddBoxGeometry(botthk.get(), ChVector<>(t2_smarticle, w_smarticle, t2_smarticle), bucket_bott->GetPos() + botthk->GetPos(), QUNIT, true);
	utils::AddBoxGeometry(botthk.get(), ChVector<>(w_smarticle, t2_smarticle, t2_smarticle), botthk->GetPos(), QUNIT, true);
	utils::AddBoxGeometry(botthk.get(), ChVector<>(t2_smarticle, w_smarticle, t2_smarticle), botthk->GetPos(), QUNIT, true);


	botthk->GetCollisionModel()->BuildModel();


	botthk->AddAsset(hookTexture);
	return botthk;
}



std::shared_ptr<ChBody> SystemGeometry::create_ChordRing(int num_boxes, double half_h, double t, double r, double sagitta, ChVector<> pos, std::shared_ptr<ChTexture> texture, double m)
{



	//https://en.wikipedia.org/wiki/Circular_segment
	auto chrdCirc = std::make_shared<ChBody>();
	chrdCirc->SetIdentifier(bucketID);
	//cyl_container->SetMass(mass);
	chrdCirc->SetPos(pos);
	chrdCirc->SetRot(QUNIT);
	chrdCirc->SetBodyFixed(false);
	chrdCirc->SetCollide(true);
	chrdCirc->SetName("ring");

	chrdCirc->SetMaterialSurface(mat_wall);
	double ang = 2.0 * PPI / num_boxes;
	ChQuaternion<> quat = QUNIT; //rotation of each plate
	double box_side = r * 2.0 * sin(PPI / num_boxes) / (cos(PPI / num_boxes));//side length of cyl
	double wallt = t / 5; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
	double o_lap = t * 2;
	ChVector<> pPos = (0, 0, 0);  //position of each plate
	double centAngle = 2 * acos(1 - sagitta / r);
	int roundPart = (2 * PPI - centAngle) / ang;
	ChVector<> box_size = ChVector<>((box_side + wallt) / 2.0, //size of plates
		wallt,
		half_h + o_lap);
	ChVector<> gyr(0, 0, 0);
	ChVector<> xxIner(0, 0, 0);
	ChVector<> cmRel(0, 0, 0);
	ChVector<> rel_loc(0, 0, 0);
	double mEach = m / (num_boxes);



	double clen = 2 * r*sin(centAngle / 2.0);//chord length
	double d = clen*clen / (8 * sagitta) - sagitta / 2;
	double boxChordAng = 2 * PPI - centAngle / 2;

	ChVector<> bigside_size(clen / 2.0,
		wallt,
		half_h + o_lap);


	double vSmall = box_size.x()*box_size.y()*box_size.z();
	double rho = mEach / vSmall;
	double vLarge = bigside_size.x()*bigside_size.y()*bigside_size.z();


	//#################################################################################
	//TODO change this to non-magic number!!!
	//#################################################################################
	double mExtra = 0;
	double mBig = m - mEach*(roundPart + 1) + mExtra;
	double mTot = m + mExtra;

	ChVector<>bigWallPos = pos + ChVector<>(sin(boxChordAng)*(d + box_side / 2),
		cos(boxChordAng)*(d + box_side / 2),
		0 - 1.8*t);

	//get cmREL
	for (int i = 0; i < roundPart + 1; i++)
	{
		pPos = pos + ChVector<>(sin(ang * i) * (wallt + r),
			cos(ang*i)*(wallt + r),
			0 - 1.8*t);
		cmRel = cmRel + (mEach*pPos);
	}
	cmRel = cmRel + mBig*bigside_size;
	cmRel = cmRel / mTot;

	for (int i = 0; i < roundPart + 1; i++)
	{

		pPos = pos + ChVector<>(sin(ang * i) * (wallt + r),
			cos(ang*i)*(wallt + r),
			0 - 1.8*t);

		quat = Angle_to_Quat(ANGLE, ChVector<>(0, 0, ang*i));
		rel_loc = pPos - cmRel;
		gyr = utils::CalcBoxGyration(box_size, pPos, quat).Get_Diag();
		xxIner.x() = xxIner.x() + mEach * (gyr.x() + ChVector<>(0, rel_loc.y(), rel_loc.z()).Length2());
		xxIner.y() = xxIner.y() + mEach * (gyr.y() + ChVector<>(rel_loc.x(), 0, rel_loc.z()).Length2());
		xxIner.z() = xxIner.z() + mEach * (gyr.z() + ChVector<>(rel_loc.x(), rel_loc.y(), 0).Length2());


		chrdCirc->AddAsset(texture);
		chrdCirc->GetCollisionModel()->SetEnvelope(collisionEnvelope);
		utils::AddBoxGeometry(chrdCirc.get(), box_size, pPos, quat, true);
	}
	///big side stuff

	quat = Angle_to_Quat(ANGLE, ChVector<>(0, 0, boxChordAng));
	rel_loc = bigWallPos - cmRel;
	gyr = utils::CalcBoxGyration(bigside_size, bigWallPos, quat).Get_Diag();

	xxIner.x() = xxIner.x() + mBig * (gyr.x() + ChVector<>(0, rel_loc.y(), rel_loc.z()).Length2());
	xxIner.y() = xxIner.y() + mBig * (gyr.y() + ChVector<>(rel_loc.x(), 0, rel_loc.z()).Length2());
	xxIner.z() = xxIner.z() + mBig * (gyr.z() + ChVector<>(rel_loc.x(), rel_loc.y(), 0).Length2());
	chrdCirc->AddAsset(texture);
	chrdCirc->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(chrdCirc.get(), bigside_size, bigWallPos, quat, true);

	chrdCirc->GetCollisionModel()->BuildModel();

	chrdCirc->SetMass(mTot);
	//chrdCirc->SetInertiaXX(xxIner);
	quat = Angle_to_Quat(ANGLE, ChVector<>(0, 0, 2 * PPI - boxChordAng - PPI / 2));
	//auto ro= chrdCirc->GetRot().Rotate(ChVector<>(0, 0, boxChordAng));
	chrdCirc->SetRot(quat);
	return chrdCirc;
}
std::shared_ptr<ChBody> SystemGeometry::create_EmptyEllipsev1(int num_boxes, bool overlap, bool createVector, double half_h, double t, double r, ChVector<> pos, bool halfVis, std::shared_ptr<ChTexture> texture, double m, double ax, double by)
{
	auto cyl_container = std::make_shared<ChBody>();
	cyl_container->SetIdentifier(bucketID);
	//cyl_container->SetMass(mass);
	cyl_container->SetPos(pos);
	cyl_container->SetRot(QUNIT);
	cyl_container->SetBodyFixed(false);
	cyl_container->SetCollide(true);
	ChVector<> box_size = (0, 0, 0); //size of plates
	//double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double wallt = t / 5; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
	//double half_height = bucket_interior_halfDim.z();
	double box_side = r * 2.0 * by*sin(PPI / num_boxes) / ax*(cos(PPI / num_boxes));//side length of cyl
	double o_lap = 0;
	if (overlap) { o_lap = t * 2; }
	double ang = 2.0 * PPI / num_boxes;
	ChVector<> pPos = (0, 0, 0);  //position of each plate
	ChQuaternion<> quat = QUNIT; //rotation of each plate


	double h = half_h + o_lap;
	//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
	//cyl_container->GetCollisionModel()->BuildModel();
	cyl_container->GetCollisionModel()->ClearModel();
	//https://en.wikipedia.org/wiki/List_of_moments_of_inertia
	double Ixx = m / 12.0*(3 * (r*r + (r - 2 * wallt)*(r - 2 * wallt)) + (h * 2)*(h + o_lap * 2));
	double Iyy = m / 12.0*(3 * (r*r + (r - 2 * wallt)*(r - 2 * wallt)) + (h * 2)*(h * 2));
	double Izz = m / 2 * (r*r + (r - 2 * wallt)*(r - 2 * wallt));
	ChMatrix33<double> iner(Ixx, 0.0, 0.0, 0.0, Iyy, 0.0, 0.0, 0.0, Izz);
	cyl_container->SetInertia(iner);



	cyl_container->SetMaterialSurface(mat_wall);
	cyl_container->SetName("ring");
	for (int i = 0; i < num_boxes; i++)
	{

		box_size = ChVector<>((box_side + wallt) / 2.0,
			wallt,
			half_h + o_lap);

		if (createVector)
		{
			//sin for x and cos for y because we want them *tangent* to angle specified!
			pPos = pos + ChVector<>(sin(ang * i) * (-wallt + r*ax),
				cos(ang*i)*(wallt + r*by),
				half_h);
		}
		else
		{
			pPos = pos + ChVector<>(sin(ang * i) * (-4 * wallt + ax*r),
				cos(ang*i)*(-4 * wallt + by*r),
				0 - 1.8*t);

			//TODO ######take into account angle of box!!
		}
		quat = Angle_to_Quat(ANGLE, ChVector<>(0, 0, ang*i));

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
		
		utils::AddBoxGeometry(cyl_container.get(), box_size, pPos, quat, m_visualization);
	}
	cyl_container->GetCollisionModel()->SetEnvelope(collisionEnvelope);

	cyl_container->SetMass(m);

	cyl_container->GetCollisionModel()->BuildModel();

	if (createVector)
	{
		for (int i = 0; i < num_boxes; i++)
		{
			auto wallPiece = std::make_shared<ChBody>();
			box_size = ChVector<>((box_side + wallt) / 2.0,
				wallt,
				half_h + o_lap);

			pPos = pos + ChVector<>(ax*sin(ang * i) * (-wallt + r),
				by*cos(ang*i)*(-wallt + r),
				half_h);

			quat = Angle_to_Quat(ANGLE, ChVector<>(0, 0, ang*i));

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
			wallPiece->SetMass(m);
			//wallPiece->SetPos(ChVector<>(0,0,0));


			//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
			wallPiece->GetCollisionModel()->BuildModel();
			wallPiece->GetCollisionModel()->SetFamily(envFamily); ////#############
			wallPiece->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
			sys->AddBody(wallPiece);
			wallPiece->SetMaterialSurface(mat_wall);
			wallPiece->SetRot(QUNIT);
			wallPiece->SetBodyFixed(true);
			wallPiece->SetCollide(true);
			bucket_bod_vec.emplace_back(wallPiece);

		}
	}
	GetLog() << cyl_container->GetMass();
	return cyl_container;
}

std::shared_ptr<ChBody> SystemGeometry::create_EmptyEllipse(int num_boxes, bool overlap, bool createVector, double half_h, double t, double r, ChVector<> pos, bool halfVis, std::shared_ptr<ChTexture> texture, double m, double ax, double by)
{
	auto cyl_container = std::make_shared<ChBody>();
	cyl_container->SetIdentifier(bucketID);
	//cyl_container->SetMass(mass);
	//cyl_container->SetPos(pos);
	//cyl_container->SetRot(QUNIT);
	cyl_container->SetBodyFixed(false);
	cyl_container->SetCollide(true);
	ChVector<> box_size = (0, 0, 0); //size of plates
																	 //double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double wallt = t / 5; //made this to disallow particles from sitting on thickness part of container, but keep same thickness for rest of system
												//double half_height = bucket_interior_halfDim.z(); (t/5= .0005)

	//if all other vars remain the same
	//abs(cos(pi)*(-wallt + r) * 2) - 2 * wallt represents inner diameter of the cylinder
	//if wallt was .0007 (currently .0005), inner diameter would be equal to 0.044 which is nicks exp.


	double box_side = r * 4.0 * by*sin(PPI / num_boxes) / ax*(cos(PPI / num_boxes));//side length of cyl
	double o_lap = 0;
	if (overlap) { o_lap = t * 2; }
	double ang = 2.0 * PPI / num_boxes;
	ChVector<> pPos = (0, 0, 0);  //position of each plate
	ChQuaternion<> quat = QUNIT; //rotation of each plate


	double h = half_h + o_lap;
	cyl_container->SetMaterialSurface(mat_wall);
	cyl_container->SetName("ring");
	cyl_container->GetCollisionModel()->ClearModel();
	cyl_container->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	cyl_container->SetMass(m);

	for (int i = 0; i < num_boxes; i++)
	{
		auto wallPiece = std::make_shared<ChBody>();
		box_size = ChVector<>((box_side + wallt) / 2.0,
			wallt,
			half_h + o_lap);

		pPos = ChVector<>(ax*sin(ang * i) * (-wallt + r),
			by*cos(ang*i)*(-wallt + r),
			half_h);

		quat = Angle_to_Quat(ANGLE, ChVector<>(0, 0, ang*i));

		//this is here to make half the cylinder invisible.

		//wallPiece->SetPos(pPos);
		wallPiece->SetRot(QUNIT);
		wallPiece->GetCollisionModel()->ClearModel();
		
		//wallPiece->GetCollisionModel()->AddBox(box_size.x(), box_size.y(), box_size.z(), pPos, ChMatrix33<>(quat));

		wallPiece->GetCollisionModel()->SetEnvelope(wallt);
		wallPiece->GetCollisionModel()->SetSafeMargin(wallt / 2);
		bool m_visualization = false;
		if (ang*i < 3 * PPI / 4 || ang*i > 5 * PPI / 4)
		{
			m_visualization = true;
			//auto box = std::make_shared<ChBoxShape>();
			//box->GetBoxGeometry().Size = box_size;
			//box->Pos = pPos;
			//box->Rot = ChMatrix33<>(quat);
			//wallPiece->GetAssets().push_back(box);
			//wallPiece->AddAsset(texture);
			
		}
		utils::AddBoxGeometry(cyl_container.get(), box_size, pPos, quat, m_visualization);
		wallPiece->SetMass(m);
		//wallPiece->SetPos(ChVector<>(0,0,0));


		//cyl_container->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelope);
		//wallPiece->GetCollisionModel()->BuildModel();
		wallPiece->GetCollisionModel()->SetFamily(envFamily); ////#############
		wallPiece->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	
		wallPiece->SetMaterialSurface(mat_wall);
		wallPiece->SetBodyFixed(true);
		wallPiece->SetCollide(true);
		bucket_bod_vec.emplace_back(wallPiece);
	}
	cyl_container->GetCollisionModel()->BuildModel();
	GetLog() << cyl_container->GetMass();
	sys->AddBody(cyl_container);
	return cyl_container;
}

std::shared_ptr<ChBody> SystemGeometry::create_EmptyCylinder(int num_boxes, bool overlap, bool createVector, double half_h, double t, double r, ChVector<> pos, bool halfVis, std::shared_ptr<ChTexture> texture, double m)
{
	return create_EmptyEllipse(num_boxes, overlap, createVector, half_h, t, r, pos, halfVis, texture, m, 1, 1);
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

	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hdim.x() + o_lap, 2.5*hdim.y() + o_lap, hthick), ChVector<>(0, 0, -hthick)); //bottom
	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hthick, hdim.y() + o_lap, hdim.z() + o_lap), //-x
		ChVector<>(-hdim.x() - hthick, 0, hdim.z()));
	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hthick, hdim.y() + o_lap, hdim.z() + o_lap), //+x
		ChVector<>(hdim.x() + hthick, 0, hdim.z()));
	//utils::AddBoxGeometry(flathopper.get(), ChVector<>(hdim.x() + o_lap, hthick, hdim.z() + o_lap), //-y
	//	ChVector<>(0, -hdim.y() - hthick, hdim.z()));
	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hdim.x() + o_lap, hthick, hdim.z() + o_lap), //+y
		ChVector<>(0, hdim.y() + hthick, hdim.z()));

	double hopperWallWidth = (hdim.x() - hole) / 2;
	double hopperWallHeight = hopperWallWidth / sin(wallAngle);
	double hopperWallAngleHeight = hopperWallHeight*cos(wallAngle);
	///////////hopper walls////////
	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hthick, hopperWallHeight, hdim.z() + o_lap), //left side
		ChVector<>(hdim.x() - hopperWallWidth, -hdim.y() - hopperWallHeight*cos(wallAngle), hdim.z()), Angle_to_Quat(ANGLE, ChVector<>(0, 0, wallAngle)));

	utils::AddBoxGeometry(flathopper.get(), ChVector<>(hthick, hopperWallHeight, hdim.z() + o_lap), //right side
		ChVector<>(-(hdim.x() - hopperWallWidth), -hdim.y() - hopperWallAngleHeight, hdim.z()), Angle_to_Quat(ANGLE, ChVector<>(0, 0, -wallAngle)));
	flathopper->GetCollisionModel()->BuildModel();

	flathopper->SetRot(Angle_to_Quat(ANGLE, ChVector<>(box_ang, 0, 0)));
	sys->AddBody(flathopper);
	ChVector<> buckBottPos = ChVector<>(bucket_ctr.x(), -hdim.y() - 2 * hopperWallAngleHeight + hthick, hdim.z());

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
	floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));//custom file
	bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddBoxGeometry(bucket_bott.get(), ChVector<>(hdim.x() + o_lap, hthick, hdim.z() + o_lap), //-y
		//ChVector<>(0, -hdim.y() - hthick, hdim.z()));
		buckBottPos);

	bucket_bott->AddAsset(floorTexture);
	bucket_bott->GetCollisionModel()->SetFamily(envFamily);
	bucket_bott->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	bucket_bott->GetCollisionModel()->BuildModel();
	bucket_bott->SetRot(Angle_to_Quat(ANGLE, ChVector<>(box_ang, 0, 0)));
	sys->AddBody(bucket_bott);

	return flathopper;
}

std::shared_ptr<ChBody> SystemGeometry::create_Hopper(double theta, bool overlap)
{
	auto hopper = std::make_shared<ChBody>();
	double t = bucket_half_thick; //bucket thickness redefined here for easier to read code
	double r = bucket_rad;
	double ang = theta;
	double h = bucket_interior_halfDim.z();
	double sH = (h - t) / sin(ang); //side height
	//double w = hole_size / 2 + cos(ang)*h;
	double w = cos(ang)*h - t / 2;
	hopper->SetPos(bucket_ctr);
	hopper->SetRot(QUNIT);
	hopper->SetBodyFixed(true);
	hopper->SetCollide(true);



	double o_lap = 0;
	if (overlap) { o_lap = 2 * t; }

	hopper->GetCollisionModel()->ClearModel();

	//bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_black_bordersBlack.png"));
	//hopper->AddAsset(bucketTexture);
	hopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);


	ChVector<> hop_pos = bucket_ctr;
	hopper->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	//bucket_rad + 2 * bucket_half_thick, bucket_rad + 2 * bucket_half_thick, bucket_half_thick
	//utils::AddBoxGeometry(ramp.get(), ChVector<>(w, w, t), rampPos, QUNIT, true);//bucket_bottom


	utils::AddBoxGeometry(hopper.get(), ChVector<>(t, r + t, h), hop_pos - VECT_X*(-sin(ang)*h + hole_size + cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLE, ChVector<>(0, -ang, 0)), true);// -x negAngle
	//utils::AddBoxGeometry(hopper.get(), ChVector<>(t, r + t, h), hop_pos - VECT_X*(sin(ang)*(h)+hole_size+cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLE, ChVector<>(0, ang, 0)), true);// -x negAngle +theta
	utils::AddBoxGeometry(hopper.get(), ChVector<>(t, r + t, h), hop_pos + VECT_X*(-sin(ang)*h + hole_size + cos(ang)*t) + VECT_Z*(h - o_lap), Angle_to_Quat(ANGLE, ChVector<>(0, ang, 0)), true);// +x  -theta = /


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
	bucket_bott->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	//original size
	//utils::AddBoxGeometry(bucket_bott.get(), Vector(bucket_rad + 2 * bucket_half_thick, bucket_rad + 2 * bucket_half_thick, bucket_half_thick), Vector(0, 0, -bucket_half_thick), QUNIT, true);
	utils::AddBoxGeometry(bucket_bott.get(), Vector(3*bucket_rad + 2 * bucket_half_thick, 3*bucket_rad + 2 * bucket_half_thick, bucket_half_thick), Vector(0, 0, -bucket_half_thick), QUNIT, true);
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
	double half_height = bucket_interior_halfDim.z() / (radMult * 2);
	double box_side = bucket_rad * radMult * 2 * tan(PPI / num_boxes);//side length of cyl
	double o_lap = 0;
	if (overlap) { o_lap = t * 2; }
	double ang = 2.0 * PPI / num_boxes;
	ChVector<> box_size = (0, 0, 0); //size of plates
	ChVector<> ridge_size = (0, 0, 0); //size of plates
	ChVector<> pPos = (0, 0, 0);  //position of each plate
	ChQuaternion<> quat = QUNIT; //rotation of each plate
	drum->GetCollisionModel()->ClearModel();

	int ridgeNum = num_boxes / ridges;
	for (int i = 0; i < num_boxes; i++)
	{

		box_size = ChVector<>((box_side + wallt) / 2.0,
			wallt,
			half_height + o_lap);
		pPos = bucket_ctr + ChVector<>(sin(ang * i) * (wallt + bucket_rad*radMult),
			cos(ang*i)*(wallt + bucket_rad*radMult),
			0);

		quat = Angle_to_Quat(ANGLE, ChVector<>(0, 0, ang*i));

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


	//utils::AddBoxGeometry(convexShape.get(), rampSize, rampPos, Angle_to_Quat(ANGLE, ChVector<>(box_ang, 0, 0)), true);



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

	case STRESSSTICK:
	{
		create_Bucket_Bott();
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
		bucket = create_bucketShell(25, true);
		
		double rodLen = bucket_interior_halfDim.z()*1.5;
		create_CentralColumn(rodLen);
		create_Truss();
		create_Prismatic(stick);

	}
	break;
	case HOOKRAISE:
	{
		create_Bucket_Bott();
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
		bucket = create_bucketShell(25, true);

		double rodLen = bucket_interior_halfDim.z()*1.5;
		create_CentralColumn(rodLen);
		create_Truss();
		//create_Prismatic(stick);
		setupBucketActuator(Q_from_AngAxis(PI_2, VECT_Y));
		
	}
	break;
	case HOOKRAISE2:
	{
	
		create_Bucket_Bott();
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
		bucket = create_bucketShell(25, true);

		//bottomHook = create_bottomHook();
		topHook = create_topHook();
		hookVol = /*bottomHookVol +*/ topHookVol;

		create_Truss();
		create_Prismatic(topHook);
		//setupBucketActuator(bucket->GetRot());

		
	}
	break;

	case HOOKFRACTURE:
	{

		create_Bucket_Bott();
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
		bucket = create_bucketShell(25, true);

		bottomHook = create_bottomHook();
		topHook = create_topHook();
		hookVol = bottomHookVol + topHookVol;

		create_Truss();
		create_Prismatic(topHook);
		//setupBucketActuator(bucket->GetRot());


	}
	break;
	case BOX:
	{
		//FUTNOTE uncomment for old box
		bucket = create_Box2();

		//Maze
		//bucket = create_BoxBig();
		//std::shared_ptr<ChBody> maze = create_Maze();
		bucket_bott = create_Bucket_Bott();
		bucket_bott->SetCollide(false);
		bucket_bott->SetPos(ChVector<>(5, 5, 5));
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_red_borderRed.png"));
		
		setupBucketActuator(Q_from_AngAxis(PI_2, VECT_Y));
		bucket->SetBodyFixed(true);
		break;
	}
	case BOXDROP:
	{
			create_Truss();
			create_VibrateLink(omega_bucket, actuation_amp, actuationStart, bucket_bott);
			break;

		break;
	}
	case CYLINDER:
	{
		create_Bucket_Bott();
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
		bucket = create_bucketShell(25, true);

		create_Truss();
		create_VibrateLink(omega_bucket, actuation_amp, actuationStart, bucket_bott);


		break;
	}
	case KNOBCYLINDER:
	{
		create_Bucket_Bott();
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
		bucket = create_bucketShell(25, true);
		create_Truss();
		create_VibrateLink(omega_bucket, actuation_amp, actuationStart, bucket_bott);

		unsigned int kpr = 4;//knobs per row
		unsigned int rows = 15; //knob per z
		double rodlen = bucket_interior_halfDim.z()*2.0;
		create_CentralColumn(rodlen);
		create_Truss();
		create_Knobs(kpr, rows, rodlen);
		setupBucketActuator(bucket->GetRot());
		break;
	}
	case FLATHOPPER:
	{
		bucket = create_FlatHopper(boxdim);
		floorTexture->SetTextureFilename(GetChronoDataFile("cubetexture_brown_bordersBlack.png"));
		GetLog() << "\n\nfriction: " << bucket->GetMaterialSurfaceNSC()->GetKfriction() << " " << bucket->GetMaterialSurfaceNSC()->GetSfriction() << "\n";
		break;
	}
	case HOPPER:
	{
		bucket = create_Hopper(box_ang, true);
		bucketTexture->SetTextureFilename(GetChronoDataFile("cubetexture_black_bordersBlack.png"));
		create_Bucket_Bott();

		create_Truss();
		create_VibrateLink(omega_bucket, actuation_amp, actuationStart, bucket);
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

		setupBucketActuator(bucket->GetRot());
		bucket->SetBodyFixed(true);
		break;
	}
	}
	bucket->SetMaterialSurface(mat_wall);
	bucket_bott->SetMaterialSurface(mat_wall);
	bucket->AddAsset(bucketTexture);
	bucket->SetBodyFixed(true);
	bucket->SetCollide(true);
	bucket->GetCollisionModel()->SetFamily(envFamily);
	bucket->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

}

void SystemGeometry::create_Ground()
{

	ChVector<> boxDim = sizeScale * ChVector<>(0.1, 0.1, .002);
	ChVector<> boxLoc = sizeScale * ChVector<>(0, 0, -5.0*this->bucket_interior_halfDim.z());

	ground = std::make_shared<ChBody>();

	ground->SetMaterialSurface(mat_wall);
	ground->SetPos(boxLoc);

	// ground->SetIdentifier(-1);
	ground->SetBodyFixed(true);
	ground->SetCollide(true);

	ground->GetCollisionModel()->ClearModel();
	ground->GetCollisionModel()->SetEnvelope(collisionEnvelope);
	utils::AddCylinderGeometry(ground.get(), boxDim.x(), boxDim.z(), ChVector<>(0, 0, 0), Q_from_AngAxis(PI_2, VECT_X));
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
	vibrate_link = std::make_shared<ChLinkLockLock>();
	vibrate_link->Initialize(truss, body, ChCoordsys<>(ChVector<>(0, 0, 0)));
	double phase = -w*t_0;
	//auto vibMotion = std::make_shared<ChFunction_Sine>();


	//ChFunction_Sine* vibMotion = new ChFunction_Sine();  // phase freq ampl
	auto vibMotion = std::make_shared<ChFunction_Sine>();  // phase freq ampl

	vibMotion->Set_phase(phase);
	vibMotion->Set_amp(A);
	vibMotion->Set_w(w);
	vibrate_link->SetMotion_Z(vibMotion);
	vibrate_link->SetDisabled(true);
	sys->Add(vibrate_link);

}
void SystemGeometry::create_Prismatic(std::shared_ptr<ChBody> body)
{

	////////////////

	////
	//ChFunctionCustom* pris_motion = new ChFunctionCustom(0, 1.5, 0);  // phase freq ampl

//	
	//auto pris_motion = std::make_shared<ChFunctionCustom>(0, actuationSpd+.1, 0);  // phase freq ampl
	//pris_motion->Set_y(.1);
	//pris_motion->Set_y_dx(.1);
	//pris_motion->Set_y_dxdx(.1);
	////

	auto func = std::make_shared<ChFunctionCustom>();
	pris_link = std::make_shared<ChLinkLockLock>();


	link_prismatic = std::make_shared<ChLinkLockPrismatic>();
	pris_engine = std::make_shared<ChLinkLinActuator>();
	double offset = 0.5;//need offset for operation of prismatic in both directions
	

	//auto sinefunc = std::make_shared<ChFunction_Sine>();
	if (bucketType == HOOKRAISE2)
	{
		link_prismatic->Initialize(body, truss, 
			ChCoordsys<>(body->GetPos(), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Z)));
		pris_engine->Initialize(body, truss,false,
			ChCoordsys<>(body->GetPos() + ChVector<>(0,0, offset), QUNIT), //need to offset body marker by offset or it will move when sim starts
			ChCoordsys<>(body->GetPos(), QUNIT));
		pris_engine->Set_lin_offset(offset);
		//pris_engine->Initialize(body, truss, false, ChCoordsys<>(truss->GetPos() + (0, 0, .2), QUNIT), ChCoordsys<>(body->GetPos(), QUNIT));
		//pris_engine->Initialize(truss, body, ChCoordsys<>(VNULL, QUNIT));
		//pris_engine->GetMask()->Constr_N(3).SetMode(CONSTRAINT_LOCK);

		//for (int i = 0; i < 7; i++)
		//	pris_engine->GetMask()->Constr_N(i).SetMode(CONSTRAINT_LOCK);
		
		//pris_engine->GetMask()->Constr_N(2).SetMode(CONSTRAINT_FREE); //set z to free

	}
	if (bucketType == HOOKFRACTURE)
	{
		link_prismatic->Initialize(body, truss,
			ChCoordsys<>(body->GetPos(), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Z)));
		pris_engine->Initialize(body, truss, false,
			ChCoordsys<>(body->GetPos() + ChVector<>(0, 0, offset), QUNIT), //need to offset body marker by offset or it will move when sim starts
			ChCoordsys<>(body->GetPos(), QUNIT));
		pris_engine->Set_lin_offset(offset);
		//pris_engine->Initialize(body, truss, false, ChCoordsys<>(truss->GetPos() + (0, 0, .2), QUNIT), ChCoordsys<>(body->GetPos(), QUNIT));
		//pris_engine->Initialize(truss, body, ChCoordsys<>(VNULL, QUNIT));
		//pris_engine->GetMask()->Constr_N(3).SetMode(CONSTRAINT_LOCK);

		//for (int i = 0; i < 7; i++)
		//	pris_engine->GetMask()->Constr_N(i).SetMode(CONSTRAINT_LOCK);

		//pris_engine->GetMask()->Constr_N(2).SetMode(CONSTRAINT_FREE); //set z to free

	}
	if (bucketType == STRESSSTICK || bucketType == HOOKRAISE)
	{
		//auto pris_motion = std::make_shared<ChFunctionCustom>(0, 1.5, .1);  // phase freq ampl
		//pris_link->Initialize(body, truss, ChCoordsys<>(ChVector<>(0, 0, 0)));
		//pris_link->SetMotion_Z(pris_motion);

		link_prismatic->Initialize(truss, body, ChCoordsys<>(ChVector<>(0, 0, 0)));
		//link_prismatic->SetMotion_Z(pris_motion);

		pris_engine->Initialize(truss, body, false, ChCoordsys<>(body->GetPos() + ChVector<>(0, 0, -stickLen), QUNIT), ChCoordsys<>(body->GetPos() + ChVector<>(0, 0, stickLen), QUNIT));
		func->Set_y(0);
		func->Set_y_dx(2.5 - .5); //the value in this is always -2.5+(value specified), dont know where -2.5 comes from....
		//pris_engine->Set_dist_funct(func);

	}


	sys->AddLink(link_prismatic);
	sys->AddLink(pris_engine);
	
	
	//pris_engine->SetDisabled(false);
	//link_prismatic->SetDisabled(false);
	//pris_motion = std::dynamic_pointer_cast<ChFunctionCustom>(pris_engine->Get_dist_funct());
	
}

void SystemGeometry::create_Truss()
{

	truss = std::make_shared<ChBody>();
	sys->AddBody(truss);
	truss->SetPos((0, 0, 0));
	//truss->SetMass(100);
	//truss->SetInertiaXX((30, 30, 30));
	truss->SetBodyFixed(true);
	truss->SetCollide(false);
	//truss->GetCollisionModel()->ClearModel();
	//utils::AddCylinderGeometry(truss.get(), t2_smarticle / 2, bucket_interior_halfDim.z() * 1, bucket_ctr + ChVector<>(0, 0, bucket_interior_halfDim.z()), Angle_to_Quat(ANGLE, ChVector<>(-PI_2, 0, 0)), true);
	//truss->GetCollisionModel()->BuildModel();
	//truss->AddAsset(sphereTexture);
	

	
}

void SystemGeometry::create_Knobs(double kpr, double rows, double length)
{
	//utils should be -PI_2?
	//utils::AddCylinderGeometry(truss.get(), t2_smarticle / 2, bucket_interior_halfDim.z() * 1, bucket_ctr + ChVector<>(0, 0, bucket_interior_halfDim.z()), Angle_to_Quat(ANGLE, ChVector<>(PI_2, 0, 0)), true);
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
			utils::AddBoxGeometry(stick.get(), ChVector<>(knobRad*1.5, rad / 4, rad / 8), bucket_ctr + ChVector<>(rad*cos(theta), rad*sin(theta), hp*(row + 1)), Angle_to_Quat(ANGLE, ChVector<>(0, 0, theta + row % 2 * PI_2)), true);
			sphereStick.emplace_back(stick);
		}
	}
	stick->SetMaterialSurface(mat_wall);
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

	columnEngine = std::make_shared<ChLinkEngine>();
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

	//bucket_interior_halfDim.z()*1.5;
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
		utils::AddSphereGeometry(stick.get(), rad, bucket_ctr + ChVector<>(0, 0, stickLen / sphereNum * (i)), Angle_to_Quat(ANGLE, ChVector<double>(0, 0, PPI)), true);
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
				utils::AddSphereGeometry(stick.get(), rad, bucket_ctr + ChVector<>(rad*(i + 1), 0, stickLen / sphereNum), Angle_to_Quat(ANGLE, ChVector<double>(PPI, 0, 0)), true);

			}
			else
			{
				utils::AddSphereGeometry(stick.get(), rad, bucket_ctr + ChVector<>(rad*(i + 1), 0, stickLen / sphereNum), Angle_to_Quat(ANGLE, ChVector<double>(PPI, 0, 0)), true);
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
}
void SystemGeometry::vibrate_body(double t, double w, double A, double t_0, std::shared_ptr<ChBody> body)  //!!!!!!!!!!!USE CREATE_VIBRATELINK INSTEAD!!!!!!!!!!!!!!
{
	//A=amplitude
	//w=omega
	//t_0 vibration starting time
	double phase = -w*t_0;
	double x_bucket = A*sin(w * t + phase);
	double xDot_bucket = A*w*cos(w * t + phase);
	double xDDot_bucket = -1 * A*w*w*sin(w * t + phase);
	body->SetPos(body->GetPos() + ChVector<>(0, 0, x_bucket));
	body->SetPos_dt(ChVector<>(0, 0, xDot_bucket));
	body->SetPos_dtdt(ChVector<>(0, 0, xDDot_bucket));
	body->SetRot(QUNIT);
}

void SystemGeometry::setupBucketActuator(ChQuaternion<double> rot)
{

	std::shared_ptr<ChFunction_Const> mfun2; //needs to be declared outside switch
	ChVector<> pR01(0, 0, 0);
	ChQuaternion<> qx;
	bucket_actuator = std::make_shared<ChLinkEngine>();
	
	switch (bucketType)
	{
		case DRUM:
		{

			bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), rot));
			bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
			//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
			sys->AddLink(bucket_actuator);
			mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_spe_funct());
			mfun2->Set_yconst(0);
			//mfun2->Set_yconst(drum_omega);
			break;
		}
		case BOX:
		{
			/*qx = Q_from_AngAxis(PI_2, VECT_Y);*/
			bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), rot));
			bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
			//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
			sys->AddLink(bucket_actuator);
			mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
			mfun2->Set_yconst(0);
			break;

		}
		case FLATHOPPER:
		{
			//qx = Q_from_AngAxis(PI_2, VECT_Y);
			bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), rot));
			bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
			//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
			//sys->AddLink(bucket_actuator);
			mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
			mfun2->Set_yconst(0);
			break;

		}
		case KNOBCYLINDER:
		{
			
			bucket_actuator->Initialize(stick, truss, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), rot));
			bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
			//drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
			//sys->AddLink(bucket_actuator);
			mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_spe_funct());
			mfun2->Set_yconst(0);
			//mfun2->Set_yconst(drum_omega);
		}
		case HOOKRAISE:
		{
			bucket_actuator->Initialize(bucket_bott, bucket, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), rot));
			bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_ROTATION);
			bucket_actuator->SetMotion_axis(ChVector<>(0, 0, 1));
			//sys->AddLink(bucket_actuator);
			mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
			mfun2->Set_yconst(0);
			break;
		}
		case HOOKRAISE2:
		{

			//auto func = std::make_shared<ChFunctionCustom>();
			//pris_link = std::make_shared<ChLinkLockLock>();
	

			//auto prismatic1 = std::make_shared<ChLinkLockPrismatic>();
			//prismatic1->Initialize(ground, slider1, ChCoordsys<>(ChVector<>(0, 0, -1), Q_from_AngY(CH_C_PI_2)));
			//prismatic1->GetLimit_Z()->Set_active(true);
			//prismatic1->GetLimit_Z()->Set_min(-6);
			//system.AddLink(prismatic1);

			// .. create the prismatic joint between the fork and arm
			//bucket_actuator = std::make_shared<ChLinkLinActuator>();
			//link_prismaticFork->Initialize(
			//	fork, arm,
			//	ChCoordsys<>(
			//		POS_prismatic,
			//		chrono::Q_from_AngAxis(CH_C_PI / 2,
			//			VECT_X)));  // set prism as vertical (default would be aligned to z, horizontal
			//app->GetSystem()->AddLink(link_prismaticFork);

			//// .. create the linear actuator that pushes upward the fork
			//link_actuatorFork = std::make_shared<ChLinkLinActuator>();
			//link_actuatorFork->Initialize(fork, arm, false,
			//	ChCoordsys<>(POS_prismatic + ChVector<>(0, 0.01, 0), QUNIT),
			//	ChCoordsys<>(POS_prismatic, QUNIT));
			//app->GetSystem()->AddLink(link_actuatorFork);

			bucket_actuator->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_PRISM);
			bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
			
			//link_prismatic = std::make_shared<ChLinkLockPointLine>();
			//bucket_actuator->Initialize(topHook, truss, ChCoordsys<>(ChVector<>(0, 0, 0)));
			//pris_link->SetMotion_Z(pris_motion);

			//bucket_actuator->Initialize(topHook, truss, true, ChCoordsys<>(topHook->GetPos() + ChVector<>(0, 0, -stickLen), QUNIT), ChCoordsys<>(truss->GetPos(), QUNIT));
			bucket_actuator->Initialize(topHook, truss, true, ChCoordsys<>(topHook->GetPos(), QUNIT), ChCoordsys<>(truss->GetPos(), QUNIT));
			mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_spe_funct());
			mfun2->Set_yconst(3);
			bucket_actuator->SetMotion_Z(mfun2);
			//func->Set_y(0);
			//func->Set_y_dx(2.5 - .5); //the value in this is always -2.5+(value specified), dont know where -2.5 comes from....
			//pris_engine->Set_dist_funct(func);
			//	break;

			////qx = Q_from_AngAxis(PI_2, VECT_Y);
			//bucket_actuator->Initialize(topHook, truss, ChCoordsys<>(bucket->GetRot().Rotate(pR01) + bucket->GetPos(), rot));
			//bucket_actuator->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
			////drum_actuator->SetMotion_axis(ChVector<>(1, 0, 0));
			////sys->AddLink(bucket_actuator);
			//mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_spe_funct());
			//mfun2->Set_yconst(3);

			break;
		}
	}
	sys->AddLink(bucket_actuator);
	bucket_actuator->SetDisabled(true);
}

void SystemGeometry::rotate_body_rot(double t, std::shared_ptr<ChBody> body, std::shared_ptr<ChLinkEngine> actuator, double ang)
{
	//static std::shared_ptr<ChFunction_Const> mfun2;
	std::shared_ptr<ChFunction_Const> mfun2;
	body->SetBodyFixed(false);
	mfun2 = std::dynamic_pointer_cast<ChFunction_Const>(bucket_actuator->Get_rot_funct());
	//set rotation amount equal to box_ang-currentangle, therefore no rotation if 2 values are equal
	mfun2->Set_yconst(box_ang - ang);
	//bucket->GetRot()).x()
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
