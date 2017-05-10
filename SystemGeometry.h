#ifndef SYSTEMGEOMETRY_H_
#define SYSTEMGEOMETRY_H_

#include <stdlib.h>
#include "common.h"

#include "Smarticle.h"
#include <irrlicht.h>
#include "chrono_irrlicht/ChBodySceneNode.h"  //changed path from unit to chrono to reflect changes in updated chrono
#include "chrono_irrlicht/ChBodySceneNodeTools.h"
#include "chrono_irrlicht/ChIrrWizard.h"
#include "core/ChRealtimeStep.h"
#include "assets/ChTexture.h"
#include "chrono/core/ChFileutils.h"
#include <iostream>
#include "assets/ChTexture.h"
#include "chrono_irrlicht/ChIrrApp.h"//changed path from unit to chrono to reflect changes in updated chrono
#include "chrono_irrlicht/ChIrrTools.h"

namespace chrono {

			class SystemGeometry{

			public:
				
				SystemGeometry(std::shared_ptr <CH_SYSTEM> msys, BucketType sysType,double collisionEnv,double l_smart, double w_smart, double t_smart, double t2_smart);
				~SystemGeometry(); 

				
				//methods
				BucketType getGeom();
				void create_Ground();
				void create_Container();
				std::shared_ptr<ChBody> create_Maze();
				std::shared_ptr<ChBody> create_Box();
				std::shared_ptr<ChBody> create_BoxBig();
				std::shared_ptr<ChBody> create_EmptyEllipse(int num_boxes, bool overlap, bool createVector, double half_height, double thick, double rad, ChVector<> pos, bool halfVis, std::shared_ptr<ChTexture> texture, double mass, double ax, double by);
				std::shared_ptr<ChBody> create_EmptyCylinder(int num_boxes, bool overlap, bool createVector, double half_height, double thick, double rad, ChVector<> pos, bool halfVis, std::shared_ptr<ChTexture> texture, double mass);
				std::shared_ptr<ChBody> create_bucketShell(int num_boxes, bool overlap);
				std::shared_ptr<ChBody> create_FlatHopper(ChVector<> hdim);
				std::shared_ptr<ChBody> create_Hopper(double theta, bool overlap);
				std::shared_ptr<ChBody> create_Bucket_Bott();
				std::shared_ptr<ChBody> create_Drum(int num_boxes, bool overlap, int ridges = 5);
				std::shared_ptr<ChBody> create_Hull(double numBoxes);
				std::shared_ptr<ChBody> create_ChordRing(int num_boxes, double half_h, double t, double r, double sagitta, ChVector<> pos, std::shared_ptr<ChTexture> texture, double m);
				void										create_CentralColumn(double length);
				void										create_Prismatic(std::shared_ptr<ChBody> body);
				void										create_Knobs(double kpr, double rows, double length);
				void										create_CentralColumnEngine(double t_0);
				void										create_Truss();
				void										create_VibrateLink(double w, double A, double t_0, std::shared_ptr<ChBody> body);

				//void vibrate_body(double t, std::shared_ptr<ChBody> mainbody, std::shared_ptr<ChBody> body);
				void setUpBucketActuator(ChQuaternion<double> rot);
				void setUpBucketActuator();
				void vibrate_body(double t, double w, double A, double t_0, std::shared_ptr<ChBody> body);
				void rotate_body_rot(double t, std::shared_ptr<ChBody> body, std::shared_ptr<ChLinkEngine> actuator,double ang);
				void rotate_body_sp(double t, std::shared_ptr<ChBody> body, std::shared_ptr<ChLinkEngine> actuator, double w);
				//vars
				std::shared_ptr<CH_SYSTEM> sys;

				
				std::shared_ptr<ChBody> bucket_bott;
				std::shared_ptr<ChBody> ground;
				std::shared_ptr<ChBody> bucket;
				std::shared_ptr<ChBody> truss;
				std::shared_ptr<ChBody> stick;
				

				std::vector<std::shared_ptr<ChBody>> sphereStick;
				
				std::shared_ptr<ChLinkLock> vibrate_link;
				std::shared_ptr<ChLinkLock> pris_link;
				std::shared_ptr<ChLinkEngine> columnEngine;
				std::shared_ptr<ChLinkEngine> bucket_actuator;
				std::shared_ptr<ChLinkLinActuator> pris_engine;
				std::shared_ptr<ChLinkLockPrismatic> link_prismatic;


				static std::shared_ptr <ChTexture> bucketTexture;
				static std::shared_ptr <ChTexture> sphereTexture;
				static std::shared_ptr <ChTexture> floorTexture;
				static std::shared_ptr <ChTexture> groundTexture;
				static std::shared_ptr<SOLVER(ChMaterialSurface)> mat_wall;

				double collisionEnvelope;

				double bucket_rad;
				double bucket_half_thick;
				double stickLen;
				int envFamily;
				double rad;
				double w_smarticle;
				double l_smarticle;
				double t_smarticle;
				double t2_smarticle;
	
				double rho_cylinder;
				double hole_size;

				ChVector<> bucket_interior_halfDim;
				ChVector<> boxdim;
				ChVector<> bucket_ctr;
			private:

				//vars
				int bucketID;
				double wall_fric;//.3814; //keyboard box friction = .3814



				//system type
				BucketType sType;
			};
}
	
#endif /* SMARTICLE_H_ */
