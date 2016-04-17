#ifndef INCLUDE_CONTROLLER_H_
#define INCLUDE_CONTROLLER_H_

#include <queue>
#include <vector>
#include <motion_functions/ChFunction_Base.h>
#include <physics/ChSystem.h>
#include <physics/ChLinkEngine.h>
#include "ChFunction_Controller.h"


namespace chrono {


	//class ChReportContactCallback2;
	class Smarticle;

	class Controller {
	public:
		Controller(chrono::ChSystem *ch_system, Smarticle *smarticle);
		~Controller();
		// Step the controller
		bool Step(double dt);
		// get the toruqe for joint i
		//size_t GetNumEngines();

		std::shared_ptr<ChLinkEngine> GetEngine(size_t i);
		std::shared_ptr<ChFunctionController> engine_funct0;
		std::shared_ptr<ChFunctionController> engine_funct1;
		
		double GetDesiredAngle(size_t i, double t);
		double GetExpAngle(size_t idx, double t);
		double GetActuatorOmega(size_t i, double t);
		double GetCurrOmega(size_t i, double t);
		double GetCurrTorque(size_t i, double t);
		double LinearInterpolate(size_t i, double curr, double desired);
		double OmegaLimiter(size_t i, double omega);
		double GetAngle(size_t i, double t);
		double GetAngularSpeed(size_t i, double t);
		void CalcCurr_Omega(size_t i, double t);
		bool OT();
		void UseForceControl(size_t i);
		void UseSpeedControl(size_t i);
		double velOld[2];
		double DD[2];
		double II[2];
		double yold[2];
		double omegaLimit = 5;
		double outputLimit = 0;
		Smarticle *smarticle_;
		bool resetCumError = false;
		//void setMoveVector(unsigned int guiState);
		std::vector <double> prevError_;
		std::vector <double> prevOmegError_;
		std::vector <double> cumError_;
		std::vector <double> cumOmegError_;
		std::vector <double> prevAngle_;
		std::vector <bool> successfulMove_;
	protected:
		chrono::ChSystem *ch_system_;

		double num_waves_ = 2.0;
		double default_amplitude_ = 0.1;
		double command_amplitude_ = default_amplitude_;
		int command_count_down_ = 0;
		// reference counting
		
		double friction = 1;

	};
}
#endif // INCLUDE_CHFUNCTION_CONTROLLER_H_
