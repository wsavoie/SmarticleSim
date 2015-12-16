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
		ChSharedPtr<ChLinkEngine> GetEngine(size_t i);
		ChSharedPtr<ChFunctionController> engine_funct0;
		ChSharedPtr<ChFunctionController> engine_funct1;
		void SetDesiredAngle(size_t i, double desOmeg);
		void SetDesiredAngularSpeed(size_t i, double desangle);
		//double GetCurrAngle(size_t i, double t);
		//double SetGetCurrAngle(size_t i, double t);
		double GetDesiredAngle(size_t i, double t);
		//void SetCurrAngle(size_t i, double ang);
		double GetActuatorOmega(size_t i, double t);
		double GetCurrOmega(size_t i, double t);
		double GetDesiredAngularSpeed(size_t i, double t);
		double GetDesiredAngularSpeed2(size_t i, double t, double error);
		double GetDesiredAngularSpeedForFunction(size_t i, double t);
		double GetCurrTorque(size_t i, double t);
		double LinearInterpolate(size_t i, double curr, double desired);
		double OmegaLimiter(size_t i, double omega);
		double GetAngle(size_t i, double t);
		double GetAngularSpeed(size_t i, double t);
		void CalcCurr_Omega(size_t i, double t);
		bool OT();
		void UseForceControl();
		void UseSpeedControl();

		double omegaLimit = 5;
		double outputLimit = 0;

		bool resetCumError = false;
		//void setMoveVector(unsigned int guiState);
		std::vector <double> prevError_;
		std::vector <double> cumError_;
		std::vector <double> prevAngle_;
		std::vector <bool> successfulMove_;
	protected:
		// Process the commands in the queue, one at a time
		void ProcessCommandQueue(double dt);
		chrono::ChSystem *ch_system_;
		Smarticle *smarticle_;
		//ChSharedPtr<ChFunctionController> engine_funct0;
		//ChSharedPtr<ChFunctionController> engine_funct1;
		// Contact force on each of the robot segment.
		//chrono::ChReportContactCallback2 *contact_reporter_;
		//std::vector<chrono::ChVector<> > contact_force_list_;

		//std::vector <double> desiredOmega_;
		//std::vector <double> desiredAngle_;
		//std::vector <double> currAngle_;
		//std::vector <double> currOmega_;
		//std::vector <double> currTorque_;
//		std::vector <double> prevError_;
		double num_waves_ = 2.0;
		double default_amplitude_ = 0.1;
		double command_amplitude_ = default_amplitude_;
		int command_count_down_ = 0;
		// reference counting
		size_t steps_ = 0;
		double friction = 1;

	};
}
#endif // INCLUDE_CHFUNCTION_CONTROLLER_H_
