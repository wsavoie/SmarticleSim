#ifndef INCLUDE_CONTROLLER_H_
#define INCLUDE_CONTROLLER_H_

#include <queue>
#include <vector>
#include <motion_functions/ChFunction_Base.h>
#include <physics/ChLinkEngine.h>
#include "ChFunction_Controller.h"
#include "common.h"
//#include "chrono/physics/ChSystemSMC.h"
//#include "chrono/physics/ChSystemNSC.h"
namespace chrono {


	//class ChReportContactCallback2;
	class Smarticle;

	class Controller {
	public:
		Controller(std::shared_ptr<CH_SYSTEM> ch_system, Smarticle* smarticle);
		~Controller();
		// Step the controller
		bool Step(double dt);
		// get the toruqe for joint i
		//size_t GetNumEngines();

		std::shared_ptr<ChLinkEngine> GetEngine(size_t i);
		std::shared_ptr<ChFunctionController> engine_funct0;
		std::shared_ptr<ChFunctionController> engine_funct1;

		double smartStr;

		void setSmartStr();
		double GetDesiredAngle(size_t i, double t);
		double GetExpAngle(size_t idx, double t);
		double GetActuatorOmega(size_t i, double t);
		double GetCurrOmega(size_t i, double t);
		double GetCurrTorque(size_t i, double t);
		double GetCurrReactTorque(size_t i, double t);
		double LinearInterpolate(size_t i, double curr, double desired);
		double OmegaLimiter(size_t i, double omega);
		double GetAngle(size_t i, double t);
		double GetAngularSpeed(size_t i, double t);
		void CalcCurr_Omega(size_t i, double t);
		bool OT();
		void UseForceControl(size_t i);
		void UseSpeedControl(size_t i);
		void UsePositionControl(size_t i);

		double torOld[2];
		double torCur[2];

		double velOld[2];
		double velCur[2];

		double posOld[2];
		double posCur[2];

		double prevAngle[2];//previous directed angle used to check trajectory
		double desPrev[2];

		double DD[2];
		double II[2];
		double yold[2];
		double ycurr[2];


		double omegaLimit = 5;
		double outputLimit = 0;
		std::shared_ptr<Smarticle> smarticle_;
		bool resetCumError = false;
		//void setMoveVector(unsigned int guiState);
		std::vector <double> prevError_;
		int maxDAvgSize = 100;
		std::deque <double> dAvg0_;
		std::deque <double> dAvg1_;
		std::vector <double> prevOmegError_;
		std::vector <double> cumError_;
		std::vector <double> cumOmegError_;
		std::vector <double> prevAngle_;
		std::vector <bool> successfulMove_;

		// last input
		std::vector <double>  mLastError;
		// accumulated input
		double mAccuError[2];
		// last time the function is called
		std::vector <double>  mLastCalled;
		// last output
		std::vector <double> mLastValue;
		std::shared_ptr<CH_SYSTEM> ch_system_;
	//protected:
		


	};
}
#endif // INCLUDE_CHFUNCTION_CONTROLLER_H_
