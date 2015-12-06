#ifndef INCLUDE_CONTROLLER_H_
#define INCLUDE_CONTROLLER_H_

#include <queue>
#include <vector>
#include <motion_functions/ChFunction_Base.h>
#include <physics/ChSystem.h>
#include <physics/ChLinkEngine.h>
#include <utility>      // std::pair, std::make_pair

namespace chrono {

	//class ChReportContactCallback2;
	class Smarticle;

	class Controller {
	public:
		//Controller(chrono::ChSystem *ch_system, Smarticle *smarticle) : ch_system_(ch_system), smarticle_(smarticle){}
		Controller(chrono::ChSystem *ch_system, Smarticle *smarticle);
		~Controller();
		// Step the controller
		bool Step(double dt);
		// get the toruqe for joint i
		//size_t GetNumEngines();
		ChSharedPtr<ChLinkEngine> GetEngine(size_t i);
		
		void SetDesiredAngle(size_t i, double desangle);
		double GetCurrAngle(size_t i, double t);
		double GetDesiredAngle(size_t i, double t);

		double GetCurrAngularSpeed(size_t i, double t);
		double GetDesiredAngularSpeed(size_t i, double t);
		double GetDesiredAngularSpeed2(size_t i, double t, double error);
		double GetDesiredAngularSpeedForFunction(size_t i, double t);
		double GetCurrTorque(size_t i, double t);
		
		double GetMediaTorque(size_t i, double t);
		double GetContactTorque(size_t i, double t);
		double GetAngle(size_t i, double t);
		double GetAngularSpeed(size_t i, double t);
		bool OT();
		void UseForceControl();
		void UseSpeedControl();
		// Change the undultation amplitude of the snake
		//void SetCommandAmplitude(double amp);
		//void PushCommandToQueue(const Json::Value &command);
		//std::vector<std::pair<double, double>>* getMoveVector(unsigned int guiState);

		double outputLimit = 0;
		//double p_gain = 0.0000350;
		//double i_gain = 0.0000300;
		//double d_gain = 0.000005;
		double p_gain = 0;
		double i_gain = 0;
		double d_gain = 0;
		void setMoveVector(unsigned int guiState);
	protected:
		// Process the commands in the queue, one at a time
		void ProcessCommandQueue(double dt);
		chrono::ChSystem *ch_system_;
		Smarticle *smarticle_;
		// Contact force on each of the robot segment.
		//chrono::ChReportContactCallback2 *contact_reporter_;
		//std::vector<chrono::ChVector<> > contact_force_list_;
		// the torques at joints, computed from contact forces.
		// Parametr for the CPG

		std::vector <double> desiredOmega_;
		std::vector <double> desiredAngle_;
		std::vector <double> currAngle_;
		std::vector <double> currOmega_;
		std::vector <double> currTorque_;


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
