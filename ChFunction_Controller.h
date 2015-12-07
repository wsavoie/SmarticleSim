#ifndef INCLUDE_CHFUNCTION_CONTROLLER_H_
#define INCLUDE_CHFUNCTION_CONTROLLER_H_

#include <cmath>
#include <algorithm>

#include <motion_functions/ChFunction_Base.h>


namespace chrono {
	//class Smarticle;
	class Controller;

	class ChFunctionController : public ChFunction {
	public:

		//ChFunctionController(size_t index, Smarticle* smarticle)
		//	: index_(index), smarticle_(smarticle) {}
		ChFunctionController(size_t index, Controller* controller)
			: index_(index), controller_(controller) {}
		virtual ~ChFunctionController(){};
		ChFunction *new_Duplicate() {
			return new ChFunctionController(index_, controller_);

		}

		
		void ResetCumulative();
		//~ChFunctionController();
		int Get_Type() { return 9527; }
		double Get_y(double curr_t);
		double Get_y_dx(double new_t) { return 0; }


		double angle_limit = 120*chrono::CH_C_PI/180;
		//statefeedback gains that give a second order closed loop response with natural undamped
		//frequency of 3 rads/sec and damping factor of 0.5 are: KCSS = [kcss1 kcss2]'
		double kcss1 = 6.75;  
		double kcss2 = 3.5;

		double friction = 1; //friction
		double u = 1;//motor input= (kcss1*voltage-KCSS'[angPos angVel]     [x1 x2]
		double r = 1;//motor  reference signal r(t)?
	protected:
		// The low level PID controller in motor.
		double ComputeOutput(double t);
		double cum_error_ = 0;
		double prevError = 0;
		Controller *controller_;
		size_t index_;
	};
} // namespace chrono

#endif // INCLUDE_CHFUNCTION_CONTROLLER_H_
