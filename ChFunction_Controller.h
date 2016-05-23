#ifndef INCLUDE_CHFUNCTION_CONTROLLER_H_
#define INCLUDE_CHFUNCTION_CONTROLLER_H_

#include <cmath>
#include <algorithm>

#include <motion_functions/ChFunction_Base.h>

namespace chrono {
	class Controller;

	class ChFunctionController : public ChFunction {
	public:
		ChFunctionController(size_t index, Controller* controller)
			: index_(index), controller_(controller) {}
		virtual ~ChFunctionController(){};
		virtual ChFunction *new_Duplicate() {
			return 0;
			//new ChFunctionController(index_, controller_);

		}
		virtual ChFunctionController* Clone() const override { return new ChFunctionController(*this); }
		virtual double Get_y(double curr_t) const override;

		void ResetCumulative(double t);
		int Get_Type() { return 9527; }
		double Get_y2(double curr_t);
		double Get_y(double curr_t);
		double Get_y_dx(double new_t) { return 0;};

	protected:
		double ComputeOutput(double t);
		double OutputToOmega(double t, double out);
		double OmegaToTorque(double t, double out);
		//double cum_error_ = 0;
		//double prevError = 0;
		Controller *controller_;
		size_t index_;
	};
} // namespace chrono

#endif // INCLUDE_CHFUNCTION_CONTROLLER_H_
