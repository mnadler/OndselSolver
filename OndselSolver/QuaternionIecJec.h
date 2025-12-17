#pragma once

#include <memory>

#include "KinematicIeJe.h"

namespace MbD {

	class QuaternionIecJec : public KinematicIeJe
	{
		//aQijIeJe axis
		//Computes the imaginary component [axis] of conj(qE_I) * qE_J
		//This equals zero when frames I and J have identical orientation
	public:
		QuaternionIecJec();
		QuaternionIecJec(EndFrmsptr frmi, EndFrmsptr frmj, size_t axis);

		void calcPostDynCorrectorIteration() override;
		double value() override;

		size_t axis{};   //0, 1, 2 = x, y, z imaginary component
		double aQijIeJe{};  //Constraint value: Im[axis](conj(qI) * qJ)
	};
}
