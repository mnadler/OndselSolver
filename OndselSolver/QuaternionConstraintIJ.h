#pragma once

#include "ConstraintIJ.h"

namespace MbD {
	class QuaternionIecJec;

	class QuaternionConstraintIJ : public ConstraintIJ
	{
		//axis aQijIeJe
		//Constrains imaginary component [axis] of relative quaternion to zero
	public:
		QuaternionConstraintIJ(EndFrmsptr frmi, EndFrmsptr frmj, size_t axis);

		void calcPostDynCorrectorIteration() override;
		virtual void initaQijIeJe();
		void initialize() override;
		void initializeGlobally() override;
		void initializeLocally() override;
		void postInput() override;
		void postPosICIteration() override;
		void preAccIC() override;
		void prePosIC() override;
		void preVelIC() override;
		void simUpdateAll() override;
		ConstraintType type() override;

		size_t axis;
		std::shared_ptr<QuaternionIecJec> aQijIeJe;
	};
}
