#pragma once

#include <cstdint>

#include "QuaternionConstraintIJ.h"

namespace MbD {
	class QuaternionConstraintIqcJc : public QuaternionConstraintIJ
	{
		//pGpEI ppGpEIpEI iqEI
		//Frame I has Euler parameters (can vary), frame J is fixed
	public:
		QuaternionConstraintIqcJc(EndFrmsptr frmi, EndFrmsptr frmj, size_t axis);

		void addToJointTorqueI(FColDsptr col) override;
		void calcPostDynCorrectorIteration() override;
		void fillAccICIterError(FColDsptr col) override;
		void fillPosICError(FColDsptr col) override;
		void fillPosICJacob(SpMatDsptr mat) override;
		void fillPosKineJacob(SpMatDsptr mat) override;
		void fillVelICJacob(SpMatDsptr mat) override;
		void initaQijIeJe() override;
		void useEquationNumbers() override;

		FRowDsptr pGpEI;
		FMatDsptr ppGpEIpEI;
		size_t iqEI = SIZE_MAX;
	};
}
