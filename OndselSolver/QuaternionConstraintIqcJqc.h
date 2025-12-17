#pragma once

#include <cstdint>

#include "QuaternionConstraintIqcJc.h"

namespace MbD {
	class QuaternionConstraintIqcJqc : public QuaternionConstraintIqcJc
	{
		//pGpEJ ppGpEIpEJ ppGpEJpEJ iqEJ
		//Both frames have Euler parameters (can vary)
	public:
		QuaternionConstraintIqcJqc(EndFrmsptr frmi, EndFrmsptr frmj, size_t axis);

		void calcPostDynCorrectorIteration() override;
		void fillAccICIterError(FColDsptr col) override;
		void fillPosICError(FColDsptr col) override;
		void fillPosICJacob(SpMatDsptr mat) override;
		void fillPosKineJacob(SpMatDsptr mat) override;
		void fillVelICJacob(SpMatDsptr mat) override;
		void initaQijIeJe() override;
		void useEquationNumbers() override;
		std::string constraintSpec() override;

		FRowDsptr pGpEJ;
		FMatDsptr ppGpEIpEJ;
		FMatDsptr ppGpEJpEJ;
		size_t iqEJ = SIZE_MAX;
	};
}
