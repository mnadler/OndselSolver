#pragma once

#include "QuaternionIeqcJec.h"

namespace MbD {

	class QuaternionIeqcJeqc : public QuaternionIeqcJec
	{
		//pQijIeJepEJ ppQijIeJepEIpEJ ppQijIeJepEJpEJ
		//Adds frame J Euler parameter derivatives
	public:
		QuaternionIeqcJeqc();
		QuaternionIeqcJeqc(EndFrmsptr frmi, EndFrmsptr frmj, size_t axis);

		void initialize() override;
		void initializeGlobally() override;
		void calcPostDynCorrectorIteration() override;
		FRowDsptr pvaluepEJ() override;
		FMatDsptr ppvaluepEIpEJ() override;
		FMatDsptr ppvaluepEJpEJ() override;

		FRowDsptr pQijIeJepEJ;      // 1x4: ∂C/∂qE_J
		FMatDsptr ppQijIeJepEIpEJ;  // 4x4: ∂²C/∂qE_I∂qE_J (constant)
		FMatDsptr ppQijIeJepEJpEJ;  // 4x4: ∂²C/∂qE_J² (zero for bilinear)
	};
}
