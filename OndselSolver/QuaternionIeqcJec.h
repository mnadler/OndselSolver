#pragma once

#include "QuaternionIecJec.h"

namespace MbD {

	class QuaternionIeqcJec : public QuaternionIecJec
	{
		//pQijIeJepEI ppQijIeJepEIpEI
		//Adds frame I Euler parameter derivatives
	public:
		QuaternionIeqcJec();
		QuaternionIeqcJec(EndFrmsptr frmi, EndFrmsptr frmj, size_t axis);

		void initialize() override;
		void calcPostDynCorrectorIteration() override;
		FRowDsptr pvaluepEI() override;
		FMatDsptr ppvaluepEIpEI() override;

		FRowDsptr pQijIeJepEI;      // 1x4: ∂C/∂qE_I
		FMatDsptr ppQijIeJepEIpEI;  // 4x4: ∂²C/∂qE_I² (zero for bilinear)
	};
}
