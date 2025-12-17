#include <iostream>

#include "QuaternionConstraintIqcJc.h"
#include "QuaternionIeqcJec.h"
#include "EndFrameqc.h"
#include "CREATE.h"

using namespace MbD;

QuaternionConstraintIqcJc::QuaternionConstraintIqcJc(EndFrmsptr frmi, EndFrmsptr frmj, size_t axs) :
	QuaternionConstraintIJ(frmi, frmj, axs)
{
	std::cerr << "QuaternionConstraintIqcJc constructor DONE axis=" << axs << std::endl;
}

void QuaternionConstraintIqcJc::initaQijIeJe()
{
	aQijIeJe = CREATE<QuaternionIeqcJec>::With(frmI, frmJ, axis);
}

void QuaternionConstraintIqcJc::calcPostDynCorrectorIteration()
{
	QuaternionConstraintIJ::calcPostDynCorrectorIteration();
	auto aQijIeqJe = std::static_pointer_cast<QuaternionIeqcJec>(aQijIeJe);
	pGpEI = aQijIeqJe->pQijIeJepEI;
	ppGpEIpEI = aQijIeqJe->ppQijIeJepEIpEI;
}

void QuaternionConstraintIqcJc::useEquationNumbers()
{
	iqEI = std::static_pointer_cast<EndFrameqc>(frmI)->iqE();
}

void QuaternionConstraintIqcJc::fillPosICError(FColDsptr col)
{
	Constraint::fillPosICError(col);
	col->atiplusFullVectortimes(iqEI, pGpEI, lam);
}

void QuaternionConstraintIqcJc::fillPosICJacob(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqEI, pGpEI);
	mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
	mat->atijplusFullMatrixtimes(iqEI, iqEI, ppGpEIpEI, lam);
}

void QuaternionConstraintIqcJc::fillPosKineJacob(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqEI, pGpEI);
}

void QuaternionConstraintIqcJc::fillVelICJacob(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqEI, pGpEI);
	mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
}

void QuaternionConstraintIqcJc::fillAccICIterError(FColDsptr col)
{
	col->atiplusFullVector(iqEI, pGpEI->times(lam));
	auto efrmIqc = std::static_pointer_cast<EndFrameqc>(frmI);
	auto qEdotI = efrmIqc->qEdot();
	auto sum = pGpEI->timesFullColumn(efrmIqc->qEddot());
	sum += qEdotI->transposeTimesFullColumn(ppGpEIpEI->timesFullColumn(qEdotI));
	col->atiplusNumber(iG, sum);
}

void QuaternionConstraintIqcJc::addToJointTorqueI(FColDsptr jointTorque)
{
	auto aBOIp = frmI->aBOp();
	auto lampGpE = pGpEI->transpose()->times(lam);
	auto c2Torque = aBOIp->timesFullColumn(lampGpE);
	jointTorque->equalSelfPlusFullColumntimes(c2Torque, 0.5);
}
