#include "QuaternionConstraintIqcJqc.h"
#include "QuaternionIeqcJeqc.h"
#include "EndFrameqc.h"
#include "CREATE.h"

using namespace MbD;

QuaternionConstraintIqcJqc::QuaternionConstraintIqcJqc(EndFrmsptr frmi, EndFrmsptr frmj, size_t axs) :
	QuaternionConstraintIqcJc(frmi, frmj, axs)
{
}

void QuaternionConstraintIqcJqc::initaQijIeJe()
{
	aQijIeJe = CREATE<QuaternionIeqcJeqc>::With(frmI, frmJ, axis);
}

void QuaternionConstraintIqcJqc::calcPostDynCorrectorIteration()
{
	QuaternionConstraintIqcJc::calcPostDynCorrectorIteration();
	auto aQijIeqJqe = std::static_pointer_cast<QuaternionIeqcJeqc>(aQijIeJe);
	pGpEJ = aQijIeqJqe->pQijIeJepEJ;
	ppGpEIpEJ = aQijIeqJqe->ppQijIeJepEIpEJ;
	ppGpEJpEJ = aQijIeqJqe->ppQijIeJepEJpEJ;
}

void QuaternionConstraintIqcJqc::useEquationNumbers()
{
	QuaternionConstraintIqcJc::useEquationNumbers();
	iqEJ = std::static_pointer_cast<EndFrameqc>(frmJ)->iqE();
}

std::string QuaternionConstraintIqcJqc::constraintSpec()
{
	return "QuaternionConstraint" + MbDMath::xyzFromInt(axis);
}

void QuaternionConstraintIqcJqc::fillPosICError(FColDsptr col)
{
	QuaternionConstraintIqcJc::fillPosICError(col);
	col->atiplusFullVectortimes(iqEJ, pGpEJ, lam);
}

void QuaternionConstraintIqcJqc::fillPosICJacob(SpMatDsptr mat)
{
	QuaternionConstraintIqcJc::fillPosICJacob(mat);
	mat->atijplusFullRow(iG, iqEJ, pGpEJ);
	mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
	auto ppGpEIpEJlam = ppGpEIpEJ->times(lam);
	mat->atijplusFullMatrix(iqEI, iqEJ, ppGpEIpEJlam);
	mat->atijplusTransposeFullMatrix(iqEJ, iqEI, ppGpEIpEJlam);
	mat->atijplusFullMatrixtimes(iqEJ, iqEJ, ppGpEJpEJ, lam);
}

void QuaternionConstraintIqcJqc::fillPosKineJacob(SpMatDsptr mat)
{
	QuaternionConstraintIqcJc::fillPosKineJacob(mat);
	mat->atijplusFullRow(iG, iqEJ, pGpEJ);
}

void QuaternionConstraintIqcJqc::fillVelICJacob(SpMatDsptr mat)
{
	QuaternionConstraintIqcJc::fillVelICJacob(mat);
	mat->atijplusFullRow(iG, iqEJ, pGpEJ);
	mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
}

void QuaternionConstraintIqcJqc::fillAccICIterError(FColDsptr col)
{
	QuaternionConstraintIqcJc::fillAccICIterError(col);
	col->atiplusFullVectortimes(iqEJ, pGpEJ, lam);
	auto efrmIqc = std::static_pointer_cast<EndFrameqc>(frmI);
	auto efrmJqc = std::static_pointer_cast<EndFrameqc>(frmJ);
	auto qEdotI = efrmIqc->qEdot();
	auto qEdotJ = efrmJqc->qEdot();
	double sum = pGpEJ->timesFullColumn(efrmJqc->qEddot());
	sum += (qEdotI->transposeTimesFullColumn(ppGpEIpEJ->timesFullColumn(qEdotJ))) * 2.0;
	sum += qEdotJ->transposeTimesFullColumn(ppGpEJpEJ->timesFullColumn(qEdotJ));
	col->atiplusNumber(iG, sum);
}
