#include "QuaternionConstraintIJ.h"
#include "QuaternionIecJec.h"
#include "EndFramec.h"
#include "CREATE.h"

using namespace MbD;

QuaternionConstraintIJ::QuaternionConstraintIJ(EndFrmsptr frmi, EndFrmsptr frmj, size_t axs) :
	ConstraintIJ(frmi, frmj), axis(axs)
{
}

void QuaternionConstraintIJ::initialize()
{
	ConstraintIJ::initialize();
	initaQijIeJe();
}

void QuaternionConstraintIJ::initializeLocally()
{
	aQijIeJe->initializeLocally();
}

void QuaternionConstraintIJ::initializeGlobally()
{
	aQijIeJe->initializeGlobally();
}

void QuaternionConstraintIJ::initaQijIeJe()
{
	aQijIeJe = CREATE<QuaternionIecJec>::With(frmI, frmJ, axis);
}

void QuaternionConstraintIJ::postInput()
{
	aQijIeJe->postInput();
	ConstraintIJ::postInput();
}

void QuaternionConstraintIJ::calcPostDynCorrectorIteration()
{
	aG = aQijIeJe->aQijIeJe - aConstant;
}

void QuaternionConstraintIJ::prePosIC()
{
	aQijIeJe->prePosIC();
	ConstraintIJ::prePosIC();
}

void QuaternionConstraintIJ::postPosICIteration()
{
	aQijIeJe->postPosICIteration();
	ConstraintIJ::postPosICIteration();
}

ConstraintType QuaternionConstraintIJ::type()
{
	return perpendicular;
}

void QuaternionConstraintIJ::preVelIC()
{
	aQijIeJe->preVelIC();
	ConstraintIJ::preVelIC();
}

void QuaternionConstraintIJ::simUpdateAll()
{
	aQijIeJe->simUpdateAll();
	ConstraintIJ::simUpdateAll();
}

void QuaternionConstraintIJ::preAccIC()
{
	aQijIeJe->preAccIC();
	ConstraintIJ::preAccIC();
}
