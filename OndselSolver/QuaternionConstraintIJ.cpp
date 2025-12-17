#include <iostream>

#include "QuaternionConstraintIJ.h"
#include "QuaternionIecJec.h"
#include "EndFramec.h"
#include "CREATE.h"

using namespace MbD;

QuaternionConstraintIJ::QuaternionConstraintIJ(EndFrmsptr frmi, EndFrmsptr frmj, size_t axs) :
	ConstraintIJ(frmi, frmj), axis(axs)
{
	std::cerr << "QuaternionConstraintIJ constructor DONE axis=" << axs << std::endl;
}

void QuaternionConstraintIJ::initialize()
{
	std::cerr << "QuaternionConstraintIJ::initialize() axis=" << axis << std::endl;
	ConstraintIJ::initialize();
	initaQijIeJe();
	std::cerr << "QuaternionConstraintIJ::initialize() DONE" << std::endl;
}

void QuaternionConstraintIJ::initializeLocally()
{
	std::cerr << "QuaternionConstraintIJ::initializeLocally() axis=" << axis << std::endl;
	aQijIeJe->initializeLocally();
	std::cerr << "QuaternionConstraintIJ::initializeLocally() DONE" << std::endl;
}

void QuaternionConstraintIJ::initializeGlobally()
{
	std::cerr << "QuaternionConstraintIJ::initializeGlobally() axis=" << axis << std::endl;
	aQijIeJe->initializeGlobally();
	std::cerr << "QuaternionConstraintIJ::initializeGlobally() DONE" << std::endl;
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
	std::cerr << "QuaternionConstraintIJ::calcPostDynCorrectorIteration() axis=" << axis << std::endl;
	aG = aQijIeJe->aQijIeJe - aConstant;
	std::cerr << "QuaternionConstraintIJ::calcPostDynCorrectorIteration() DONE aG=" << aG << std::endl;
}

void QuaternionConstraintIJ::prePosIC()
{
	aQijIeJe->prePosIC();
	ConstraintIJ::prePosIC();
}

void QuaternionConstraintIJ::postPosICIteration()
{
	std::cerr << "QuaternionConstraintIJ::postPosICIteration() axis=" << axis << std::endl;
	aQijIeJe->postPosICIteration();
	std::cerr << "QuaternionConstraintIJ::postPosICIteration() after aQijIeJe" << std::endl;
	ConstraintIJ::postPosICIteration();
	std::cerr << "QuaternionConstraintIJ::postPosICIteration() DONE" << std::endl;
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
