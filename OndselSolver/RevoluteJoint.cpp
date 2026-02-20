/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "RevoluteJoint.h"
#include "System.h"
#include "AtPointConstraintIJ.h"
#include "QuaternionConstraintIqcJqc.h"
#include "CREATE.h"

using namespace MbD;

RevoluteJoint::RevoluteJoint() 
{
}

RevoluteJoint::RevoluteJoint(const std::string& str) : AtPointJoint(str)
{
}

void RevoluteJoint::initializeGlobally()
{
	if (constraints->empty())
	{
		createAtPointConstraints();
		// Use quaternion-based orientation constraints instead of direction cosine constraints.
		// Quaternion constraints ensure Z axes are ALIGNED (same direction), not just parallel.
		// DirectionCosine constraints allowed anti-parallel (180° flipped) Z axes.
		// We constrain X and Y imaginary components to zero, allowing rotation around Z.
		addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 0));
		addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 1));
		this->root()->hasChanged = true;
	}
	else {
		Joint::initializeGlobally();
	}
}
