/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include "FixedJoint.h"
#include "QuaternionConstraintIqcJqc.h"
#include "EndFrameqc.h"
#include "MarkerFrame.h"
#include "PartFrame.h"
#include "System.h"
#include "CREATE.h"

using namespace MbD;

MbD::FixedJoint::FixedJoint()
{
}

MbD::FixedJoint::FixedJoint(const std::string& str) : AtPointJoint(str)
{
}

void MbD::FixedJoint::initializeGlobally()
{
	if (constraints->empty())
	{
		createAtPointConstraints();
		// Use quaternion-based orientation constraints instead of direction cosine constraints.
		// Quaternion constraints have only ONE valid solution (aligned orientation),
		// whereas direction cosine perpendicularity constraints have TWO solutions
		// (aligned and anti-parallel/180° rotated).
		addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 0));
		addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 1));
		addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 2));
		this->root()->hasChanged = true;
	}
	else {
		Joint::initializeGlobally();
	}
}
