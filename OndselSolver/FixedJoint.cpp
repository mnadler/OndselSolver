/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <iostream>

#include "FixedJoint.h"
#include "QuaternionConstraintIqcJqc.h"
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
	std::cerr << "FixedJoint::initializeGlobally() START" << std::endl;
	if (constraints->empty())
	{
		std::cerr << "FixedJoint::initializeGlobally() creating constraints" << std::endl;
		createAtPointConstraints();
		std::cerr << "FixedJoint::initializeGlobally() AtPoint constraints created" << std::endl;
		// Use quaternion-based orientation constraints instead of direction cosine constraints.
		// Quaternion constraints have only ONE valid solution (aligned orientation),
		// whereas direction cosine perpendicularity constraints have TWO solutions
		// (aligned and anti-parallel/180° rotated).
		std::cerr << "FixedJoint::initializeGlobally() creating QuaternionConstraint 0" << std::endl;
		addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 0));
		std::cerr << "FixedJoint::initializeGlobally() creating QuaternionConstraint 1" << std::endl;
		addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 1));
		std::cerr << "FixedJoint::initializeGlobally() creating QuaternionConstraint 2" << std::endl;
		addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 2));
		std::cerr << "FixedJoint::initializeGlobally() all constraints created" << std::endl;
		this->root()->hasChanged = true;
	}
	else {
		std::cerr << "FixedJoint::initializeGlobally() using existing constraints" << std::endl;
		Joint::initializeGlobally();
	}
	std::cerr << "FixedJoint::initializeGlobally() DONE" << std::endl;
}
