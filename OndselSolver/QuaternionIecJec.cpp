#include <memory>

#include "QuaternionIecJec.h"
#include "EndFrameqc.h"

using namespace MbD;

QuaternionIecJec::QuaternionIecJec() = default;

QuaternionIecJec::QuaternionIecJec(EndFrmsptr frmi, EndFrmsptr frmj, size_t axs) :
	KinematicIeJe(frmi, frmj), axis(axs)
{
}

void QuaternionIecJec::calcPostDynCorrectorIteration()
{
	// Compute imaginary component [axis] of conj(qI) * qJ
	// This uses the Hamilton product formula.
	// qE = [e0, e1, e2, e3] where e0,e1,e2 = imaginary (x,y,z), e3 = scalar (w)
	// conj(q) = [-q.x, -q.y, -q.z, q.w] = [-q[0], -q[1], -q[2], q[3]]
	//
	// IMPORTANT: Use qEO() to get the END FRAME's world orientation quaternion,
	// which accounts for the marker frame's relative orientation (aApm).
	// qEO = qEpart * qEmarker

	auto efrmI = std::static_pointer_cast<EndFrameqc>(frmI);
	auto efrmJ = std::static_pointer_cast<EndFrameqc>(frmJ);
	auto qEI = efrmI->qEO();  // Use world quaternion, not part quaternion
	auto qEJ = efrmJ->qEO();  // Use world quaternion, not part quaternion

	double qI0 = qEI->at(0);
	double qI1 = qEI->at(1);
	double qI2 = qEI->at(2);
	double qI3 = qEI->at(3);

	double qJ0 = qEJ->at(0);
	double qJ1 = qEJ->at(1);
	double qJ2 = qEJ->at(2);
	double qJ3 = qEJ->at(3);

	// Hamilton product p = conj(qI) * qJ
	// where conj(qI) = (-qI0, -qI1, -qI2, qI3)
	// p[0] = qI3*qJ0 - qI0*qJ3 - qI1*qJ2 + qI2*qJ1  (x imaginary)
	// p[1] = qI3*qJ1 + qI0*qJ2 - qI1*qJ3 - qI2*qJ0  (y imaginary)
	// p[2] = qI3*qJ2 - qI0*qJ1 + qI1*qJ0 - qI2*qJ3  (z imaginary)
	// p[3] = qI3*qJ3 + qI0*qJ0 + qI1*qJ1 + qI2*qJ2  (scalar)

	if (axis == 0) {
		aQijIeJe = qI3*qJ0 - qI0*qJ3 - qI1*qJ2 + qI2*qJ1;
	}
	else if (axis == 1) {
		aQijIeJe = qI3*qJ1 + qI0*qJ2 - qI1*qJ3 - qI2*qJ0;
	}
	else { // axis == 2
		aQijIeJe = qI3*qJ2 - qI0*qJ1 + qI1*qJ0 - qI2*qJ3;
	}
}

double MbD::QuaternionIecJec::value()
{
	return aQijIeJe;
}
