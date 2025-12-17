#include "QuaternionIeqcJec.h"
#include "EndFrameqc.h"
#include "MarkerFrame.h"

using namespace MbD;

QuaternionIeqcJec::QuaternionIeqcJec() = default;

QuaternionIeqcJec::QuaternionIeqcJec(EndFrmsptr frmi, EndFrmsptr frmj, size_t axs) :
	QuaternionIecJec(frmi, frmj, axs)
{
}

void QuaternionIeqcJec::initialize()
{
	QuaternionIecJec::initialize();
	pQijIeJepEI = std::make_shared<FullRow<double>>(4);
	ppQijIeJepEIpEI = std::make_shared<FullMatrix<double>>(4, 4);
}

FRowDsptr QuaternionIeqcJec::pvaluepEI()
{
	return pQijIeJepEI;
}

FMatDsptr QuaternionIeqcJec::ppvaluepEIpEI()
{
	return ppQijIeJepEIpEI;
}

void QuaternionIeqcJec::calcPostDynCorrectorIteration()
{
	QuaternionIecJec::calcPostDynCorrectorIteration();

	// Compute ∂C/∂qI_part using chain rule:
	// C = f(qI_world, qJ_world) where qI_world = qI_part * qI_marker
	// ∂C/∂qI_part = (∂C/∂qI_world) * (∂qI_world/∂qI_part)

	auto efrmI = std::static_pointer_cast<EndFrameqc>(frmI);
	auto efrmJ = std::static_pointer_cast<EndFrameqc>(frmJ);

	// Get world quaternion of J (accounts for marker frame)
	auto qJ_world = efrmJ->qEO();
	double qJw0 = qJ_world->at(0);
	double qJw1 = qJ_world->at(1);
	double qJw2 = qJ_world->at(2);
	double qJw3 = qJ_world->at(3);

	// Get marker's relative quaternion for I
	auto qI_marker = efrmI->markerFrame->qEpm;
	double qIm0 = qI_marker->at(0);
	double qIm1 = qI_marker->at(1);
	double qIm2 = qI_marker->at(2);
	double qIm3 = qI_marker->at(3);

	// First compute ∂C/∂qI_world (1x4 row vector)
	// C = conj(qI_world) * qJ_world, taking imaginary component [axis]
	// p[0] = qIw3*qJw0 - qIw0*qJw3 - qIw1*qJw2 + qIw2*qJw1
	// p[1] = qIw3*qJw1 + qIw0*qJw2 - qIw1*qJw3 - qIw2*qJw0
	// p[2] = qIw3*qJw2 - qIw0*qJw1 + qIw1*qJw0 - qIw2*qJw3

	double pCpqIw[4];  // ∂C/∂qI_world
	if (axis == 0) {
		pCpqIw[0] = -qJw3; pCpqIw[1] = -qJw2; pCpqIw[2] = qJw1; pCpqIw[3] = qJw0;
	}
	else if (axis == 1) {
		pCpqIw[0] = qJw2; pCpqIw[1] = -qJw3; pCpqIw[2] = -qJw0; pCpqIw[3] = qJw1;
	}
	else { // axis == 2
		pCpqIw[0] = -qJw1; pCpqIw[1] = qJw0; pCpqIw[2] = -qJw3; pCpqIw[3] = qJw2;
	}

	// Now compute ∂qI_world/∂qI_part (4x4 matrix)
	// qI_world = qI_part * qI_marker (Hamilton product)
	// For p = q1 * q2 with q2 constant:
	// ∂p/∂q1 where p[i] depends on q1[j]:
	// p[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1]
	// p[1] = q1[3]*q2[1] - q1[0]*q2[2] + q1[1]*q2[3] + q1[2]*q2[0]
	// p[2] = q1[3]*q2[2] + q1[0]*q2[1] - q1[1]*q2[0] + q1[2]*q2[3]
	// p[3] = q1[3]*q2[3] - q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2]
	//
	// M[i][j] = ∂p[i]/∂q1[j]
	// M[0] = [q2[3], q2[2], -q2[1], q2[0]]
	// M[1] = [-q2[2], q2[3], q2[0], q2[1]]
	// M[2] = [q2[1], -q2[0], q2[3], q2[2]]
	// M[3] = [-q2[0], -q2[1], -q2[2], q2[3]]

	double M[4][4];
	M[0][0] = qIm3;  M[0][1] = qIm2;  M[0][2] = -qIm1; M[0][3] = qIm0;
	M[1][0] = -qIm2; M[1][1] = qIm3;  M[1][2] = qIm0;  M[1][3] = qIm1;
	M[2][0] = qIm1;  M[2][1] = -qIm0; M[2][2] = qIm3;  M[2][3] = qIm2;
	M[3][0] = -qIm0; M[3][1] = -qIm1; M[3][2] = -qIm2; M[3][3] = qIm3;

	// Apply chain rule: ∂C/∂qI_part = (∂C/∂qI_world) * (∂qI_world/∂qI_part)
	// Result is 1x4: (1x4) * (4x4)
	// pCpqIp[j] = sum_i pCpqIw[i] * M[i][j]
	for (size_t j = 0; j < 4; j++) {
		double sum = 0.0;
		for (size_t i = 0; i < 4; i++) {
			sum += pCpqIw[i] * M[i][j];
		}
		pQijIeJepEI->at(j) = sum;
	}

	// Second derivatives: Since C is bilinear in qI_world and qJ_world,
	// and qI_world is linear in qI_part, ∂²C/∂qI_part² = 0
	// ppQijIeJepEIpEI stays zero (initialized)
}
