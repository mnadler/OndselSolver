#include "QuaternionIeqcJeqc.h"
#include "EndFrameqc.h"
#include "MarkerFrame.h"

using namespace MbD;

QuaternionIeqcJeqc::QuaternionIeqcJeqc() = default;

QuaternionIeqcJeqc::QuaternionIeqcJeqc(EndFrmsptr frmi, EndFrmsptr frmj, size_t axs) :
	QuaternionIeqcJec(frmi, frmj, axs)
{
}

void QuaternionIeqcJeqc::initialize()
{
	QuaternionIeqcJec::initialize();
	pQijIeJepEJ = std::make_shared<FullRow<double>>(4);
	ppQijIeJepEIpEJ = std::make_shared<FullMatrix<double>>(4, 4);
	ppQijIeJepEJpEJ = std::make_shared<FullMatrix<double>>(4, 4);
}

void QuaternionIeqcJeqc::initializeGlobally()
{
	// Compute the constant cross-derivative matrix ∂²C/∂qI_part∂qJ_part
	// With chain rule: ∂²C/∂qI_part∂qJ_part = Mi^T * H_world * Mj
	// Where:
	// - Mi = ∂qI_world/∂qI_part (4x4)
	// - Mj = ∂qJ_world/∂qJ_part (4x4)
	// - H_world = ∂²C/∂qI_world∂qJ_world (4x4, constant)

	auto efrmI = std::static_pointer_cast<EndFrameqc>(frmI);
	auto efrmJ = std::static_pointer_cast<EndFrameqc>(frmJ);

	// Get marker quaternions
	auto qI_marker = efrmI->markerFrame->qEpm;
	auto qJ_marker = efrmJ->markerFrame->qEpm;

	double qIm0 = qI_marker->at(0), qIm1 = qI_marker->at(1), qIm2 = qI_marker->at(2), qIm3 = qI_marker->at(3);
	double qJm0 = qJ_marker->at(0), qJm1 = qJ_marker->at(1), qJm2 = qJ_marker->at(2), qJm3 = qJ_marker->at(3);

	// Build Mi = ∂qI_world/∂qI_part (4x4)
	double Mi[4][4];
	Mi[0][0] = qIm3;  Mi[0][1] = qIm2;  Mi[0][2] = -qIm1; Mi[0][3] = qIm0;
	Mi[1][0] = -qIm2; Mi[1][1] = qIm3;  Mi[1][2] = qIm0;  Mi[1][3] = qIm1;
	Mi[2][0] = qIm1;  Mi[2][1] = -qIm0; Mi[2][2] = qIm3;  Mi[2][3] = qIm2;
	Mi[3][0] = -qIm0; Mi[3][1] = -qIm1; Mi[3][2] = -qIm2; Mi[3][3] = qIm3;

	// Build Mj = ∂qJ_world/∂qJ_part (4x4)
	double Mj[4][4];
	Mj[0][0] = qJm3;  Mj[0][1] = qJm2;  Mj[0][2] = -qJm1; Mj[0][3] = qJm0;
	Mj[1][0] = -qJm2; Mj[1][1] = qJm3;  Mj[1][2] = qJm0;  Mj[1][3] = qJm1;
	Mj[2][0] = qJm1;  Mj[2][1] = -qJm0; Mj[2][2] = qJm3;  Mj[2][3] = qJm2;
	Mj[3][0] = -qJm0; Mj[3][1] = -qJm1; Mj[3][2] = -qJm2; Mj[3][3] = qJm3;

	// Build H_world = ∂²C/∂qI_world∂qJ_world (constant 4x4)
	double H[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

	if (axis == 0) {
		H[3][0] = 1.0; H[0][3] = -1.0; H[1][2] = -1.0; H[2][1] = 1.0;
	}
	else if (axis == 1) {
		H[3][1] = 1.0; H[0][2] = 1.0; H[1][3] = -1.0; H[2][0] = -1.0;
	}
	else { // axis == 2
		H[3][2] = 1.0; H[0][1] = -1.0; H[1][0] = 1.0; H[2][3] = -1.0;
	}

	// Compute result = Mi^T * H * Mj
	// First: temp = H * Mj (4x4)
	double temp[4][4];
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			temp[i][j] = 0.0;
			for (size_t k = 0; k < 4; k++) {
				temp[i][j] += H[i][k] * Mj[k][j];
			}
		}
	}

	// Then: result = Mi^T * temp (4x4)
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			double val = 0.0;
			for (size_t k = 0; k < 4; k++) {
				val += Mi[k][i] * temp[k][j];  // Mi^T: swap indices
			}
			ppQijIeJepEIpEJ->at(i)->at(j) = val;
		}
	}

	// ppQijIeJepEJpEJ stays zero (linear in qJ_world, and qJ_world linear in qJ_part)
}

FRowDsptr QuaternionIeqcJeqc::pvaluepEJ()
{
	return pQijIeJepEJ;
}

FMatDsptr QuaternionIeqcJeqc::ppvaluepEIpEJ()
{
	return ppQijIeJepEIpEJ;
}

FMatDsptr QuaternionIeqcJeqc::ppvaluepEJpEJ()
{
	return ppQijIeJepEJpEJ;
}

void QuaternionIeqcJeqc::calcPostDynCorrectorIteration()
{
	QuaternionIeqcJec::calcPostDynCorrectorIteration();

	// Compute ∂C/∂qJ_part using chain rule:
	// C = f(qI_world, qJ_world) where qJ_world = qJ_part * qJ_marker
	// ∂C/∂qJ_part = (∂C/∂qJ_world) * (∂qJ_world/∂qJ_part)

	auto efrmI = std::static_pointer_cast<EndFrameqc>(frmI);
	auto efrmJ = std::static_pointer_cast<EndFrameqc>(frmJ);

	// Get world quaternion of I (accounts for marker frame)
	auto qI_world = efrmI->qEO();
	double qIw0 = qI_world->at(0);
	double qIw1 = qI_world->at(1);
	double qIw2 = qI_world->at(2);
	double qIw3 = qI_world->at(3);

	// Get marker's relative quaternion for J
	auto qJ_marker = efrmJ->markerFrame->qEpm;
	double qJm0 = qJ_marker->at(0);
	double qJm1 = qJ_marker->at(1);
	double qJm2 = qJ_marker->at(2);
	double qJm3 = qJ_marker->at(3);

	// First compute ∂C/∂qJ_world (1x4 row vector)
	// C = conj(qI_world) * qJ_world, taking imaginary component [axis]
	// p[0] = qIw3*qJw0 - qIw0*qJw3 - qIw1*qJw2 + qIw2*qJw1
	// p[1] = qIw3*qJw1 + qIw0*qJw2 - qIw1*qJw3 - qIw2*qJw0
	// p[2] = qIw3*qJw2 - qIw0*qJw1 + qIw1*qJw0 - qIw2*qJw3

	double pCpqJw[4];  // ∂C/∂qJ_world
	if (axis == 0) {
		pCpqJw[0] = qIw3; pCpqJw[1] = qIw2; pCpqJw[2] = -qIw1; pCpqJw[3] = -qIw0;
	}
	else if (axis == 1) {
		pCpqJw[0] = -qIw2; pCpqJw[1] = qIw3; pCpqJw[2] = qIw0; pCpqJw[3] = -qIw1;
	}
	else { // axis == 2
		pCpqJw[0] = qIw1; pCpqJw[1] = -qIw0; pCpqJw[2] = qIw3; pCpqJw[3] = -qIw2;
	}

	// Now compute ∂qJ_world/∂qJ_part (4x4 matrix)
	// qJ_world = qJ_part * qJ_marker (Hamilton product)
	// Same structure as for I
	double Mj[4][4];
	Mj[0][0] = qJm3;  Mj[0][1] = qJm2;  Mj[0][2] = -qJm1; Mj[0][3] = qJm0;
	Mj[1][0] = -qJm2; Mj[1][1] = qJm3;  Mj[1][2] = qJm0;  Mj[1][3] = qJm1;
	Mj[2][0] = qJm1;  Mj[2][1] = -qJm0; Mj[2][2] = qJm3;  Mj[2][3] = qJm2;
	Mj[3][0] = -qJm0; Mj[3][1] = -qJm1; Mj[3][2] = -qJm2; Mj[3][3] = qJm3;

	// Apply chain rule: ∂C/∂qJ_part = (∂C/∂qJ_world) * (∂qJ_world/∂qJ_part)
	// Result is 1x4: (1x4) * (4x4)
	// pCpqJp[j] = sum_i pCpqJw[i] * Mj[i][j]
	for (size_t j = 0; j < 4; j++) {
		double sum = 0.0;
		for (size_t i = 0; i < 4; i++) {
			sum += pCpqJw[i] * Mj[i][j];
		}
		pQijIeJepEJ->at(j) = sum;
	}

	// ppQijIeJepEJpEJ stays zero (linear in qJ_world, and qJ_world linear in qJ_part)
}
