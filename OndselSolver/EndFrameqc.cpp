/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include <memory>

#include "EndFrameqc.h"
#include "EndFrameqct.h"
#include "Variable.h"
#include "MarkerFrame.h"
#include "CREATE.h"
#include "EndFrameqct2.h"

using namespace MbD;

EndFrameqc::EndFrameqc() {
}

EndFrameqc::EndFrameqc(const std::string& str) : EndFramec(str) {
}

void EndFrameqc::initialize()
{
	prOeOpE = std::make_shared<FullMatrix<double>>(3, 4);
	pprOeOpEpE = std::make_shared<FullMatrix<FColDsptr>>(4, 4);
	pAOepE = std::make_shared<FullColumn<FMatDsptr>>(4);
	ppAOepEpE = std::make_shared<FullMatrix<FMatDsptr>>(4, 4);
}

void EndFrameqc::initializeGlobally()
{
	pprOeOpEpE = markerFrame->pprOmOpEpE;
	ppAOepEpE = markerFrame->ppAOmpEpE;
}

void EndFrameqc::initEndFrameqct()
{
	endFrameqct = CREATE<EndFrameqct>::With(this->name.data());
	endFrameqct->prOeOpE = prOeOpE;
	endFrameqct->pprOeOpEpE = pprOeOpEpE;
	endFrameqct->pAOepE = pAOepE;
	endFrameqct->ppAOepEpE = ppAOepEpE;
	endFrameqct->setMarkerFrame(markerFrame);
}

void MbD::EndFrameqc::initEndFrameqct2()
{
	endFrameqct = CREATE<EndFrameqct2>::With(this->name.data());
	endFrameqct->prOeOpE = prOeOpE;
	endFrameqct->pprOeOpEpE = pprOeOpEpE;
	endFrameqct->pAOepE = pAOepE;
	endFrameqct->ppAOepEpE = ppAOepEpE;
	endFrameqct->setMarkerFrame(markerFrame);
}

FMatFColDsptr EndFrameqc::ppAjOepEpE(size_t jj)
{
	auto answer = std::make_shared<FullMatrix<FColDsptr>>(4, 4);
	for (size_t i = 0; i < 4; i++) {
		auto& answeri = answer->at(i);
		auto& ppAOepEipE = ppAOepEpE->at(i);
		for (size_t j = i; j < 4; j++) {
			answeri->at(j) = ppAOepEipE->at(j)->column(jj);
		}
	}
	answer->symLowerWithUpper();
	return answer;
}

void EndFrameqc::calcPostDynCorrectorIteration()
{
	EndFramec::calcPostDynCorrectorIteration();
	prOeOpE = markerFrame->prOmOpE;
	pAOepE = markerFrame->pAOmpE;
}

FMatDsptr EndFrameqc::pAjOepET(size_t axis)
{
	auto answer = std::make_shared<FullMatrix<double>>(4, 3);
	for (size_t i = 0; i < 4; i++) {
		auto& answeri = answer->at(i);
		auto& pAOepEi = pAOepE->at(i);
		for (size_t j = 0; j < 3; j++) {
			auto& answerij = pAOepEi->at(j)->at(axis);
			answeri->at(j) = answerij;
		}
	}
	return answer;
}

FMatDsptr EndFrameqc::ppriOeOpEpE(size_t ii)
{
	auto answer = std::make_shared<FullMatrix<double>>(4, 4);
	for (size_t i = 0; i < 4; i++) {
		auto& answeri = answer->at(i);
		auto& pprOeOpEipE = pprOeOpEpE->at(i);
		for (size_t j = 0; j < 4; j++) {
			auto& answerij = pprOeOpEipE->at(j)->at(ii);
			answeri->at(j) = answerij;
		}
	}
	return answer;
}

size_t EndFrameqc::iqX()
{
	return markerFrame->iqX();
}

size_t EndFrameqc::iqE()
{
	return markerFrame->iqE();
}

FRowDsptr EndFrameqc::priOeOpE(size_t i)
{
	return prOeOpE->at(i);
}

FColDsptr EndFrameqc::qXdot()
{
	return markerFrame->qXdot();
}

std::shared_ptr<EulerParametersDot<double>> EndFrameqc::qEdot()
{
	return markerFrame->qEdot();
}

std::shared_ptr<EulerParameters<double>> EndFrameqc::qE()
{
	return markerFrame->qE();
}

std::shared_ptr<EulerParameters<double>> EndFrameqc::qEO()
{
	// Compute world orientation quaternion: qEworld = qEpart * qEmarker
	// This accounts for the marker frame's relative orientation to the part.
	auto qEpart = markerFrame->qE();
	auto qEmarker = markerFrame->qEpm;

	auto result = std::make_shared<EulerParameters<double>>(4);

	// Hamilton product: p = q1 * q2
	// OndselSolver convention: [e0, e1, e2, e3] = [x, y, z, w]
	// p[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1]
	// p[1] = q1[3]*q2[1] - q1[0]*q2[2] + q1[1]*q2[3] + q1[2]*q2[0]
	// p[2] = q1[3]*q2[2] + q1[0]*q2[1] - q1[1]*q2[0] + q1[2]*q2[3]
	// p[3] = q1[3]*q2[3] - q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2]

	double q1_0 = qEpart->at(0), q1_1 = qEpart->at(1), q1_2 = qEpart->at(2), q1_3 = qEpart->at(3);
	double q2_0 = qEmarker->at(0), q2_1 = qEmarker->at(1), q2_2 = qEmarker->at(2), q2_3 = qEmarker->at(3);

	result->at(0) = q1_3*q2_0 + q1_0*q2_3 + q1_1*q2_2 - q1_2*q2_1;
	result->at(1) = q1_3*q2_1 - q1_0*q2_2 + q1_1*q2_3 + q1_2*q2_0;
	result->at(2) = q1_3*q2_2 + q1_0*q2_1 - q1_1*q2_0 + q1_2*q2_3;
	result->at(3) = q1_3*q2_3 - q1_0*q2_0 - q1_1*q2_1 - q1_2*q2_2;

	return result;
}

FColDsptr EndFrameqc::qXddot()
{
	return markerFrame->qXddot();
}

FColDsptr EndFrameqc::qEddot()
{
	return markerFrame->qEddot();
}

FColDsptr EndFrameqc::rpep()
{
	return markerFrame->rpmp;
}

FColFMatDsptr EndFrameqc::pAOppE()
{
	return markerFrame->pAOppE();
}

FMatDsptr EndFrameqc::aBOp()
{
	return markerFrame->aBOp();
}

bool MbD::EndFrameqc::isEndFrameqc()
{
	return true;
}
