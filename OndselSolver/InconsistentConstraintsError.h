#pragma once

#include <stdexcept>
#include <memory>
#include <vector>
#include <string>
#include <array>

#include "FullColumn.h"

namespace MbD {

	struct LCSDiagnostic {
		std::string name;
		std::array<double, 3> positionOnPart;
		std::array<double, 3> worldPosition;
		std::array<double, 4> worldQuaternion = {0, 0, 0, 1};  // [x, y, z, w]
	};

	struct ConstraintDiagnostic {
		size_t equationNumber;
		std::string type;
		double violation;
	};

	struct JointDiagnostic {
		std::string name;
		std::string type;
		std::string partIName;
		std::string partJName;
		LCSDiagnostic lcsI;
		LCSDiagnostic lcsJ;
		double relativeAngleDegrees = 0.0;  // Angle between frames
		std::array<double, 4> relativeQuaternion = {0, 0, 0, 1};  // conj(qI) * qJ: [x, y, z, w]
		std::vector<ConstraintDiagnostic> inconsistentConstraints;
	};

	struct PartConstraintDiagnostic {
		std::string partName;
		std::string constraintType;
		size_t equationNumber;
		double violation;
	};

	struct InconsistencyDiagnostic {
		std::string affectedPartName;
		double totalViolation = 0.0;
		std::vector<JointDiagnostic> joints;
		std::vector<PartConstraintDiagnostic> partConstraints;
	};

	class InconsistentConstraintsError : virtual public std::runtime_error
	{
	protected:
		std::shared_ptr<std::vector<size_t>> inconsistentEqnNos;
		std::shared_ptr<std::vector<double>> rhsValues;
		std::shared_ptr<InconsistencyDiagnostic> diagnostic;

	public:
		explicit
			InconsistentConstraintsError(const std::string& msg, std::shared_ptr<FullColumn<size_t>> eqnNos) :
			std::runtime_error(msg), inconsistentEqnNos(eqnNos)
		{
		}

		explicit
			InconsistentConstraintsError(const std::string& msg, std::shared_ptr<FullColumn<size_t>> eqnNos,
				std::shared_ptr<std::vector<double>> rhs) :
			std::runtime_error(msg), inconsistentEqnNos(eqnNos), rhsValues(rhs)
		{
		}

		explicit InconsistentConstraintsError(const std::string& msg) : std::runtime_error(msg)
		{
		}

		virtual ~InconsistentConstraintsError() noexcept {}

		virtual std::shared_ptr<std::vector<size_t>> getInconsistentEqnNos() const noexcept {
			return inconsistentEqnNos;
		}

		virtual std::shared_ptr<std::vector<double>> getRhsValues() const noexcept {
			return rhsValues;
		}

		void setDiagnostic(std::shared_ptr<InconsistencyDiagnostic> diag) {
			diagnostic = diag;
		}

		std::shared_ptr<InconsistencyDiagnostic> getDiagnostic() const noexcept {
			return diagnostic;
		}
	};
}
