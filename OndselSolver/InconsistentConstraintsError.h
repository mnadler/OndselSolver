#pragma once

#include <stdexcept>
#include <memory>
#include <vector>

#include "FullColumn.h"

namespace MbD {
	class InconsistentConstraintsError : virtual public std::runtime_error
	{
	protected:
		std::shared_ptr<std::vector<size_t>> inconsistentEqnNos;

	public:
		explicit
			InconsistentConstraintsError(const std::string& msg, std::shared_ptr<FullColumn<size_t>> eqnNos) :
			std::runtime_error(msg), inconsistentEqnNos(eqnNos)
		{
		}
		explicit InconsistentConstraintsError(const std::string& msg) : std::runtime_error(msg)
		{
		}

		virtual ~InconsistentConstraintsError() noexcept {}

		virtual std::shared_ptr<std::vector<size_t>> getInconsistentEqnNos() const noexcept {
			return inconsistentEqnNos;
		}
	};
}
