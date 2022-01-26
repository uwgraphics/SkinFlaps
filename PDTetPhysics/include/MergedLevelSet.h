#pragma once
#include <vector>
#include <string>

#include "PhysBAM_Tools/Vectors/VECTOR.h"

namespace PhysBAM {
	template<class VectorType>
		class LEVELSET_IMPLICIT_OBJECT;

	template <class VectorType>
	class MergedLevelSet
	{
		//using VectorType = VECTOR<T, d>;
	private:
		using T = typename VectorType::SCALAR;
		static constexpr int d = VectorType::m;
		std::vector<PhysBAM::LEVELSET_IMPLICIT_OBJECT<VectorType>*>m_levelSet;
		T gridDX;
	public:
		T Extended_Phi(const VectorType& pos) const;
		VectorType Extended_Normal(const VectorType& pos) const;

		void initializeLevelSet(std::vector<std::string>& paths, const T dx = .025) {
			gridDX = dx;
			for (const auto& s : paths) addLevelSet(s);
		}

		~MergedLevelSet();

	private:
		void addLevelSet(const std::string& collisionObjPath);
	};
}

