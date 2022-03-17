#include "GridDeformerTet.h"
#include "Add_Force.h"


#include <omp.h>

#include "dumper.h"


namespace {

}

namespace PhysBAM {

    template <class dataType, int dim>
    void GridDeformerTet<std::vector<VECTOR<dataType,dim>>>::initializeDeformer() {
        //LOG::SCOPE scope("GridDeformerTet::initializeDeformer()");
        m_elementFlags.resize(m_elements.size());
        for (auto & f:m_elementFlags)
            f=ElementFlag::unCollisionEl;
    }

    template <class dataType, int dim>
    void GridDeformerTet<std::vector<VECTOR<dataType,dim>>>::initializeUndeformedState() {
        // LOG::SCOPE scope("GridDeformerTet::initializeUndeformedState()");
		m_gradientMatrix.clear();
		m_elementRestVolume.clear();
        for (int e = 0; e < m_elements.size(); e++) {
            const auto &element = m_elements[e];
            GradientMatrixType gradientMatrix;
			T elementRestVolume = 0;
            DiscretizationType::computeGradientMatrixAndRestVolume(element, m_X, gradientMatrix, elementRestVolume);
            m_gradientMatrix.push_back(gradientMatrix);
            m_elementRestVolume.push_back(elementRestVolume);
        }
#if 0
        dumper::writeElementsByte<4>(m_elements);
        dumper::writePositionsByte(m_X);
#endif
    }

	template<class dataType, int dim>
	void GridDeformerTet<std::vector<VECTOR<dataType, dim>>>::initializeAuxiliaryStructures()
	{
		// allocate memory for reshaped data
		int uncollisionSize = 0;
		int collisionSize = 0;
		for (const auto & f : m_elementFlags) {
			if (f == ElementFlag::unCollisionEl) uncollisionSize++;
			else if (f == ElementFlag::CollisionEl) collisionSize++;
			else throw std::logic_error("elements must be either unCollisionEl or CollisionEl");
		}

		m_nUncollisionBlocks = (uncollisionSize + (BlockWidth - 1)) / BlockWidth;
		m_nCollisionBlocks = (collisionSize + (BlockWidth - 1)) / BlockWidth;


#ifdef _WIN32
		m_reshapeUncollisionX = reinterpret_cast<BlockedShapeMatrixType>(_aligned_malloc(m_nUncollisionBlocks*BlockWidth*(d + 1)*d * sizeof(T), Alignment));
		m_reshapeCollisionX = reinterpret_cast<BlockedShapeMatrixType>(_aligned_malloc(m_nCollisionBlocks*BlockWidth*(d + 1)*d * sizeof(T), Alignment));

		m_reshapeUncollisionGradientMatrix = reinterpret_cast<BlockedMatrixType>(_aligned_malloc(m_nUncollisionBlocks*BlockWidth*d*d * sizeof(T), Alignment));
		m_reshapeCollisionGradientMatrix = reinterpret_cast<BlockedMatrixType>(_aligned_malloc(m_nCollisionBlocks*BlockWidth*d*d * sizeof(T), Alignment));

		m_reshapeUncollisionElementRestVolume = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nUncollisionBlocks*BlockWidth * sizeof(T), Alignment));
		m_reshapeCollisionElementRestVolume = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nCollisionBlocks*BlockWidth * sizeof(T), Alignment));


		m_reshapeUncollisionElementRestVolume = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nUncollisionBlocks * BlockWidth * sizeof(T), Alignment));
		m_reshapeCollisionElementRestVolume = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nCollisionBlocks * BlockWidth * sizeof(T), Alignment));

		m_reshapeUncollisionMuLow = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nUncollisionBlocks * BlockWidth * sizeof(T), Alignment));
		m_reshapeCollisionMuLow = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nCollisionBlocks * BlockWidth * sizeof(T), Alignment));

		m_reshapeUncollisionMuHigh = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nUncollisionBlocks * BlockWidth * sizeof(T), Alignment));
		m_reshapeCollisionMuHigh = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nCollisionBlocks * BlockWidth * sizeof(T), Alignment));

		m_reshapeUncollisionRangeMin = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nUncollisionBlocks * BlockWidth * sizeof(T), Alignment));
		m_reshapeCollisionRangeMin = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nCollisionBlocks * BlockWidth * sizeof(T), Alignment));

		m_reshapeUncollisionRangeMax = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nUncollisionBlocks * BlockWidth * sizeof(T), Alignment));
		m_reshapeCollisionRangeMax = reinterpret_cast<BlockedScalarType>(_aligned_malloc(m_nCollisionBlocks * BlockWidth * sizeof(T), Alignment));

#else

		m_reshapeUncollisionX = reinterpret_cast<BlockedShapeMatrixType>(aligned_alloc(Alignment, m_nUncollisionBlocks*BlockWidth*(d + 1)*d * sizeof(T)));
		m_reshapeCollisionX = reinterpret_cast<BlockedShapeMatrixType>(aligned_alloc(Alignment, m_nCollisionBlocks*BlockWidth*(d + 1)*d * sizeof(T)));

		m_reshapeUncollisionGradientMatrix = reinterpret_cast<BlockedMatrixType>(aligned_alloc(Alignment, m_nUncollisionBlocks*BlockWidth*d*d * sizeof(T)));
		m_reshapeCollisionGradientMatrix = reinterpret_cast<BlockedMatrixType>(aligned_alloc(Alignment, m_nCollisionBlocks*BlockWidth*d*d * sizeof(T)));

		m_reshapeUncollisionElementRestVolume = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nUncollisionBlocks*BlockWidth * sizeof(T)));
		m_reshapeCollisionElementRestVolume = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nCollisionBlocks*BlockWidth * sizeof(T)));

		m_reshapeUncollisionMuLow = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nUncollisionBlocks * BlockWidth * sizeof(T)));
		m_reshapeCollisionMuLow = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nCollisionBlocks * BlockWidth * sizeof(T)));

		m_reshapeUncollisionMuHigh = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nUncollisionBlocks * BlockWidth * sizeof(T)));
		m_reshapeCollisionMuHigh = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nCollisionBlocks * BlockWidth * sizeof(T)));

		m_reshapeUncollisionRangeMin = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nUncollisionBlocks * BlockWidth * sizeof(T)));
		m_reshapeCollisionRangeMin = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nCollisionBlocks * BlockWidth * sizeof(T)));

		m_reshapeUncollisionRangeMax = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nUncollisionBlocks * BlockWidth * sizeof(T)));
		m_reshapeCollisionRangeMax = reinterpret_cast<BlockedScalarType>(aligned_alloc(Alignment, m_nCollisionBlocks * BlockWidth * sizeof(T)));
#endif
		if (m_reshapeUncollisionX == nullptr || m_reshapeCollisionX == nullptr ||
			m_reshapeUncollisionGradientMatrix == nullptr || m_reshapeCollisionGradientMatrix == nullptr ||
			m_reshapeUncollisionElementRestVolume == nullptr || m_reshapeCollisionElementRestVolume == nullptr)
			throw std::logic_error("fail to allocate memory for m_reshapeX");

		// initialize reshaped data
		for (int e = 0, numOfUncollision = 0, numOfCollision = 0; e < m_elements.size(); e++) {
			if (m_elementFlags[e] == ElementFlag::unCollisionEl) {
				for (int i = 0; i < d + 1; i++)
					for (int j = 0; j < d; j++)
						m_reshapeUncollisionX[numOfUncollision / BlockWidth][i][j][numOfUncollision%BlockWidth] = m_X[m_elements[e][i]](j + 1);

				for (int i = 0; i < d; i++)
					for (int j = 0; j < d; j++)
						m_reshapeUncollisionGradientMatrix[numOfUncollision / BlockWidth][i + 3 * j][numOfUncollision%BlockWidth] = m_gradientMatrix[e](i + 1, j + 1);

				m_reshapeUncollisionElementRestVolume[numOfUncollision / BlockWidth][numOfUncollision%BlockWidth] = m_elementRestVolume[e];
				m_reshapeUncollisionMuLow[numOfUncollision / BlockWidth][numOfUncollision % BlockWidth] = m_muLow[e];
				m_reshapeUncollisionMuHigh[numOfUncollision / BlockWidth][numOfUncollision % BlockWidth] = m_muHigh[e];
				m_reshapeUncollisionRangeMin[numOfUncollision / BlockWidth][numOfUncollision % BlockWidth] = m_rangeMin[e];
				m_reshapeUncollisionRangeMax[numOfUncollision / BlockWidth][numOfUncollision % BlockWidth] = m_rangeMax[e];

				numOfUncollision++;
			}
			else if (m_elementFlags[e] == ElementFlag::CollisionEl) {
				for (int i = 0; i < d + 1; i++)
					for (int j = 0; j < d; j++)
						m_reshapeCollisionX[numOfCollision / BlockWidth][i][j][numOfCollision%BlockWidth] = m_X[m_elements[e][i]](j + 1);

				for (int i = 0; i < d; i++)
					for (int j = 0; j < d; j++)
						m_reshapeCollisionGradientMatrix[numOfCollision / BlockWidth][i + 3 * j][numOfCollision%BlockWidth] = m_gradientMatrix[e](i + 1, j + 1);

				m_reshapeCollisionElementRestVolume[numOfCollision / BlockWidth][numOfCollision % BlockWidth] = m_elementRestVolume[e];
				m_reshapeCollisionMuLow[numOfCollision / BlockWidth][numOfCollision % BlockWidth] = m_muLow[e];
				m_reshapeCollisionMuHigh[numOfCollision / BlockWidth][numOfCollision % BlockWidth] = m_muHigh[e];
				m_reshapeCollisionRangeMax[numOfCollision / BlockWidth][numOfCollision % BlockWidth] = m_rangeMax[e];
				m_reshapeCollisionRangeMin[numOfCollision / BlockWidth][numOfCollision % BlockWidth] = m_rangeMax[e];
				numOfCollision++;
			}
			else throw std::logic_error("elements must be either unCollisionEl or CollisionEl");
		}

		// initialize auxiliary structure
		std::vector<std::vector<int>> reshapeUncollisionIndices(m_X.size());
		std::vector<std::vector<int>> reshapeCollisionIndices(m_X.size());
		for (int e = 0, numOfUncollision = 0, numOfCollision = 0; e < m_elements.size(); e++) {
			if (m_elementFlags[e] == ElementFlag::unCollisionEl) {
				int blockIndex = numOfUncollision / BlockWidth;
				int blockOffset = numOfUncollision % BlockWidth;
				for (int v = 0; v < d + 1; v++) {
					int p = m_elements[e][v];
					int offset = (blockIndex * (d + 1) + v) * d * BlockWidth + blockOffset;
					reshapeUncollisionIndices[p].push_back(offset);
				}
				numOfUncollision++;
			}
			else if (m_elementFlags[e] == ElementFlag::CollisionEl) {
				int blockIndex = numOfCollision / BlockWidth;
				int blockOffset = numOfCollision % BlockWidth;
				for (int v = 0; v < d + 1; v++) {
					int p = m_elements[e][v];
					int offset = (blockIndex * (d + 1) + v) * d * BlockWidth + blockOffset;
					reshapeCollisionIndices[p].push_back(offset);
				}
				numOfCollision++;
			}
			else throw std::logic_error("elements must be either unCollisionEl or CollisionEl");
		}

		m_reshapeUncollisionIndicesOffsets.resize(m_X.size() + 1);
		m_reshapeCollisionIndicesOffsets.resize(m_X.size() + 1);
		m_reshapeUncollisionIndicesOffsets[0] = 0;
		m_reshapeCollisionIndicesOffsets[0] = 0;
		m_reshapeUncollisionIndicesValues.clear();
		m_reshapeCollisionIndicesValues.clear();
		for (size_t i = 0; i < m_X.size(); i++) {
			m_reshapeUncollisionIndicesOffsets[i + 1] = m_reshapeUncollisionIndicesOffsets[i] + (int)reshapeUncollisionIndices[i].size();
			m_reshapeCollisionIndicesOffsets[i + 1] = m_reshapeCollisionIndicesOffsets[i] + (int)reshapeCollisionIndices[i].size();
			for (const auto offset : reshapeUncollisionIndices[i])
				m_reshapeUncollisionIndicesValues.push_back(offset);
			for (const auto offset : reshapeCollisionIndices[i])
				m_reshapeCollisionIndicesValues.push_back(offset);
		}


#ifdef _WIN32
		m_reshapeUncollisionElement = reinterpret_cast<BlockedElementType>(_aligned_malloc(m_nUncollisionBlocks*BlockWidth*(d + 1) * sizeof(int), Alignment));
		m_reshapeCollisionElement = reinterpret_cast<BlockedElementType>(_aligned_malloc(m_nCollisionBlocks*BlockWidth*(d + 1) * sizeof(int), Alignment));
#else
		m_reshapeUncollisionElement = reinterpret_cast<BlockedElementType>(aligned_alloc(Alignment, m_nUncollisionBlocks*BlockWidth*(d + 1) * sizeof(int)));
		m_reshapeCollisionElement = reinterpret_cast<BlockedElementType>(aligned_alloc(Alignment, m_nCollisionBlocks*BlockWidth*(d + 1) * sizeof(int)));
#endif

		for (int b = 0; b < m_nUncollisionBlocks; b++)
			for (int v = 0; v < d + 1; v++)
				for (int e = 0; e < BlockWidth; e++)
					m_reshapeUncollisionElement[b][v][e] = 0;

		for (int b = 0; b < m_nCollisionBlocks; b++)
			for (int v = 0; v < d + 1; v++)
				for (int e = 0; e < BlockWidth; e++)
					m_reshapeCollisionElement[b][v][e] = 0;

		for (int e = 0, numOfUncollision = 0, numOfCollision = 0; e < m_elements.size(); e++) {
			if (m_elementFlags[e] == ElementFlag::unCollisionEl) {
				for (int v = 0; v < d + 1; v++)
					m_reshapeUncollisionElement[numOfUncollision / BlockWidth][v][numOfUncollision%BlockWidth] = m_elements[e][v];
				numOfUncollision++;
			}
			else if (m_elementFlags[e] == ElementFlag::CollisionEl) {
				for (int v = 0; v < d + 1; v++)
					m_reshapeCollisionElement[numOfCollision / BlockWidth][v][numOfCollision%BlockWidth] = m_elements[e][v];
				numOfCollision++;
			}
			else throw std::logic_error("elements must be either unCollisionEl or CollisionEl");
		}
	}

    template <class dataType, int dim>
    void GridDeformerTet<std::vector<VECTOR<dataType,dim>>>::updatePositionBasedState(const ElementFlag flag/*, const dataType rangeMin, const dataType rangeMax*/)
    {
        //LOG::SCOPE scope("GridDeformerTet::updatePositionBasedState()");
        /*
          ShapeMatrixType Ds;
          // #pragma omp parallel for
          for (int e = 0; e < m_elements.size(); e++) {
          if (m_elementFlags[e] == flag) {
          const auto &element = m_elements[e];
          computeShapeMatrix(Ds, element, m_X);
          auto F = computeGradient(Ds, m_gradientMatrix[e]);
          auto &U = m_U[e];
          auto &Sigma = m_Sigma[e];
          auto &V = m_V[e];

          F.Fast_Singular_Value_Decomposition(U, Sigma, V);

          for (int i = 0; i < d; i++) {
          if (Sigma(i + 1) < rangeMin)  Sigma(i + 1) = rangeMin;
          else if (Sigma(i + 1) > rangeMax)  Sigma(i + 1) = rangeMax;
          }
          }
          }
        */

        if (flag == ElementFlag::unCollisionEl) {
            blockX<T, BlockWidth>(&m_X[0](1), &m_reshapeUncollisionElement[0][0][0], m_nUncollisionBlocks, &m_reshapeUncollisionX[0][0][0][0]);
        }
        else if (flag == ElementFlag::CollisionEl) {
            blockX<T, BlockWidth>(&m_X[0](1), &m_reshapeCollisionElement[0][0][0], m_nCollisionBlocks, &m_reshapeCollisionX[0][0][0][0]);
        }
        else
            throw std::logic_error("Invalid Elementflag");
    }

    template <class dataType, int dim>
    void GridDeformerTet<std::vector<VECTOR<dataType,dim>>>::addCollisionForce(StateVariableType &f) const
    {
        {
            for (int c = 0; c < m_collisionConstraints.size(); c++) {
                const auto &constraint = m_collisionConstraints[c];
                VectorType x;
                x = DiscretizationType::interpolateX(constraint.m_elementIndex, constraint.m_weights, m_X);
                x -= constraint.m_xT;
                x *= -constraint.m_stiffness;

                DiscretizationType::distributeForces(x, constraint.m_elementIndex, constraint.m_weights, f);
            }
        }

        {
            for (int c = 0; c < m_collisionSutures.size(); c++) {
                const auto &suture = m_collisionSutures[c];
                VectorType x1;
                x1 = DiscretizationType::interpolateX(suture.m_elementIndex1, suture.m_weights1, m_X);
                auto x2 = DiscretizationType::interpolateX(suture.m_elementIndex2, suture.m_weights2, m_X);
//                x1 = x1 - x2 - (x1 - x2).Normalized()*suture.m_restLength;
				x1 = (x1 - x2).Projected(suture.m_normal);
                x1 *= -suture.m_stiffness;
                DiscretizationType::distributeForces(x1, suture.m_elementIndex1, suture.m_weights1, f);
                DiscretizationType::distributeForces(-x1, suture.m_elementIndex2, suture.m_weights2, f);
            }
        }
    }

    template <class dataType, int dim>
    void GridDeformerTet<std::vector<VECTOR<dataType,dim>>>:: addConstraintForce(StateVariableType &f) const
    {
        StateVariableType fLocal;
        {
            for (int c = 0; c < m_constraints.size(); c++) {
				const auto &constraint = m_constraints[c];
				if (constraint.m_stiffness) {
					VectorType x;
					x = DiscretizationType::interpolateX(constraint.m_elementIndex, constraint.m_weights, m_X);
					x -= constraint.m_xT;
					const T length = x.Lp_Norm(2);
					if (length > constraint.m_stressLimit)
						x *= constraint.m_stressLimit / length;
					x *= -constraint.m_stiffness;

					DiscretizationType::distributeForces(x, constraint.m_elementIndex, constraint.m_weights, f);
				}
            }

        }

        {
            for (int c = 0; c < m_sutures.size(); c++) {
                const auto &suture = m_sutures[c];
                VectorType x1;
                x1 = DiscretizationType::interpolateX(suture.m_elementIndex1, suture.m_weights1, m_X);
                auto x2 = DiscretizationType::interpolateX(suture.m_elementIndex2, suture.m_weights2, m_X);
                x1 = x1 - x2 - (x1 - x2).Normalized()*suture.m_restLength;
                x1 *= -suture.m_stiffness;
                DiscretizationType::distributeForces(x1, suture.m_elementIndex1, suture.m_weights1, f);
                DiscretizationType::distributeForces(-x1, suture.m_elementIndex2, suture.m_weights2, f);
            }
        }

		{
			for (int c = 0; c < m_fakeSutures.size(); c+=2) {
				VectorType x0, x1, x;
				const auto &constraint0 = m_fakeSutures[c];
				const auto &constraint1 = m_fakeSutures[c+1];

				x0 = DiscretizationType::interpolateX(constraint0.m_elementIndex, constraint0.m_weights, m_X);
				x1 = DiscretizationType::interpolateX(constraint1.m_elementIndex, constraint1.m_weights, m_X);
				x = (x1 - x0)/2;
				// x -= constraint.m_xT;
				x *= constraint0.m_stiffness;

				DiscretizationType::distributeForces(x, constraint0.m_elementIndex, constraint0.m_weights, f);
				DiscretizationType::distributeForces(-x, constraint1.m_elementIndex, constraint1.m_weights, f);
			}
		}
    }

	template<class dataType, int dim>
	void GridDeformerTet<std::vector<VECTOR<dataType, dim>>>::initializeCollisionElements()
	{
		for (int i = 0; i < m_elements.size(); i++) {
			const auto& e = DiscretizationType::getElementIndex(m_elements[i]);
			bool isR2 = true;
			for (const auto& idx : e)
				if (IteratorType::at(m_nodeType, idx) != NodeType::Collision) {
					isR2 = false;
					break;
				}
			if (isR2)
				m_elementFlags[i] = ElementFlag::CollisionEl;
		}
	}

    template <class dataType, int dim>
    void GridDeformerTet<std::vector<VECTOR<dataType,dim>>>::addElasticForce(std::vector<VECTOR<dataType, dim>> &SIMDf, const ElementFlag flag /*, const dataType rangeMin, const dataType rangeMax, const dataType weightProportion */ ) const
    {
        //alignas(sizeof(T)*BlockWidth) T muLow[BlockWidth];
        //alignas(sizeof(T)*BlockWidth) T muHigh[BlockWidth];
        //alignas(sizeof(T)*BlockWidth) T strainMin[BlockWidth];
        //alignas(sizeof(T)*BlockWidth) T strainMax[BlockWidth];

        // for (int i = 0; i < BlockWidth; i++) muLow[i] = weightProportion*weightProportion*m_uniformMu;
		//for (int i = 0; i < BlockWidth; i++) muHigh[i] = m_muHigh[0];//m_uniformMu;
        //for (int i = 0; i < BlockWidth; i++) strainMin[i] = rangeMin;
        //for (int i = 0; i < BlockWidth; i++) strainMax[i] = rangeMax;

        if (flag == ElementFlag::unCollisionEl) {
#ifdef _WIN32
			BlockedShapeMatrixType reshapeUncollisionf = reinterpret_cast<BlockedShapeMatrixType>(_aligned_malloc(m_nUncollisionBlocks*BlockWidth*(d + 1)*d * sizeof(T), Alignment));
#else
			BlockedShapeMatrixType reshapeUncollisionf = reinterpret_cast<BlockedShapeMatrixType>(aligned_alloc(Alignment, m_nUncollisionBlocks*BlockWidth*(d+1)*d*sizeof(T)));
#endif
			if (reshapeUncollisionf) {
				for (int b = 0; b < m_nUncollisionBlocks; b++)
					for (int v = 0; v < d + 1; v++)
						for (int i = 0; i < d; i++)
							for (int e = 0; e < BlockWidth; e++)
								reshapeUncollisionf[b][v][i][e] = 0;

#pragma omp parallel for
				for (int be = 0; be < m_nUncollisionBlocks; be++) {
					for (int ee = 0; ee < BlockWidth; ee += Tarch::Width)
						Add_Force<Tarch, T[BlockWidth]>(reinterpret_cast<T(&)[d + 1][d][BlockWidth]>(m_reshapeUncollisionX[be][0][0][ee]),
							reinterpret_cast<T(&)[d * d][BlockWidth]>(m_reshapeUncollisionGradientMatrix[be][0][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeUncollisionElementRestVolume[be][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeUncollisionMuLow[be][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeUncollisionMuHigh[be][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeUncollisionRangeMin[be][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeUncollisionRangeMax[be][ee]),
							reinterpret_cast<T(&)[d + 1][d][BlockWidth]>(reshapeUncollisionf[be][0][0][ee]));
				}

				unblockAddForce<T, BlockWidth>(&reshapeUncollisionf[0][0][0][0], &m_reshapeUncollisionIndicesOffsets[0], &m_reshapeUncollisionIndicesValues[0], (int)m_X.size(), &SIMDf[0](1));
#ifdef _WIN32
				_aligned_free(reshapeUncollisionf);
#else
				free(reshapeUncollisionf);
#endif
			}
			else {
				std::cerr << "Not enough space for reshapeUncollisionf" << std::endl;
				exit(1);
			}
        }
        else if (flag == ElementFlag::CollisionEl) {
#ifdef _WIN32
			BlockedShapeMatrixType reshapeCollisionf = reinterpret_cast<BlockedShapeMatrixType>(_aligned_malloc(m_nCollisionBlocks*BlockWidth*(d + 1)*d * sizeof(T), Alignment));
#else
			BlockedShapeMatrixType reshapeCollisionf = reinterpret_cast<BlockedShapeMatrixType>(aligned_alloc(Alignment, m_nCollisionBlocks*BlockWidth*(d+1)*d*sizeof(T)));
#endif
			if (reshapeCollisionf) {
				for (int b = 0; b < m_nCollisionBlocks; b++)
					for (int v = 0; v < d + 1; v++)
						for (int i = 0; i < d; i++)
							for (int e = 0; e < BlockWidth; e++)
								reshapeCollisionf[b][v][i][e] = 0;

#pragma omp parallel for
				for (int be = 0; be < m_nCollisionBlocks; be++) {
					for (int ee = 0; ee < BlockWidth; ee += Tarch::Width)
						Add_Force<Tarch, T[BlockWidth]>(reinterpret_cast<T(&)[d + 1][d][BlockWidth]>(m_reshapeCollisionX[be][0][0][ee]),
							reinterpret_cast<T(&)[d * d][BlockWidth]>(m_reshapeCollisionGradientMatrix[be][0][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeCollisionElementRestVolume[be][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeCollisionMuLow[be][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeCollisionMuHigh[be][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeCollisionRangeMin[be][ee]),
							reinterpret_cast<T(&)[BlockWidth]>(m_reshapeCollisionRangeMax[be][ee]),
							reinterpret_cast<T(&)[d + 1][d][BlockWidth]>(reshapeCollisionf[be][0][0][ee]));
				}

				unblockAddForce<T, BlockWidth>(&reshapeCollisionf[0][0][0][0], &m_reshapeCollisionIndicesOffsets[0], &m_reshapeCollisionIndicesValues[0], (int)m_X.size(), &SIMDf[0](1));
#ifdef _WIN32
				_aligned_free(reshapeCollisionf);
#else
				free(reshapeCollisionf);
#endif
			}
			else {
				std::cerr << "Not enough space for reshapeCollisionf" << std::endl;
				exit(1);
			}
        }
    }

    template <class dataType, int dim>
    void GridDeformerTet<std::vector<VECTOR<dataType,dim>>>::deallocateAuxiliaryStructures() {
#ifdef _WIN32
		if (m_reshapeUncollisionX) _aligned_free(m_reshapeUncollisionX);	
		if (m_reshapeCollisionX) _aligned_free(m_reshapeCollisionX);
		if (m_reshapeUncollisionGradientMatrix) _aligned_free(m_reshapeUncollisionGradientMatrix);
		if (m_reshapeCollisionGradientMatrix) _aligned_free(m_reshapeCollisionGradientMatrix);
		if (m_reshapeUncollisionElementRestVolume) _aligned_free(m_reshapeUncollisionElementRestVolume);
		if (m_reshapeCollisionElementRestVolume) _aligned_free(m_reshapeCollisionElementRestVolume);
		if (m_reshapeUncollisionElement) _aligned_free(m_reshapeUncollisionElement);
		if (m_reshapeCollisionElement) _aligned_free(m_reshapeCollisionElement);

		if (m_reshapeUncollisionMuLow) _aligned_free(m_reshapeUncollisionMuLow);
		if (m_reshapeCollisionMuLow) _aligned_free(m_reshapeCollisionMuLow);
		if (m_reshapeUncollisionMuHigh) _aligned_free(m_reshapeUncollisionMuHigh);
		if (m_reshapeCollisionMuHigh) _aligned_free(m_reshapeCollisionMuHigh);
		if (m_reshapeUncollisionRangeMin) _aligned_free(m_reshapeUncollisionRangeMin);
		if (m_reshapeCollisionRangeMin) _aligned_free(m_reshapeCollisionRangeMin);
		if (m_reshapeUncollisionRangeMax) _aligned_free(m_reshapeUncollisionRangeMax);
		if (m_reshapeCollisionRangeMax) _aligned_free(m_reshapeCollisionRangeMax);
#else
        free(m_reshapeUncollisionX);
        free(m_reshapeCollisionX);
        free(m_reshapeUncollisionGradientMatrix);
        free(m_reshapeCollisionGradientMatrix);
        free(m_reshapeUncollisionElementRestVolume);
        free(m_reshapeCollisionElementRestVolume);
        free(m_reshapeUncollisionElement);
        free(m_reshapeCollisionElement);

		free(m_reshapeUncollisionMuHigh);
		free(m_reshapeCollisionMuHigh);
		free(m_reshapeUncollisionMuLow);
		free(m_reshapeCollisionMuLow);
		free(m_reshapeUncollisionRangeMin);
		free(m_reshapeCollisionRangeMin);
		free(m_reshapeUncollisionRangeMax);
		free(m_reshapeCollisionRangeMax);

#endif
		m_reshapeUncollisionX = nullptr;
		m_reshapeCollisionX = nullptr;
		m_reshapeUncollisionGradientMatrix = nullptr;
		m_reshapeCollisionGradientMatrix = nullptr;
		m_reshapeUncollisionElementRestVolume = nullptr;
		m_reshapeCollisionElementRestVolume = nullptr;
		m_reshapeUncollisionElement = nullptr;
		m_reshapeCollisionElement = nullptr;

		m_reshapeUncollisionMuLow = nullptr;
		m_reshapeCollisionMuLow = nullptr;
		m_reshapeUncollisionMuHigh = nullptr;
		m_reshapeCollisionMuHigh = nullptr;
		m_reshapeUncollisionRangeMax = nullptr;
		m_reshapeCollisionRangeMax = nullptr;
		m_reshapeUncollisionRangeMin = nullptr;
		m_reshapeCollisionRangeMin = nullptr;


		
    }

/*void addElasticForce(StateVariableType &f, const ElementFlag flag, const T weightProportion) const {
            //LOG::SCOPE scope("GridDeformerTet::addElasticForce()");
            {
                for (int e = 0; e < m_elements.size(); e++)
                    if (m_elementFlags[e] == flag) {
                        const auto &element = m_elements[e];
                        ShapeMatrixType H;
			            ShapeMatrixType BioH;
                        ShapeMatrixType Ds;

                        const auto &gradientMatrix = m_gradientMatrix[e];
                        computeShapeMatrix(Ds, element, m_X);
                        auto F = computeGradient(Ds, gradientMatrix);

                        auto &U = m_U[e];
                        auto &V = m_V[e];
			            auto &Sigma = m_Sigma[e];

                        MatrixType R = U.Times_Transpose(V);
			            MatrixType BioR = U * (Sigma.Times_Transpose(V));
			            // computeDivergence need to a linear operation
                        MatrixType P = (F - R) * 2. * m_uniformMu * weightProportion * weightProportion + (F - BioR) * 2. * m_uniformMu;

                        const auto elementRestVolume = m_elementRestVolume[e];
                        H = computeDivergence(P, gradientMatrix, -elementRestVolume);
                        distributeForces(H, element, f);
                    }
            }
        }*/



}



// explicit template instantiation
namespace PhysBAM {
    template struct GridDeformerTet<std::vector<VECTOR<float, 3>>>;
}
