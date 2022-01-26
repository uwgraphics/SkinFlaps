#pragma once

#include <vector>

template<class T, int CoordinateStride>
void unblockAddForce(const T* fReshapedBasePtr, const int* reshapeIndicesOffsets, const int* reshapeIndicesValues, const int nParticles, T* f);

template<class T, int CoordinateStride>
void blockX(const T* X, const int* elementsPtr, const int nBlocks, T* XBasePtr);
