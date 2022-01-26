#include <omp.h>

#include "ReshapeDataStructure.h"

template<class T, int CoordinateStride>
void unblockAddForce(const T* fReshapedBasePtr, const int* reshapeIndicesOffsets, const int* reshapeIndicesValues, const int nParticles, T* f) {
    #pragma omp parallel for
    for (int i = 0; i < nParticles; i++) {
        T fX = 0., fY = 0., fZ = 0.;
        const int* offsetPtr = &reshapeIndicesValues[reshapeIndicesOffsets[i]];
        for (int j = reshapeIndicesOffsets[i]; j < reshapeIndicesOffsets[i+1]; j++, offsetPtr++){
            fX += fReshapedBasePtr[*offsetPtr];
            fY += fReshapedBasePtr[*offsetPtr + CoordinateStride];
            fZ += fReshapedBasePtr[*offsetPtr + 2 * CoordinateStride];
        }
        f[3*i] += fX;
        f[3*i + 1] += fY;
        f[3*i + 2] += fZ;
    }
}

template<class T, int CoordinateStride>
void blockX(const T* X, const int* elementsPtr, const int nBlocks, T* XBasePtr) {
    // Assume everything in elementsPtr can access valid location in X without seg fault
    using WideType = T (&)[CoordinateStride];
    constexpr int d = 3;
    auto reshapeElements = reinterpret_cast<const int (*) [d+1][CoordinateStride]>(elementsPtr);
    auto XReshapedBasePtr = reinterpret_cast<T (*) [d+1][d][CoordinateStride]>(XBasePtr);
    
    #pragma omp parallel for
    for( int b = 0; b < nBlocks; b++ )
        for ( int v = 0; v < d+1; v++) {
            WideType xBlock = XReshapedBasePtr[b][v][0];
            WideType yBlock = XReshapedBasePtr[b][v][1];
            WideType zBlock = XReshapedBasePtr[b][v][2];
            for ( int e = 0; e < CoordinateStride; e++) {
                int offset = reshapeElements[b][v][e];
                xBlock[e] = X[offset*d];
                yBlock[e] = X[offset*d + 1];
                zBlock[e] = X[offset*d + 2];
            }
        }
}

template
void unblockAddForce<float, 16>(const float* fReshapedBasePtr, const int* reshapeIndicesOffsets, const int* reshapeIndicesValues, const int nParticles, float* f);

template
void unblockAddForce<double, 16>(const double* fReshapedBasePtr, const int* reshapeIndicesOffsets, const int* reshapeIndicesValues, const int nParticles, double* f);

template
void blockX<float, 16>(const float* X, const int* elementsPtr, const int nBlocks, float* XBasePtr);

template
void blockX<double, 16>(const double* X, const int* elementsPtr, const int nBlocks, double* XBasePtr);
