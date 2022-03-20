#pragma once

#include <stdexcept>
#include <cusparse_v2.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>
#include <iostream>

// #include "dumper.h"

template <class DataType>
struct CudaSparse;

template<>
struct CudaSparse<float>
{
    // const:  W(sparse), S(sparse)
    // input:  D(sparse)
    // output: C(sparse)
    // tmp:    B = DW, WT

    using T = float;

    int m_nCollisionNodes;
    int m_nCollisionP;
    int m_nSelfCollisionP;
    int m_nnzSparse;
    int m_nnzB;
    int m_nnzC;

    // cusparse stuff
    cusparseHandle_t m_cusparse_handle;
    cusparseMatDescr_t m_descr;

    int    *m_rowPtrW,    *m_rowPtrS,    *m_rowPtrD,  *m_rowPtrWT;
    int *m_colIdxPtrW, *m_colIdxPtrS, *m_colIdxPtrD,  *m_colIdxPtrWT;
    T      *m_valPtrW,    *m_valPtrS,    *m_valPtrD, *m_valPtrC,    *m_valPtrB,    *m_valPtrWT;

    alignas (128) void *m_bufferB;
    alignas (128) void *m_bufferC;
    alignas (128) void* m_bufferConversion;


    csrgemm2Info_t m_infoB;
    csrgemm2Info_t m_infoC;

    int m_updateIdx;
    int m_updateSize;

    size_t m_bufferSizeB, m_bufferSizeC, m_bufferSizeConversion;

    const T ONE;

CudaSparse():
    m_cusparse_handle(0), m_descr(0), m_bufferB(nullptr), m_bufferC(nullptr), m_bufferConversion(nullptr), m_infoB(nullptr), m_infoC(nullptr), ONE(1), m_bufferSizeB(0), m_bufferSizeC(0), m_bufferSizeConversion(0),
	m_nCollisionNodes(0),
	m_nCollisionP(0),
	m_nSelfCollisionP(0),
	m_nnzSparse(0),
	m_nnzB(0),
	m_nnzC(0),
	m_rowPtrW(nullptr), m_rowPtrS(nullptr), m_rowPtrD(nullptr), m_rowPtrWT(nullptr),
	m_colIdxPtrW(nullptr), m_colIdxPtrS(nullptr), m_colIdxPtrD(nullptr), m_colIdxPtrWT(nullptr),
	m_valPtrW(nullptr), m_valPtrS(nullptr), m_valPtrD(nullptr), m_valPtrC(nullptr), m_valPtrB(nullptr), m_valPtrWT(nullptr),
	m_updateIdx(0), m_updateSize(0)
{	
        //cudaError_t cudaStat;
        cusparseStatus_t status;
        //cudaStream_t stream = NULL;

        status= cusparseCreate(&m_cusparse_handle);
        if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("CUSPARSE Library initialization failed");

        /* create and setup matrix descriptor */
        // sharing one descriptor for all matrices
        status = cusparseCreateMatDescr(&m_descr);
        if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Matrix descriptor initialization failed");

        cusparseSetMatType(m_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(m_descr,CUSPARSE_INDEX_BASE_ZERO);
    }

~CudaSparse() {
	cusparseStatus_t status;
	if (m_descr) {
		status = cusparseDestroyMatDescr(m_descr);
		// if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("cusparseDestroyMatDescr temp failed");
	}

	if (m_cusparse_handle) {
		status = cusparseDestroy(m_cusparse_handle);
		// if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("cusparseDestroy temp failed");
	}
}

    void initialize(const int m, const int pc, const int ps,
                    const int* rowPtrW_h, const int* colIdxPtrW_h, const T* valPtrW_h, const T* valPtrS_h) {
        // TODO: guard against double initialization
		// dumper::writeCSRbyte(pc + ps, rowPtrW_h, colIdxPtrW_h, valPtrW_h, m, "W_i.txt", "W_a.txt");
        m_nCollisionNodes = m;
        m_nCollisionP = pc;
        m_nSelfCollisionP = ps;
        m_nnzSparse = rowPtrW_h[pc+ps];

        m_updateIdx = rowPtrW_h[pc];
        m_updateSize = rowPtrW_h[pc+ps]-rowPtrW_h[pc];

        std::cout<<std::endl;
        std::cout<<"updateIdx = "<<m_updateIdx<<std::endl;
        std::cout<<"updateSize = "<<m_updateSize<<std::endl;

        cusparseStatus_t status;

        cudaError_t cudaStat;
        // spaces for const
        cudaStat = cudaMalloc((void**)&m_rowPtrW, (static_cast<size_t>(pc)+ps+1)*sizeof(int));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for rowPtrW");
        cudaStat = cudaMalloc((void**)&m_colIdxPtrW, m_nnzSparse*sizeof(int));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for colIdxPtrW");
        cudaStat = cudaMalloc((void**)&m_valPtrW, m_nnzSparse*sizeof(T));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for valPtrW");

        cudaStat = cudaMalloc((void**)&m_rowPtrS, (static_cast<size_t>(m)+1)*sizeof(int));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for rowPtrS");
        cudaStat = cudaMalloc((void**)&m_colIdxPtrS, static_cast<size_t>(m)*m*sizeof(int));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for colIdxPtrS");
        cudaStat = cudaMalloc((void**)&m_valPtrS, static_cast<size_t>(m)*m*sizeof(T));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for valPtrS");

        // space for inputs
        cudaStat = cudaMalloc((void**)&m_rowPtrD, (static_cast<size_t>(pc)+ps+1)*sizeof(int));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for rowPtrD");
        cudaStat = cudaMalloc((void**)&m_colIdxPtrD, (static_cast<size_t>(pc)+ps)*sizeof(int));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for colIdxPtrD");
        cudaStat = cudaMalloc((void**)&m_valPtrD, (static_cast<size_t>(pc)+ps)*sizeof(T));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for valPtrD");

        // space for temps
        cudaStat = cudaMalloc((void**)&m_rowPtrWT, (static_cast<size_t>(m)+1)*sizeof(int));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for rowPtrWT");
        cudaStat = cudaMalloc((void**)&m_colIdxPtrWT, m_nnzSparse*sizeof(int));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for colIdxPtrWT");
        cudaStat = cudaMalloc((void**)&m_valPtrWT, m_nnzSparse*sizeof(T));
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for valPtrWT");

        cudaStat = cudaMalloc((void**)&m_valPtrB, sizeof(T)*m_nnzSparse);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for valPtrB");

        m_nnzB = m_nnzSparse;
        m_nnzC = m_nnzSparse;

        int  *rowPtrS_h,    *rowPtrD_h;
        int *colIdxPtrS_h, *colIdxPtrD_h;

        rowPtrS_h = new int[m+1];
        colIdxPtrS_h = new int[m*m];

        rowPtrD_h = new int[pc+ps+1];
        colIdxPtrD_h = new int[pc+ps];

        rowPtrS_h[0] = 0;
        for (int i=0; i<m; i++) {
            rowPtrS_h[i+1] = (i+1)*m;
            for (int jj=rowPtrS_h[i]; jj<rowPtrS_h[i+1]; jj++)
                colIdxPtrS_h[jj] = jj-rowPtrS_h[i];
        }

        rowPtrD_h[0] = 0;
        for (int i=0; i<pc+ps; i++) {
            rowPtrD_h[i+1] = i+1;
            colIdxPtrD_h[i] = i;
        }


        cudaStat = cudaMemcpy(m_rowPtrD, rowPtrD_h, (static_cast<size_t>(pc)+ps+1)*sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for rowPtrD");
        cudaStat = cudaMemcpy(m_colIdxPtrD, colIdxPtrD_h, (static_cast<size_t>(pc)+ps)*sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for colIdxPtrD");

        // copy data in
        cudaStat = cudaMemcpy(m_rowPtrW, rowPtrW_h, (static_cast<size_t>(pc)+ps+1)*sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for rowPtrW");
        cudaStat = cudaMemcpy(m_colIdxPtrW, colIdxPtrW_h, (size_t)(m_nnzSparse*sizeof(int)), cudaMemcpyHostToDevice);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for colIdxPtrW");
        cudaStat = cudaMemcpy(m_valPtrW, valPtrW_h, (size_t)(m_nnzSparse*sizeof(T)), cudaMemcpyHostToDevice);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for valPtrW");

        cudaStat = cudaMemcpy(m_rowPtrS, rowPtrS_h, (static_cast<size_t>(m)+1)*sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for rowPtrS");
        cudaStat = cudaMemcpy(m_colIdxPtrS, colIdxPtrS_h, static_cast<size_t>(m)*m*sizeof(int), cudaMemcpyHostToDevice);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for colIdxPtrS");
        cudaStat = cudaMemcpy(m_valPtrS, valPtrS_h, static_cast<size_t>(m)*m*sizeof(T), cudaMemcpyHostToDevice);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for valPtrS");

        delete[]    rowPtrS_h; delete[]    rowPtrD_h;
        delete[] colIdxPtrS_h; delete[] colIdxPtrD_h;
        // cudaStream_t stream = NULL;
        // status = cusparseGetStream(m_cusparse_handle, &stream);
        // if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Get Stream failed");
        // 
 #if TIMING
        // cudaEvent_t start;
        // cudaEvent_t stop;
        // cudaEventCreate(&start);
        // cudaEventCreate(&stop);
        // float msecTotal = 0.0f;
#endif

        cusparseSetPointerMode(m_cusparse_handle, CUSPARSE_POINTER_MODE_HOST);

        status = cusparseCreateCsrgemm2Info(&m_infoB);
        if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("CreateCsrgemm2Info B failed");

        status = cusparseCreateCsrgemm2Info(&m_infoC);
        if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("CreateCsrgemm2Info C failed");

#if TIMING
        cudaEventRecord(start, stream);
#endif
        size_t bufferSizeConversion;

        status =
            cusparseCsr2cscEx2_bufferSize(m_cusparse_handle, pc + ps, m, m_nnzSparse,
                m_valPtrW, m_rowPtrW, m_colIdxPtrW,
                m_valPtrWT, m_rowPtrWT, m_colIdxPtrWT,
                CUDA_R_32F, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1,
                &bufferSizeConversion);
        if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Compute buffersize for conversion failed");

        if (!m_bufferConversion || bufferSizeConversion > m_bufferSizeConversion) {

            if (m_bufferConversion) {
                cudaStat = cudaFree(m_bufferConversion);
                if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for bufferConversion");
                m_bufferConversion = nullptr;
            }
            cudaStat = cudaMalloc(&m_bufferConversion, bufferSizeConversion);
            if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for temp buffer");
            m_bufferSizeConversion = bufferSizeConversion;
            std::cout << "computed buffersizeConversion = " << bufferSizeConversion << std::endl;
        }

        status =
            cusparseCsr2cscEx2
    (m_cusparse_handle, pc + ps, m, m_nnzSparse,
                m_valPtrW, m_rowPtrW, m_colIdxPtrW,
                m_valPtrWT, m_rowPtrWT, m_colIdxPtrWT,
                CUDA_R_32F, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1,
                m_bufferConversion);

#if TIMING
         cudaEventRecord(stop, stream);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&msecTotal, start, stop);
#endif
         cudaStat = cudaDeviceSynchronize();
         if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Conversion from CSR to CSC format failed");
         if (cudaStat != cudaSuccess) throw std::logic_error("couldn't synchronized device");

    }

    void reInitialize(const int m, const T* valPtrS_h) {
        cudaError_t cudaStat = cudaMemcpy(m_valPtrS, valPtrS_h, sizeof(T)*m*m, cudaMemcpyHostToDevice);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for valPtrS");
    }

    void rankKUpdate(const float* valPtrD_h, const int* colIdxPtrW_h, const float* valPtrW_h) {
        cudaEvent_t start;
        cudaEvent_t stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        float msecTotal = 0.0f;

        cudaError_t cudaStat;
        cusparseStatus_t status;

        cudaStream_t stream = NULL;
         status = cusparseGetStream(m_cusparse_handle, &stream);
         if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Get Stream failed");

         int m = m_nCollisionNodes;
         int pc = m_nCollisionP;
         int ps = m_nSelfCollisionP;

        status= cusparseCreateMatDescr(&m_descr);
        if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Matrix descriptor initialization failed");
        // update W
        if (ps) {
#if TIMING
            cudaEventRecord(start);
#endif
            cudaStat = cudaMemcpy(m_valPtrW+m_updateIdx, valPtrW_h+m_updateIdx, (size_t)(m_updateSize*sizeof(T)), cudaMemcpyHostToDevice);
            if (cudaStat != cudaSuccess) throw std::logic_error("couldn't update mem for valPtrW");

            cudaStat = cudaMemcpy(m_colIdxPtrW+m_updateIdx, colIdxPtrW_h+m_updateIdx, (size_t)(m_updateSize*sizeof(T)), cudaMemcpyHostToDevice);
#if TIMING
            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            cudaEventElapsedTime(&msecTotal, start, stop);
#endif
            if (cudaStat != cudaSuccess) throw std::logic_error("couldn't update mem for colIdxPtrW_h");
#if TIMING
            std::cout<<"            Copy W                  Time : "<<msecTotal/1000<<std::endl;

// create WT
            cudaEventRecord(start, stream);
#endif
            size_t bufferSizeConversion;

            status =
                cusparseCsr2cscEx2_bufferSize(m_cusparse_handle, pc + ps, m, m_nnzSparse,
                    m_valPtrW, m_rowPtrW, m_colIdxPtrW,
                    m_valPtrWT, m_rowPtrWT, m_colIdxPtrWT,
                    CUDA_R_32F, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1,
                    &bufferSizeConversion);
            if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Compute buffersize for conversion failed");

            if (!m_bufferConversion || bufferSizeConversion > m_bufferSizeConversion) {

                if (m_bufferConversion) {
                    cudaStat = cudaFree(m_bufferConversion);
                    if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for bufferConversion");
                    m_bufferConversion = nullptr;
                }
                cudaStat = cudaMalloc(&m_bufferConversion, bufferSizeConversion);
                if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for temp buffer");
                m_bufferSizeConversion = bufferSizeConversion;
                std::cout << "computed buffersizeConversion = " << bufferSizeConversion << std::endl;
                }

            status =
                cusparseCsr2cscEx2
                (m_cusparse_handle, pc + ps, m, m_nnzSparse,
                    m_valPtrW, m_rowPtrW, m_colIdxPtrW,
                    m_valPtrWT, m_rowPtrWT, m_colIdxPtrWT,
                    CUDA_R_32F, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1,
                    m_bufferConversion);

#if TIMING
            cudaEventRecord(stop, stream);
            cudaEventSynchronize(stop);
            cudaEventElapsedTime(&msecTotal, start, stop);
#endif
            cudaStat = cudaDeviceSynchronize();
            if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Conversion from CSR to CSC format failed");
            if (cudaStat != cudaSuccess) throw std::logic_error("couldn't synchronized device");
#if TIMING
            std::cout << "            Compute WT              Time : " << msecTotal / 1000 << std::endl;
#endif
        }

#if TIMING
        cudaEventRecord(start);
#endif
        cudaStat = cudaMemcpy(m_valPtrD, valPtrD_h, (static_cast<size_t>(pc)+ps)*sizeof(T), cudaMemcpyHostToDevice);
#if TIMING
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&msecTotal, start, stop);
        std::cout<<"            Copyin D                Time : "<<msecTotal/1000<<std::endl;
#endif
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for valPtrD");



        // compute buffer size
        cusparseStatus_t statusB;
        cusparseStatus_t statusC;
#if TIMING
         // cudaEventRecord(start, stream);
#endif
        size_t bufferSizeB, bufferSizeC;
        statusB = cusparseScsrgemm2_bufferSizeExt(m_cusparse_handle, pc+ps, m, pc+ps, &ONE,
                                                  m_descr, pc+ps, m_rowPtrD, m_colIdxPtrD,
                                                  m_descr, m_nnzSparse, m_rowPtrW, m_colIdxPtrW,
                                                  NULL,
                                                  m_descr, 0, NULL, NULL,
                                                  m_infoB,
                                                  &bufferSizeB);

        if (statusB != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Compute buffersize B failed");

        if (!m_bufferB || bufferSizeB > m_bufferSizeB) {

            if (m_bufferB) {
                cudaStat = cudaFree(m_bufferB);
                if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for bufferB");
				m_bufferB = nullptr;
            }
            cudaStat = cudaMalloc(&m_bufferB, bufferSizeB);
            if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for temp buffer");
            m_bufferSizeB = bufferSizeB;
            std::cout<<"computed buffersizeB = "<<bufferSizeB<<std::endl;
        }

#if 0
		{
		// copy out stuff for debug
			T* valPtrD_d = new T[pc + ps];
			int* rowPtrD_d = new int[pc + ps + 1];
			int * colIdxPtrD_d = new int[pc + ps];

			T* valPtrW_d = new T[m_nnzSparse];
			int* rowPtrW_d = new int[pc + ps + 1];
			int* colIdxPtrW_d = new int[m_nnzSparse];

			cudaStat = cudaMemcpy(rowPtrD_d,m_rowPtrD, (size_t)((pc + ps + 1) * sizeof(int)), cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for rowPtrD");
			cudaStat = cudaMemcpy(colIdxPtrD_d,m_colIdxPtrD, (size_t)((pc + ps) * sizeof(int)), cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for colIdxPtrD");
			cudaStat = cudaMemcpy(valPtrD_d, m_valPtrD,  (size_t)((pc + ps) * sizeof(T)), cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for valPtrD");

			// copy data in
			cudaStat = cudaMemcpy(rowPtrW_d, m_rowPtrW,  (size_t)((pc + ps + 1) * sizeof(int)), cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for rowPtrW");
			cudaStat = cudaMemcpy(colIdxPtrW_d, m_colIdxPtrW,  (size_t)(m_nnzSparse * sizeof(int)), cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for colIdxPtrW");
			cudaStat = cudaMemcpy(valPtrW_d, m_valPtrW,  (size_t)(m_nnzSparse * sizeof(T)), cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for valPtrW");

			dumper::writeCSRbyte(pc + ps, rowPtrW_d, colIdxPtrW_d, valPtrW_d, m, "Wd_i.txt", "Wd_a.txt");
			dumper::writeCSRbyte(pc + ps, rowPtrD_d, colIdxPtrD_d, valPtrD_d, pc+ps, "Dd_i.txt", "Dd_a.txt");

			delete[] rowPtrW_d;
			delete[] colIdxPtrW_d;
			delete[] valPtrW_d;

			delete[] rowPtrD_d;
			delete[] colIdxPtrD_d;
			delete[] valPtrD_d;

		}
#endif

#if TIMING
        cudaEventRecord(start, stream);
#endif
        statusB = cusparseScsrgemm2(m_cusparse_handle, pc+ps, m, pc+ps, &ONE,
                                    m_descr, pc+ps, m_valPtrD, m_rowPtrD, m_colIdxPtrD,
                                    m_descr, m_nnzSparse, m_valPtrW, m_rowPtrW, m_colIdxPtrW,
                                    NULL,
                                    m_descr, 0, NULL, NULL, NULL,
                                    m_descr, m_valPtrB, m_rowPtrW, m_colIdxPtrW,
                                    m_infoB, m_bufferB);

        if (statusB != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Compute csrgemm2 B failed");
        statusC = cusparseScsrgemm2_bufferSizeExt(m_cusparse_handle, m, m, pc+ps, &ONE,
                                                  m_descr, m_nnzSparse, m_rowPtrWT, m_colIdxPtrWT,
                                                  m_descr, m_nnzSparse, m_rowPtrW, m_colIdxPtrW,
                                                  &ONE,
                                                  m_descr, m*m, m_rowPtrS, m_colIdxPtrS,
                                                  m_infoC,
                                                  &bufferSizeC);

        if (statusC != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Compute buffersize failed");
        if (!m_bufferC || bufferSizeC > m_bufferSizeC) {
            if (m_bufferC) {
                cudaStat = cudaFree(m_bufferC);
                if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for bufferB");
				m_bufferC = nullptr;
            }
            cudaStat = cudaMalloc(&m_bufferC, bufferSizeC);
            if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for temp bufferC");
            m_bufferSizeC = bufferSizeC;
            std::cout<<"computed buffersizeC = "<<bufferSizeC<<std::endl;
        }

         statusC = cusparseScsrgemm2(m_cusparse_handle, m, m, pc+ps, &ONE,
                                    m_descr, m_nnzSparse, m_valPtrWT, m_rowPtrWT, m_colIdxPtrWT,
                                    m_descr, m_nnzSparse, m_valPtrB, m_rowPtrW, m_colIdxPtrW,
                                    &ONE,
                                    m_descr, m_nCollisionNodes*m_nCollisionNodes, m_valPtrS, m_rowPtrS, m_colIdxPtrS,
                                    m_descr, m_valPtrC, m_rowPtrS, m_colIdxPtrS,
                                    m_infoC, m_bufferC);

#if TIMING
         cudaEventRecord(stop,stream);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&msecTotal, start, stop);
         std::cout<<"            Rank K Update           Time : "<<msecTotal/1000<<std::endl;
#endif

        if (statusC != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("Compute csrgemm2 C failed");

    }

	void release() {
		cudaError_t cudaStat;
		cusparseStatus_t status;

		if (m_bufferB) {
			cudaStat = cudaFree(m_bufferB);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for bufferB");
			m_bufferB = nullptr;
		}

		m_bufferSizeB = 0;

		if (m_bufferC) {
			cudaStat = cudaFree(m_bufferC);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for bufferC");
			m_bufferC = nullptr;
		}

		m_bufferSizeC = 0;


        if (m_bufferConversion) {
            cudaStat = cudaFree(m_bufferConversion);
            if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for bufferC");
            m_bufferConversion = nullptr;
        }

        m_bufferSizeC = 0;

		if (m_infoB) {
			status = cusparseDestroyCsrgemm2Info(m_infoB);
			if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("DestroyCsrgemm2Info B failed");
			m_infoB = nullptr;
		}

		if (m_infoC) {
			status = cusparseDestroyCsrgemm2Info(m_infoC);
			if (status != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("DestroyCsrgemm2Info C failed");
			m_infoC = nullptr;
		}

		if (m_rowPtrW) {
			cudaStat = cudaFree(m_rowPtrW);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for rowPtrW");
			m_rowPtrW = nullptr;
		}
		if (m_rowPtrS) {
			cudaStat = cudaFree(m_rowPtrS);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for rowPtrS");
			m_rowPtrS = nullptr;
		}
		if (m_rowPtrD) {
			cudaStat = cudaFree(m_rowPtrD);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for rowPtrD");
			m_rowPtrD = nullptr;
		}
		if (m_rowPtrWT) {
			cudaStat = cudaFree(m_rowPtrWT);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for rowPtrWT");
			m_rowPtrWT = nullptr;
		}
		if (m_colIdxPtrW) {
			cudaStat = cudaFree(m_colIdxPtrW);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for colIdxPtrW");
			m_colIdxPtrW = nullptr;
		}
		if (m_colIdxPtrS) {
			cudaStat = cudaFree(m_colIdxPtrS);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for colIdxPtrS");
			m_colIdxPtrS = nullptr;
		}
		if (m_colIdxPtrD) {
			cudaStat = cudaFree(m_colIdxPtrD);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for colIdxPtrD");
			m_colIdxPtrD = nullptr;
		}
		if (m_colIdxPtrWT) {
			cudaStat = cudaFree(m_colIdxPtrWT);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for colIdxPtrWT");
			m_colIdxPtrWT = nullptr;
		}
		if (m_valPtrW) {
			cudaStat = cudaFree(m_valPtrW);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for valPtrW");
			m_valPtrW = nullptr;
		}
		if (m_valPtrS) {
			cudaStat = cudaFree(m_valPtrS);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for valPtrS");
			m_valPtrS = nullptr;
		}
		if (m_valPtrD) {
			cudaStat = cudaFree(m_valPtrD);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for valPtrD");
			m_valPtrD = nullptr;
		}
		if (m_valPtrB) {
			cudaStat = cudaFree(m_valPtrB);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for valPtrB");
			m_valPtrB = nullptr;
		}
		if (m_valPtrWT) {
			cudaStat = cudaFree(m_valPtrWT);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for valPtrWT");
			m_valPtrWT = nullptr;
		}

       
    }
};


template <class DataType>
struct CudaSolve;

template<>
struct CudaSolve<float> {
    using T = float;

    int m_matSize;

    cusolverDnHandle_t  m_cusolver_handle;
    const cublasFillMode_t UPLO;

    T* m_ptrA;
    T* m_ptrB;
    int* m_info;

    CudaSolve() : m_cusolver_handle(0), UPLO(CUBLAS_FILL_MODE_LOWER), m_matSize(0), m_ptrA(nullptr), m_ptrB(nullptr), m_info(nullptr) {
        cudaError_t cudaStat = cudaSuccess;
        cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
        //cudaStream_t stream = NULL;
        // create cusolver handle, bind a stream
        status = cusolverDnCreate(&m_cusolver_handle);
        if(CUSOLVER_STATUS_SUCCESS != status) throw std::logic_error("failed cusolverDnCreate");

        //cudaStat = cudaStreamCreateWithFlags(&stream, cudaStreamDefault);
        //if (cudaSuccess != cudaStat) throw std::logic_error("failed cudaStreamCreateWithFlags");
    }

	~CudaSolve() {
		cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
		status = cusolverDnDestroy(m_cusolver_handle);
		// if (CUSOLVER_STATUS_SUCCESS != status) throw std::logic_error("failed cusolverDnDestroy");
	}

    void initialize(const int m) {
        m_matSize = m;
        cudaError_t cudaStat = cudaSuccess;
        cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;

        // cudaEvent_t start, stop;
        // float msecTotal = 0.0f;
        // cudaEventCreate(&start);
        // cudaEventCreate(&stop);

        // cudaStream_t stream = NULL;
        // status = cusolverDnGetStream(m_cusolver_handle, &stream);
        // if(CUSOLVER_STATUS_SUCCESS != status) throw std::logic_error("failed cusolverDnGetStream");

        cudaStat = cudaMalloc ((void**)&m_ptrB, sizeof(T) * m );
        if (cudaSuccess != cudaStat) throw std::logic_error("failed cudaMalloc ptrB");

        cudaStat = cudaMalloc ((void**)&m_info, sizeof(int));
        if (cudaSuccess != cudaStat) throw std::logic_error("failed cudaMalloc info");
    }

    void factorize() {
        cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
        cudaError_t cudaStat = cudaSuccess;

        int workspaceSize;
        T* workspace{};

#if TIMING
        cudaEvent_t start, stop;
        float msecTotal = 0.0f;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
#endif

        cudaStream_t stream = NULL;
        status = cusolverDnGetStream(m_cusolver_handle, &stream);
        if(CUSOLVER_STATUS_SUCCESS != status) throw std::logic_error("failed cusolverDnGetStream");
        status = cusolverDnSpotrf_bufferSize(m_cusolver_handle,
                                             UPLO,
                                             m_matSize,
                                             m_ptrA,
                                             m_matSize,
                                             &workspaceSize);
        if (CUSOLVER_STATUS_SUCCESS != status) throw std::logic_error("failed cusolverDnDpotrf_bufferSize");

        cudaStat = cudaMalloc ((void**)&workspace, workspaceSize*sizeof(T));
        if (cudaSuccess != cudaStat) throw std::logic_error("failed cudaMalloc Workspace");

#if TIMING
        cudaEventRecord(start, stream);
        #endif

        status = cusolverDnSpotrf(
            m_cusolver_handle,
            UPLO,
            m_matSize,
            m_ptrA,
            m_matSize,
            workspace,
            workspaceSize,
            m_info
            );
#if TIMING
        cudaEventRecord(stop, stream);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&msecTotal, start, stop);
#endif

        if (CUSOLVER_STATUS_SUCCESS != status) throw std::logic_error("failed cusolverDnDpotrf");
#if TIMING
        std::cout<<"            Dense Fact              Time : "<<msecTotal/1000<<" s"<<std::endl;
#endif

        cudaStat = cudaFree(workspace);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for workspace");
    }

    void solve(// const T* rhs_h
        ) {
        cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
        cudaError_t cudaStat = cudaSuccess;

#if TIMING
        cudaEvent_t start, stop;
        float msecTotal = 0.0f;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
#endif

        cudaStream_t stream = NULL;
        status = cusolverDnGetStream(m_cusolver_handle, &stream);
        if(CUSOLVER_STATUS_SUCCESS != status) throw std::logic_error("failed cusolverDnGetStream");

#if TIMING
        cudaEventRecord(start, stream);
#endif
        status = cusolverDnSpotrs(
            m_cusolver_handle,
            UPLO,
            m_matSize,
            1, // solving one dimension
            m_ptrA,
            m_matSize,
            m_ptrB,
            m_matSize,
            m_info);

#if TIMING
        cudaEventRecord(stop, stream);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&msecTotal, start, stop);
#endif

        if (CUSOLVER_STATUS_SUCCESS != status) throw std::logic_error("failed cusolverDnDpotrs");

#if TIMING
        std::cout<<"            Solve                   Time : "<<msecTotal/1000<<" s"<<std::endl;
#endif
    }

    void release() {
        cudaError_t cudaStat = cudaSuccess;

		if (m_ptrB) {
			cudaStat = cudaFree(m_ptrB);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for ptrB");
			m_ptrB = nullptr;
		}
		if (m_info) {
			cudaStat = cudaFree(m_info);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for info");
			m_info = nullptr;
		}
    }
};

template <class DataType>
class CudaWrapper;

template <>
class CudaWrapper<float>
{
    using T = float;
    struct CudaSparse<float> m_sparse;
    struct CudaSolve<float> m_solve;

public:
	T *m_valPtrC = nullptr;
	cudaStream_t stream = NULL;

	CudaWrapper() {
		cudaError_t cudaStat;
		cusparseStatus_t statusS;
		cusolverStatus_t statusD = CUSOLVER_STATUS_SUCCESS;

		cudaStat = cudaStreamCreateWithFlags(&stream, cudaStreamDefault);
		if (cudaSuccess != cudaStat) throw std::logic_error("failed cudaStreamCreateWithFlags");

		statusS = cusparseSetStream(m_sparse.m_cusparse_handle, stream);
		if (statusS != CUSPARSE_STATUS_SUCCESS) throw std::logic_error("failed cusparseSetStream");
		statusD = cusolverDnSetStream(m_solve.m_cusolver_handle, stream);
		if (CUSOLVER_STATUS_SUCCESS != statusD) throw std::logic_error("failed cusolverDnSetStream");

	}
	~CudaWrapper() {
		cudaError_t cudaStat;

#if 1
		if (stream) {
			cudaStat = cudaStreamDestroy(stream);
			// if (cudaStat != cudaSuccess) throw std::logic_error("couldn't destroy stream");
		}
#endif
	}

    inline void initialize(const int m, const int pc, const int ps,
                           const int* rowPtrW_h, const int* colIdxPtrW_h, const T* valPtrW_h, const T* valPtrS_h) {

        cudaError_t cudaStat;
        //cusparseStatus_t statusS;
        //cusolverStatus_t statusD = CUSOLVER_STATUS_SUCCESS;
        // cudaStream_t stream = NULL;
        // space for output
        cudaStat = cudaMalloc((void**)&m_valPtrC, sizeof(T)*m*m);
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't allocate mem for valPtrC");
        m_sparse.initialize(m, pc, ps, rowPtrW_h, colIdxPtrW_h, valPtrW_h, valPtrS_h);
        m_sparse.m_valPtrC = m_valPtrC;
        m_solve.initialize(m);
        m_solve.m_ptrA = m_valPtrC;

      

    }

    inline void reInitialize(const int m, const T* valPtrS_h) {
        m_sparse.reInitialize(m, valPtrS_h);
    }


    inline void rankKUpdate(const float* valPtrD_h, const int* colIdxPtrW_h, const float* valPtrW_h) {
        m_sparse.rankKUpdate(valPtrD_h, colIdxPtrW_h, valPtrW_h);
    }

    inline void release() {
		cudaError_t cudaStat;

		m_solve.release();
        m_sparse.release();

		if (m_valPtrC) {
			cudaStat = cudaFree(m_valPtrC);
			if (cudaStat != cudaSuccess) throw std::logic_error("couldn't free mem for valPtrC");
			m_valPtrC = nullptr;
		}
    }

    inline void factorize() {
        m_solve.factorize();
    }

    inline void solve(T* x_h, const T* rhs_h) {
        cudaError_t cudaStat;
#if TIMING
        cudaEvent_t start, stop;
        float msecTotal = 0.0f;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);

        cudaEventRecord(start);
#endif
        cudaStat = cudaMemcpy(m_solve.m_ptrB, rhs_h, sizeof(T) * m_solve.m_matSize, cudaMemcpyHostToDevice);
#if TIMING
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&msecTotal, start, stop);
#endif
        if (cudaSuccess != cudaStat)
            throw std::logic_error("failed cudaMemcpy ptrB");
#if TIMING
        std::cout<<"            Copy In rhs             Time : "<<msecTotal/1000<<" s"<<std::endl;
#endif

        m_solve.solve();
#if TIMING
        cudaEventRecord(start);
#endif
        cudaStat = cudaMemcpy(x_h, m_solve.m_ptrB, (size_t)(m_solve.m_matSize*sizeof(T)), cudaMemcpyDeviceToHost);
#if TIMING
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&msecTotal, start, stop);
#endif
        if (cudaStat != cudaSuccess) throw std::logic_error("couldn't copy mem for x_h");
#if TIMING
        std::cout<<"            Copy Out rhs            Time : "<<msecTotal/1000<<" s"<<std::endl;
#endif

    }
};
