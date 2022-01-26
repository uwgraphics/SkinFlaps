#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Matrix_Times_Transpose.h"

#include <Thread_Queueing/PTHREAD_QUEUE.h>
#include <Kernel_Serial_Base_Helper.h>

using namespace SIMD_Numeric_Kernel;

template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}

struct timeval starttime, stoptime;
void
start_timer ()
{
  gettimeofday (&starttime, NULL);
}

void
stop_timer ()
{
  gettimeofday (&stoptime, NULL);
}

double
get_time ()
{
  return (double) stoptime.tv_sec - (double) starttime.tv_sec +
    (double) 1e-6 *(double) stoptime.tv_usec -
    (double) 1e-6 *(double) starttime.tv_usec;
}


template < int SIZE > class Matrix_Times_Transpose_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_B;
  float *_local_C;

public:
    explicit Matrix_Times_Transpose_SCALAR_None (float *A_in, float *B_in,
                                                 float *C_in):_local_A (A_in),
    _local_B (B_in), _local_C (C_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][9][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rB = reinterpret_cast < fullArray2 > (*_local_B);
    fullArray3 _rC = reinterpret_cast < fullArray3 > (*_local_C);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[9][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][0][chunk_offset]);
        refArray2 Bk =
          reinterpret_cast < refArray2 > (_rB[index][0][chunk_offset]);
        refArray3 Ck =
          reinterpret_cast < refArray3 > (_rC[index][0][chunk_offset]);

        Matrix_Times_Transpose < SIMDArchitectureScalar<float>, float[16] > (Ak, Bk, Ck);
      }

  }
};


#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Matrix_Times_Transpose_AVX_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_B;
  float *_local_C;

public:
    explicit Matrix_Times_Transpose_AVX_None (float *A_in, float *B_in,
                                              float *C_in):_local_A (A_in),
    _local_B (B_in), _local_C (C_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][9][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rB = reinterpret_cast < fullArray2 > (*_local_B);
    fullArray3 _rC = reinterpret_cast < fullArray3 > (*_local_C);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[9][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][0][chunk_offset]);
        refArray2 Bk =
          reinterpret_cast < refArray2 > (_rB[index][0][chunk_offset]);
        refArray3 Ck =
          reinterpret_cast < refArray3 > (_rC[index][0][chunk_offset]);

        Matrix_Times_Transpose < SIMDArchitectureAVX2<float>, float[16] > (Ak, Bk, Ck);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Matrix_Times_Transpose_MIC_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_B;
  float *_local_C;

public:
    explicit Matrix_Times_Transpose_MIC_None (float *A_in, float *B_in,
                                              float *C_in):_local_A (A_in),
    _local_B (B_in), _local_C (C_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][9][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rB = reinterpret_cast < fullArray2 > (*_local_B);
    fullArray3 _rC = reinterpret_cast < fullArray3 > (*_local_C);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[9][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][0][chunk_offset]);
        refArray2 Bk =
          reinterpret_cast < refArray2 > (_rB[index][0][chunk_offset]);
        refArray3 Ck =
          reinterpret_cast < refArray3 > (_rC[index][0][chunk_offset]);

        Matrix_Times_Transpose < SIMDArchitectureAVX512<float>, float[16] > (Ak, Bk, Ck);
      }

  }
};

#endif


//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================
#ifdef ENABLE_AVX_INSTRUCTION_SET
template<int data_size>
void ComputeAVX( float A[data_size][9][16],  float B[data_size][9][16],  float C[data_size][9][16] )
    {
      std::cout << "	Running " << data_size << " of AVX :  " << std::endl;

      typedef float (*Mat3_ptr_type)  [data_size][9][16];
      Mat3_ptr_type A_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&A[0][0][0]);
      Mat3_ptr_type B_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&B[0][0][0]);
      Mat3_ptr_type C_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&C[0][0][0]);

      Mat3_ptr_type A2_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&A[0][0][8]);
      Mat3_ptr_type B2_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&B[0][0][8]);
      Mat3_ptr_type C2_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&C[0][0][8]);




      //Matrix_Times_Transpose_AVX_None < data_size > op ((float *) &A,
      //                                                  (float *) &B,
      //                                                  (float *) &C);

      //for (int t = threads; t <= threads_max;
      //     t += std::max < int >(((threads_max - threads) / 30), 1))
      //  {
      //    std::cout << "Running Test with " << t << " threads." << std::endl;
          //MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
          //  Matrix_Times_Transpose_AVX_None < data_size > >helper (op,
          //                                                        data_size,
          //                                                         t);

      //    double min_time = 10000000;
      //    double max_time = -1;
      //    double avg_time = 0;

      //    for (int i = 0; i < passes; i++)
      //      {
              start_timer ();
#if 0
              helper.Run_Parallel ();
#else

              for(int j=0;j<data_size;j++){
                  Matrix_Times_Transpose < SIMDArchitectureAVX2<float>, float[16]> (
                         reinterpret_cast< const float (&)[9][16] >((*A_ptr)[j]),
                         reinterpret_cast< const float (&)[9][16] >((*B_ptr)[j]),
                         reinterpret_cast< float (&)[9][16] >((*C_ptr)[j]) );
                  Matrix_Times_Transpose < SIMDArchitectureAVX2<float>, float[16] > (
                         reinterpret_cast< const float (&)[9][16] >((*A2_ptr)[0]),
                         reinterpret_cast< const float (&)[9][16] >((*B2_ptr)[0]),
                         reinterpret_cast< float (&)[9][16] >((*C2_ptr)[0]) );
              }
#endif
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              //       min_time = std::min < double >(min_time, get_time ());
              // max_time = std::max < double >(max_time, get_time ());
              //avg_time += get_time ();
              // }
              //avg_time = avg_time / passes;
              //std::cout << "Min pass time: " << min_time << std::endl;
              //std::cout << "Max pass time: " << max_time << std::endl;
              //std::cout << "Avg pass time: " << avg_time << std::endl;
              // }




    }
#endif




int
main (int argc, char *argv[])
{
  typedef float T;

  int seed = 1;
  int threads = 1;
  int threads_max = 1;
  int passes = 1;
  const int data_size = 1000000;
  if (argc >= 2){
    threads = atoi (argv[1]);
    threads_max = threads;
  }
  if (argc >= 3)
    threads_max = atoi (argv[2]);
  if (argc >= 4)
    passes = atoi (argv[3]);
  srand (seed);

  pthread_queue = new PTHREAD_QUEUE (threads_max);

  std::
    cout << "Preparing to Run " << data_size << " of all kernels with " <<
    threads << " threads." << std::endl;



  {
    std::cout << "Running Thread Test for Matrix_Times_Transpose " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&A_type)[data_size][9][16];
    A_type A =
      reinterpret_cast < A_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&B_type)[data_size][9][16];
    B_type B =
      reinterpret_cast < B_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&C_type)[data_size][9][16];
    C_type C =
      reinterpret_cast < C_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            A[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            B[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            C[__a][__b][__c] = Get_Random < float >();
          }
    stop_timer ();

    std::cout << get_time () << "s\n\n" << std::endl;


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================
#if 0
    {
      std::cout << "	Running " << data_size << " of SCALAR :  " << std::endl;


      Matrix_Times_Transpose_SCALAR_None < data_size > op ((float *) &A,
                                                           (float *) &B,
                                                           (float *) &C);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Matrix_Times_Transpose_SCALAR_None < data_size > >helper (op,
                                                                      data_size,
                                                                      t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================
#if 0
#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of AVX :  " << std::endl;

      typedef float (*Mat3_ptr_type)  [data_size][9][16];
      Mat3_ptr_type A_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&A[0][0][0]);
      Mat3_ptr_type B_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&B[0][0][0]);
      Mat3_ptr_type C_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&C[0][0][0]);

      Mat3_ptr_type A2_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&A[0][0][8]);
      Mat3_ptr_type B2_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&B[0][0][8]);
      Mat3_ptr_type C2_ptr=reinterpret_cast< float (*) [data_size][9][16] >(&C[0][0][8]);




      //Matrix_Times_Transpose_AVX_None < data_size > op ((float *) &A,
      //                                                  (float *) &B,
      //                                                  (float *) &C);

      //for (int t = threads; t <= threads_max;
      //     t += std::max < int >(((threads_max - threads) / 30), 1))
      //  {
      //    std::cout << "Running Test with " << t << " threads." << std::endl;
          //MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
          //  Matrix_Times_Transpose_AVX_None < data_size > >helper (op,
          //                                                        data_size,
          //                                                         t);

      //    double min_time = 10000000;
      //    double max_time = -1;
      //    double avg_time = 0;

      //    for (int i = 0; i < passes; i++)
      //      {
              start_timer ();
#if 0
              helper.Run_Parallel ();
#else

              for(int j=0;j<data_size;j++){
                  Matrix_Times_Transpose < __m256, float[16], int[16] > (
                         reinterpret_cast< const float (&)[9][16] >((*A_ptr)[j]),
                         reinterpret_cast< const float (&)[9][16] >((*B_ptr)[j]),
                         reinterpret_cast< float (&)[9][16] >((*C_ptr)[j]) );
                  Matrix_Times_Transpose < __m256, float[16], int[16] > (
                         reinterpret_cast< const float (&)[9][16] >((*A2_ptr)[j]),
                         reinterpret_cast< const float (&)[9][16] >((*B2_ptr)[j]),
                         reinterpret_cast< float (&)[9][16] >((*C2_ptr)[j]) );
              }
#endif
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              //       min_time = std::min < double >(min_time, get_time ());
              // max_time = std::max < double >(max_time, get_time ());
              //avg_time += get_time ();
              // }
              //avg_time = avg_time / passes;
              //std::cout << "Min pass time: " << min_time << std::endl;
              //std::cout << "Max pass time: " << max_time << std::endl;
              //std::cout << "Avg pass time: " << avg_time << std::endl;
              // }




    }
#endif
#else
    ComputeAVX<data_size>( A, B, C );
#endif


//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================
#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of MIC :  " << std::endl;


      Matrix_Times_Transpose_MIC_None < data_size > op ((float *) &A,
                                                        (float *) &B,
                                                        (float *) &C);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Matrix_Times_Transpose_MIC_None < data_size > >helper (op,
                                                                   data_size,
                                                                   t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//        FREE MEMORY USED BY ALL VARIABLES
//
//=======================================================
    std::cout << "\nFreeing all data: " << std::endl;
    std::cout.flush ();

    _mm_free (reinterpret_cast < void *>(A));
    _mm_free (reinterpret_cast < void *>(B));
    _mm_free (reinterpret_cast < void *>(C));


  }


  return 0;

}
