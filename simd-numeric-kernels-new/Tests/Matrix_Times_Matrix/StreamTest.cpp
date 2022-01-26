#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;
#include <omp.h>
#include "Matrix_Times_Matrix.h"

#define NUM_TRIALS 1000000

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

float A[9][16] __attribute__ ((aligned (64)));
float B[9][16] __attribute__ ((aligned (64)));
float C[9][16] __attribute__ ((aligned (64)));
float C_reference[9][16] __attribute__ ((aligned (64)));
float C_original[9][16] __attribute__ ((aligned (64)));

#pragma omp threadprivate(A, B, C, C_reference, C_original)

int
main (int argc, char *argv[])
{
    using namespace SIMD_Numeric_Kernel;
  typedef float T;

  std::cout << "Preparing to Run " << NUM_TRIALS << " of all kernels." << std::
    endl;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);

  omp_set_dynamic(0);

  {
    std::cout << "Running Stream Test for Matrix_Times_Matrix " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        A[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        B[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          C_original[__a][__b] = Get_Random < float >();
          C[__a][__b] = C_original[__a][__b];
          C_reference[__a][__b] = C_original[__a][__b];
        }


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
      start_timer ();
#pragma omp parallel for copyin(A,B)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[9][16];
          typedef T (&refArray2)[9][16];
          typedef T (&refArray3)[9][16];
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Bk = reinterpret_cast < refArray2 > (B[0][i]);
              refArray3 Ck = reinterpret_cast < refArray3 > (C[0][i]);
              Matrix_Times_Matrix < SIMDArchitectureScalar<float>, float[16]> (Ak, Bk, Ck);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of AVX :  ";
      start_timer ();
#pragma omp parallel for copyin(A,B)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[9][16];
          typedef T (&refArray2)[9][16];
          typedef T (&refArray3)[9][16];
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Bk = reinterpret_cast < refArray2 > (B[0][i]);
              refArray3 Ck = reinterpret_cast < refArray3 > (C[0][i]);
              Matrix_Times_Matrix < SIMDArchitectureAVX2<float>, float[16] > (Ak, Bk,
                                                                     Ck);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of MIC :  ";
      start_timer ();
#pragma omp parallel for copyin(A,B)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[9][16];
          typedef T (&refArray2)[9][16];
          typedef T (&refArray3)[9][16];
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Bk = reinterpret_cast < refArray2 > (B[0][i]);
              refArray3 Ck = reinterpret_cast < refArray3 > (C[0][i]);
              Matrix_Times_Matrix < SIMDArchitectureAVX512<float>, float[16] > (Ak, Bk,
                                                                     Ck);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  return 0;

}
