#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Singular_Value_Decomposition.h"

#define NUM_TRIALS 1000000

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

int
main (int argc, char *argv[])
{
  typedef float T;

  std::cout << "Preparing to Run " << NUM_TRIALS << " of all kernels." << std::
    endl;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);



  {
    std::
      cout << "Running Stream Test for Singular_Value_Decomposition " << std::
      endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T A[9][16] __attribute__ ((aligned (64)));
    T U[9][16] __attribute__ ((aligned (64)));
    T U_reference[9][16] __attribute__ ((aligned (64)));
    T U_original[9][16] __attribute__ ((aligned (64)));
    T S[3][16] __attribute__ ((aligned (64)));
    T S_reference[3][16] __attribute__ ((aligned (64)));
    T S_original[3][16] __attribute__ ((aligned (64)));
    T V[9][16] __attribute__ ((aligned (64)));
    T V_reference[9][16] __attribute__ ((aligned (64)));
    T V_original[9][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        A[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          U_original[__a][__b] = Get_Random < float >();
          U[__a][__b] = U_original[__a][__b];
          U_reference[__a][__b] = U_original[__a][__b];
        }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          S_original[__a][__b] = Get_Random < float >();
          S[__a][__b] = S_original[__a][__b];
          S_reference[__a][__b] = S_original[__a][__b];
        }
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          V_original[__a][__b] = Get_Random < float >();
          V[__a][__b] = V_original[__a][__b];
          V_reference[__a][__b] = V_original[__a][__b];
        }


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[9][16];
          typedef T (&refArray2)[9][16];
          typedef T (&refArray3)[3][16];
          typedef T (&refArray4)[9][16];
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Uk = reinterpret_cast < refArray2 > (U[0][i]);
              refArray3 Sk = reinterpret_cast < refArray3 > (S[0][i]);
              refArray4 Vk = reinterpret_cast < refArray4 > (V[0][i]);
              Singular_Value_Decomposition < SIMDArchitectureScalar<float>, float[16]> (Ak,
                                                                          Uk,
                                                                          Sk,
                                                                          Vk);
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
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[9][16];
          typedef T (&refArray2)[9][16];
          typedef T (&refArray3)[3][16];
          typedef T (&refArray4)[9][16];
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Uk = reinterpret_cast < refArray2 > (U[0][i]);
              refArray3 Sk = reinterpret_cast < refArray3 > (S[0][i]);
              refArray4 Vk = reinterpret_cast < refArray4 > (V[0][i]);
              Singular_Value_Decomposition < SIMDArchitectureAVX2<float>, float[16]> (Ak,
                                                                           Uk,
                                                                           Sk,
                                                                           Vk);
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
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[9][16];
          typedef T (&refArray2)[9][16];
          typedef T (&refArray3)[3][16];
          typedef T (&refArray4)[9][16];
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Uk = reinterpret_cast < refArray2 > (U[0][i]);
              refArray3 Sk = reinterpret_cast < refArray3 > (S[0][i]);
              refArray4 Vk = reinterpret_cast < refArray4 > (V[0][i]);
              Singular_Value_Decomposition < SIMDArchitectureAVX512<float>, float[16] > (Ak,
                                                                           Uk,
                                                                           Sk,
                                                                           Vk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  return 0;

}
