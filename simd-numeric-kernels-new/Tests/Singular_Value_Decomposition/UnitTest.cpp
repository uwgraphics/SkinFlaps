#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Singular_Value_Decomposition.h"
#include "Singular_Value_Decomposition_Reference.h"

template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}

int
main (int argc, char *argv[])
{
    using namespace SIMD_Numeric_Kernel;
  typedef float T;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);



  {
    T A[9] __attribute__ ((aligned (4)));
    T U[9] __attribute__ ((aligned (4)));
    T U_reference[9] __attribute__ ((aligned (4)));
    T U_original[9] __attribute__ ((aligned (4)));
    T S[3] __attribute__ ((aligned (4)));
    T S_reference[3] __attribute__ ((aligned (4)));
    T S_original[3] __attribute__ ((aligned (4)));
    T V[9] __attribute__ ((aligned (4)));
    T V_reference[9] __attribute__ ((aligned (4)));
    T V_original[9] __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 9; __a++)
      A[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      {
        U_original[__a] = Get_Random < float >();
        U[__a] = U_original[__a];
        U_reference[__a] = U_original[__a];
      }
    for (int __a = 0; __a < 3; __a++)
      {
        S_original[__a] = Get_Random < float >();
        S[__a] = S_original[__a];
        S_reference[__a] = S_original[__a];
      }
    for (int __a = 0; __a < 9; __a++)
      {
        V_original[__a] = Get_Random < float >();
        V[__a] = V_original[__a];
        V_reference[__a] = V_original[__a];
      }


    for (int __a = 0; __a < 9; __a++)
      U[__a] = U_original[__a];
    for (int __a = 0; __a < 3; __a++)
      S[__a] = S_original[__a];
    for (int __a = 0; __a < 9; __a++)
      V[__a] = V_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
          Singular_Value_Decomposition < SIMDArchitectureScalar<float>, float>(A, U, S, V);
      }

    Singular_Value_Decomposition_Reference < float >(A, U_reference,
                                                     S_reference, V_reference);
    if (!
        (Singular_Value_Decomposition_Compare <
         float >(U, S, V, U_reference, S_reference, V_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Singular_Value_Decomposition " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}
