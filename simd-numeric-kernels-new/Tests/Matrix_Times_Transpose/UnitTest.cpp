#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Matrix_Times_Transpose.h"
#include "Matrix_Times_Transpose_Reference.h"

template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}

int
main (int argc, char *argv[])
{
    using namespace SIMD_Numeric_Kernel;
    using T = float;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);



  {
    T A[9] __attribute__ ((aligned (4)));
    T B[9] __attribute__ ((aligned (4)));
    T C[9] __attribute__ ((aligned (4)));
    T C_reference[9] __attribute__ ((aligned (4)));
    T C_original[9] __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 9; __a++)
      A[__a] = Get_Random < T >();
    for (int __a = 0; __a < 9; __a++)
      B[__a] = Get_Random < T >();
    for (int __a = 0; __a < 9; __a++)
      {
        C_original[__a] = Get_Random < T >();
        C[__a] = C_original[__a];
        C_reference[__a] = C_original[__a];
      }


    for (int __a = 0; __a < 9; __a++)
      C[__a] = C_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Matrix_Times_Transpose < SIMDArchitectureScalar<T>, T>(A, B, C);
      }

    Matrix_Times_Transpose_Reference < T >(A, B, C_reference);
    if (!(Matrix_Times_Transpose_Compare < T >(C, C_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Matrix_Times_Transpose " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}
