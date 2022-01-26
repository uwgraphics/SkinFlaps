 #pragma once

template<class Tarch,class T_DATA>
void Add_Force(const T_DATA (&x_Blocked)[4][3],
               const T_DATA (&DmInverse_Blocked)[9],
               const T_DATA &restVolume,
               const T_DATA &muLow,
               const T_DATA &muHigh,
               const T_DATA &strainMin,
               const T_DATA &strainMax,
               T_DATA (&f_Blocked)[4][3]);
