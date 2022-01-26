//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#pragma once

// Or any other 16-wide arch
#if defined(ENABLE_MIC_INSTRUCTION_SET)
#define BUILD_CONSTANT(name,value) const float name[16] __attribute__((aligned(64))) = {value,value,value,value, \
                                                                                      value,value,value,value, \
                                                                                      value,value,value,value, \
                                                                                      value,value,value,value};
#define BUILD_ICONSTANT(name,value) const int name[16] __attribute__((aligned(64))) = {value,value,value,value, \
                                                                                      value,value,value,value, \
                                                                                      value,value,value,value, \
                                                                                      value,value,value,value};
#define BUILD_TDATA(name,size) T_DATA name size __attribute__((aligned(64)));
#define BUILD_IDATA(name,size) I_DATA name size __attribute__((aligned(64)));
#else

// Or any other 8-wide arch
#if defined(ENABLE_AVX_INSTRUCTION_SET)
#define BUILD_CONSTANT(name,value) const float name[8] __attribute__((aligned(32))) = {value,value,value,value, \
                                                                                     value,value,value,value};
#define BUILD_ICONSTANT(name,value) const int name[8] __attribute__((aligned(32))) = {value,value,value,value, \
                                                                                     value,value,value,value};
#define BUILD_TDATA(name,size) T_DATA name size __attribute__((aligned(32)));
#define BUILD_IDATA(name,size) I_DATA name size __attribute__((aligned(32)));
#else
// Fallback if no vector arch's are enabled.
#define BUILD_CONSTANT(name,value) const float name[1] __attribute__((aligned(4))) = {value};
#define BUILD_ICONSTANT(name,value) const int name[1] __attribute__((aligned(4))) = {value};
#define BUILD_TDATA(name,size) T_DATA name size __attribute__((aligned(4)));
#define BUILD_IDATA(name,size) I_DATA name size __attribute__((aligned(4)));

#endif // End Size 8 group
#endif // End Size 16 group

#define BUILD_CONSTANT_16(name,value) const float name[16] __attribute__((aligned(64))) = { \
    value,value,value,value, \
    value,value,value,value, \
    value,value,value,value, \
    value,value,value,value};

#define BUILD_CONSTANT_8(name,value) const float name[8] __attribute__((aligned(32))) = { \
    value,value,value,value, \
    value,value,value,value};

#define BUILD_CONSTANT_4(name,value) const float name[4] __attribute__((aligned(16))) = { \
    value,value,value,value};

#define BUILD_CONSTANT_1(name,value) const float name[1] __attribute__((aligned(4))) = { \
    value};

#define BUILD_ICONSTANT_16(name,value) const int name[16] __attribute__((aligned(64))) = { \
    value,value,value,value, \
    value,value,value,value, \
    value,value,value,value, \
    value,value,value,value};

#define BUILD_ICONSTANT_8(name,value) const int name[8] __attribute__((aligned(32))) = { \
    value,value,value,value, \
    value,value,value,value};

#define BUILD_ICONSTANT_4(name,value) const int name[4] __attribute__((aligned(16))) = { \
    value,value,value,value};

#define BUILD_ICONSTANT_1(name,value) const int name[1] __attribute__((aligned(4))) = { \
    value};
