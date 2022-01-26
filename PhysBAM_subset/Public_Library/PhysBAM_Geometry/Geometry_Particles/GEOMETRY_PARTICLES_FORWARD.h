//#####################################################################
// Copyright 2009-2013, Craig Schroeder, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GEOMETRY_PARTICLES_FORWARD__
#define __GEOMETRY_PARTICLES_FORWARD__

#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
namespace PhysBAM{
template<class TV> class GEOMETRY_PARTICLES;
const ATTRIBUTE_ID ATTRIBUTE_ID_V(2);
const ATTRIBUTE_ID ATTRIBUTE_ID_ROTATION(3);
const ATTRIBUTE_ID ATTRIBUTE_ID_ANGULAR_VELOCITY(4);
const ATTRIBUTE_ID ATTRIBUTE_ID_RIGID_GEOMETRY(5);
const ATTRIBUTE_ID ATTRIBUTE_ID_STRUCTURE_IDS(6);
const ATTRIBUTE_ID ATTRIBUTE_ID_COLLIDABLE(30);
const ATTRIBUTE_ID ATTRIBUTE_ID_COLOR(31);
const ATTRIBUTE_ID ATTRIBUTE_ID_DISPLAY_SIZE(32);
const ATTRIBUTE_ID ATTRIBUTE_ID_RADIUS(15);
const ATTRIBUTE_ID ATTRIBUTE_ID_BRADIUS(50);
const ATTRIBUTE_ID ATTRIBUTE_ID_BMASS(51);
const ATTRIBUTE_ID ATTRIBUTE_ID_PRESSURE(52);
const ATTRIBUTE_ID ATTRIBUTE_ID_RADIAL_VELOCITY(53);
const ATTRIBUTE_ID ATTRIBUTE_ID_WEIGHTS(54);
const ATTRIBUTE_ID ATTRIBUTE_ID_CUMULATIVE(55);
const ATTRIBUTE_ID ATTRIBUTE_ID_RELATIVE_VELOCITY(56);
const ATTRIBUTE_ID ATTRIBUTE_ID_SIGNED_DISTANCE(57);
const ATTRIBUTE_ID ATTRIBUTE_ID_ANGULAR(34);
const ATTRIBUTE_ID ATTRIBUTE_ID_THICKNESS(60);
const ATTRIBUTE_ID ATTRIBUTE_ID_THIN_VOLUME_THICKNESS(74);
const ATTRIBUTE_ID ATTRIBUTE_ID_REDUCED_BASIS(61);
const ATTRIBUTE_ID ATTRIBUTE_ID_REDUCED_DISPLACEMENTS(62);
const ATTRIBUTE_ID ATTRIBUTE_ID_REDUCED_VELOCITIES(63);
const ATTRIBUTE_ID ATTRIBUTE_ID_REDUCED_REST_GEOMETRY(64);
const ATTRIBUTE_ID ATTRIBUTE_ID_ANGULAR_ACCELERATION(65);
const ATTRIBUTE_ID ATTRIBUTE_ID_INCOMPRESSIBILITY(66);
const ATTRIBUTE_ID ATTRIBUTE_ID_BUBBLE(67);
const ATTRIBUTE_ID ATTRIBUTE_ID_SURFACE_TENSION_FORCE(68);
const ATTRIBUTE_ID ATTRIBUTE_ID_PRESSURE_FORCE(69);
const ATTRIBUTE_ID ATTRIBUTE_ID_BODY_FORCE(70);
const ATTRIBUTE_ID ATTRIBUTE_ID_AREA_WEIGHTED_NORMAL(71);
const ATTRIBUTE_ID ATTRIBUTE_ID_OLD_ID(72);
const ATTRIBUTE_ID ATTRIBUTE_ID_DIMENSION(73);
const ATTRIBUTE_ID ATTRIBUTE_ID_REDUCED_SCALING_FACTOR(75);
const ATTRIBUTE_ID ATTRIBUTE_ID_REDUCED_LAZY_DISPLACEMENTS(81);
const ATTRIBUTE_ID ATTRIBUTE_ID_REDUCED_LAZY_VELOCITIES(82);
const ATTRIBUTE_ID ATTRIBUTE_ID_DIMENSIONAL_MASS(83);
const ATTRIBUTE_ID ATTRIBUTE_ID_GUID(84);
const ATTRIBUTE_ID ATTRIBUTE_ID_VISCOSITY_FORCE(85);
const ATTRIBUTE_ID ATTRIBUTE_ID_DISTANCE_TO_RIM(86);
const ATTRIBUTE_ID ATTRIBUTE_ID_TID(91);
const ATTRIBUTE_ID ATTRIBUTE_ID_LOCAL_INDEX(92);
const ATTRIBUTE_ID ATTRIBUTE_ID_DIMENSIONAL_SIZE(93);
const ATTRIBUTE_ID ATTRIBUTE_ID_THREAD_BOUNDARY(94);
const ATTRIBUTE_ID ATTRIBUTE_ID_RGB_COLOR(95);
const ATTRIBUTE_ID ATTRIBUTE_ID_THREAD_BOUNDARY_VECTOR(96);
const ATTRIBUTE_ID ATTRIBUTE_ID_DISTANCE_TO_VOLUME(97);
const ATTRIBUTE_ID ATTRIBUTE_ID_SMALL_REGION_ID(98);
}
#endif
