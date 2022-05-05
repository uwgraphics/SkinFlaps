//#####################################################################
// Copyright (c) 2019, Eftychios Sifakis, Yutian Tao, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################

#pragma once

    enum class NodeType { Inactive = 0, Active = 1, Dirichlet = 2, Collision = 3 };
    enum class ElementFlag {unCollisionEl = 1, CollisionEl = 2};
