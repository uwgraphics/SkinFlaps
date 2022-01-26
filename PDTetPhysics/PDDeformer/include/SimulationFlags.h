//#####################################################################
// Copyright (c) 2019, Eftychios Sifakis, Yutian Tao, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################

#pragma once

namespace PDSimulation {

    enum NodeType { InactiveNode = 0, ActiveNode = 1, DirichletNode = 2, CollisionNode = 3 };
    enum ElementFlag {unCollisionEl = 1, CollisionEl = 2};

}
