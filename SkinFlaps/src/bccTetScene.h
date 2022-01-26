//////////////////////////////////////////////////////////////////
// File: bccTetScene.h
// Author: Court Cutting
// Date: 3/4/2019
// Purpose: Bcc tet based projective dynamics physics library interface to surgical simulator code.
//    This version uses the ShapeOp subroutine library of Sofien Bouaziz ( http://ShapeOp.org )
//    to do the physics.  So far have experimented with mass-springs, physX(NVIDIA),
//    physBAM(Fedkiw, Teran, Sifakis, et al.), corotated linear elasticity(Teran, Sifakis, Mitchell) and
//    the shape matching code of Rivers and James.
//    Copyright 2019 - All rights reserved at the present time.
///////////////////////////////////////////////////////////////////

#ifndef __BCC_TET_SCENE__
#define __BCC_TET_SCENE__

#include "surgGraphics.h"
#include "vnBccTetrahedra.h"
//#include "vnBccTetCutter.h"
#include "vnBccTetCutterTbb.h"
#include "tetCollisions.h"
#include "tetSubset.h"
#include "pdTetPhysics.h"

// forward declarations
class gl3wGraphics;
class surgicalActions;

class bccTetScene
{
public:
	bool loadScene(const char *dataDirectory, const char *sceneFileName);
	void createNewPhysicsLattice(int maximumDimensionSubdivisions);
	void updateOldPhysicsLattice();
	inline void nonTetPhysicsUpdate() {_ptp.initializePhysics();}
	void updatePhysics();
	void setFixedRegion(const float(&corners)[6]);
	void fixPeriostealPeriferalVertices();
	void updateSurfaceDraw();
	pdTetPhysics* getPdTetPhysics_2(){ return &_ptp; }
	inline void setForcesAppliedFlag(){ _forcesApplied = true; }
	inline void promoteSutures() { _ptp.promoteAllSutures(); _ptp.initializePhysics(); }
	vnBccTetrahedra* getVirtualNodedBccTetrahedra() { return &_vnTets; }
	void setVisability(char surface, char physics);	// 0=off, 1=on, 2=don't change
	void setGl3wGraphics(gl3wGraphics *gl3w) { _gl3w = gl3w; }
	void createTetLatticeDrawing();
	void drawTetLattice();
	void eraseTetLattice();
	void setSurgicalActions(surgicalActions *sa) { _surgAct = sa; }
	void setPhysicsPause(bool pause) { _physicsPaused = pause; }
	inline bool isPhysicsPaused(){ return  _physicsPaused; }
	bccTetScene();
	~bccTetScene();

private:
	gl3wGraphics *_gl3w;
	surgicalActions *_surgAct;
	vnBccTetrahedra _vnTets;
	tetCollisions _tetCol;
	tetSubset _tetSubsets;
//	vnBccTetCutter _tc;  // original serial version giving deterministic outcome, but slow.
	vnBccTetCutterTbb _tc;  // multithreaded version using Intel threaded building blocks.  Much faster, but indices of nodes and tets different each run as nondeterministic.
	pdTetPhysics _ptp;
	bool _forcesApplied;
	bool _tetsModified;

	materialTriangles *_mt;  // pointer from surgGraphics.
	struct boundingBox3{
		float corners[6];
	};
	std::vector<GLfloat> _nodeGraphicsPositions;  // homogeneous coords[4]

	std::vector<Vec3f> _firstSpatialCoords;

	bool _physicsPaused;
	void initPdPhysics();
};

#endif // __BCC_TET_SCENE__