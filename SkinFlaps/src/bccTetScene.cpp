//////////////////////////////////////////////////////////////////
// File: bccTetScene.cpp
// Author: Court Cutting
// Date: 3/4/2019
// Purpose: Bcc tet based projective dynamics physics library interface to surgical simulator code.
//    This version uses the ShapeOp subroutine library of Sofien Bouaziz ( http://ShapeOp.org )
//    to do the physics.  So far have experimented with mass-springs, physX(NVIDIA),
//    physBAM(Fedkiw, Teran, Sifakis, et al.), corotated linear elasticity(Teran, Sifakis, Mitchell) and
//    the shape matching code of Rivers and James.
//    Copyright 2019 - All rights reserved at the present time.
///////////////////////////////////////////////////////////////////

#include "bccTetScene.h"
#include <string>
#include <fstream>
#include <algorithm>
#include "gl3wGraphics.h"
#include "surgicalActions.h"
#include "boundingBox.h"
#include "json.h"
#include "closestPointOnTriangle.h"
#include "remapTetPhysics.h"
#include <iostream>
#include <chrono>
#include <ctime>

bool bccTetScene::loadScene(const char *dataDirectory, const char *sceneFileName)
{
	_physicsPaused = true;
	std::string path(dataDirectory);
	path.append(sceneFileName);
	std::ifstream istr(path.c_str());
	std::string jsonStr;
	if (!istr.is_open()) {
		path = std::string("Unable to load: ") + path;
		_surgAct->sendUserMessage(path.c_str(), "Error Message");
		istr.close();
		return false;
	}
	else {
		char ch;
		while (istr.get(ch))
			jsonStr.push_back(ch);

	}
	istr.close();
	json::Value my_data = json::Deserialize(jsonStr);  // will trim leading and trailing white space from {} pair
	if (my_data.GetType() != json::ObjectVal) {
		_surgAct->sendUserMessage("Module file not in correct JSON format-", "Error Message");
		return false;
	}
	//	_fixedMeshes.clear();
	//	_fixedRegions.clear();
	json::Object scnObj = my_data.ToObject();
	json::Object::ValueMap::iterator oit, suboit, suboit2;
	// get texture files first
	std::map<int, GLuint> txMap;
	std::string nrm, tex;
	if ((oit = scnObj.find("textureFiles")) == scnObj.end()) {
		_surgAct->sendUserMessage("No texture files in scene file-", "Error Message");
		return false;
	}
	else {
		json::Object txObj = oit->second.ToObject();
		for (suboit = txObj.begin(); suboit != txObj.end(); ++suboit) {
			path = dataDirectory + suboit->first;
			GLuint txNow = _gl3w->getTextures()->loadTexture(suboit->second.ToInt(), path.c_str());
			if (txNow > 0xfffffffe) {
				path = "Unable to load bitmap .bmp input file: " + path;
				_surgAct->sendUserMessage(path.c_str(), "Error Message");
				return false;
			}
			int txNum = suboit->second.ToInt();
			txMap.insert(std::make_pair(txNum, txNow));
		}
	}
	if ((oit = scnObj.find("staticObjects")) != scnObj.end()) {
		json::Object statObj = oit->second.ToObject();
		for (suboit = statObj.begin(); suboit != statObj.end(); ++suboit) {
			path = dataDirectory + suboit->first;
			std::map<int, GLuint>::iterator tit;
			std::vector<int> txIds;
			json::Object tmapObj = suboit->second.ToObject();
			for (suboit2 = tmapObj.begin(); suboit2 != tmapObj.end(); ++suboit2) {
				if (suboit2->first == "textureMap")
					txIds.push_back(suboit2->second.ToInt());
				else if (suboit2->first == "normalMap")
					txIds.push_back(suboit2->second.ToInt());
				else {
					_surgAct->sendUserMessage("Incorrect static object section in .smd input file-", "Error Message");
					return false;
				}
			}
			// this is a staticTriangle, not elastic so put on graphics card and clean up
			if ( _gl3w->loadStaticObjFile(path.c_str(), txIds, true) == NULL)
			{
				_surgAct->sendUserMessage("Unable to load fixed triangle .obj input file-", "Error Message");
				return false;
			}
		}
	}
	std::string deepBedFilepath;
	deepBedFilepath.clear();
	if ((oit = scnObj.find("dynamicObjects")) == scnObj.end()) {
		_surgAct->sendUserMessage("No dynamic objects in this scene file-", "Error Message");
		return false;
	}
	else {
		json::Object dynObj = oit->second.ToObject();
		std::vector<int> txIds;
		for (suboit = dynObj.begin(); suboit != dynObj.end(); ++suboit) {
			path = dataDirectory + suboit->first;
			std::map<int, GLuint>::iterator tit;
			json::Object tmapObj = suboit->second.ToObject();
			for (suboit2 = tmapObj.begin(); suboit2 != tmapObj.end(); ++suboit2) {
				if (suboit2->first == "textureMaps") {
					json::Array txArr;
					txArr = suboit2->second.ToArray();
					for (int i = 0; i < txArr.size(); ++i) {
						txIds.push_back(txArr[i].ToInt());
						if (!_gl3w->getTextures()->textureExists(txIds.back())) {
							_surgAct->sendUserMessage("Missing texture or normal map in dynamic triangle section in .smd input file-", "Error Message");
							return false;
						}
					}
				}
			}
			_mt = _surgAct->getSurgGraphics()->getMaterialTriangles();
			if (_mt->readObjFile(path.c_str())) {
				_surgAct->sendUserMessage("Unable to load fixed materialTriangle .obj input file-", "Error Message");
				return false;
			}
			_mt->collectCreateTextureSeams();  // must do this on load to register same material texture seams and create hard texture & normal seams between materials.
			_surgAct->getSurgGraphics()->setGl3wGraphics(_gl3w);
			std::string vtxShd(dataDirectory), frgShd(dataDirectory);
			vtxShd.append("mtVertexShader.txt");
			frgShd.append("mtFragmentShader.txt");
			_surgAct->getSurgGraphics()->setTextureFilesCreateProgram(txIds, vtxShd.c_str(), frgShd.c_str());  // openGL buffers ceated here
			_surgAct->getSurgGraphics()->setNewTopology();
			_surgAct->getSurgGraphics()->updatePositionsNormalsTangents();
			_surgAct->getSurgGraphics()->computeLocalBounds();
			path = suboit->first;
			size_t pos = path.rfind(".obj");
			path.erase(pos);
			_mt->setName(path.c_str());
			_surgAct->getSurgGraphics()->getSceneNode()->setName(path.c_str());
			// input new deep bed file here
			deepBedFilepath.clear();
			deepBedFilepath.append(dataDirectory);
			deepBedFilepath.append(path);
			deepBedFilepath.append(".bed");
		}
	}
	if ((oit = scnObj.find("fixedGeometry")) != scnObj.end()) {
		json::Object fixObj = oit->second.ToObject();
		for (suboit = fixObj.begin(); suboit != fixObj.end(); ++suboit) {
			if (suboit->first == "mesh") {
				json::Array meshArr = suboit->second.ToArray();
				json::Array::ValueVector::iterator ait;
				for (ait = meshArr.begin(); ait != meshArr.end(); ++ait) {
					path = dataDirectory + ait->ToString();
					//					_fixedMeshes.push_back(staticTriangle());
					//					staticTriangle *tuv = &_fixedMeshes.back();
					//					if (tuv->readObjFile(path.c_str())) {
					//						_fixedMeshes.pop_back();
					//						_surgAct->sendUserMessage("Unable to load fixed geometry mesh file-", "Error Message");
					//						return false;
					//					}
					//					tuv->setName(ait->ToString().c_str());
				}
			}
			//			if (suboit->first == "region") {  // currently unused
			//				json::Array bbArr;
			//				bbArr = suboit->second.ToArray();
			//				_fixedRegions.push_back(boundingBox3());
			//				for (int i = 0; i < 6; ++i)
			//					_fixedRegions.back().corners[i] = bbArr[i];
			//			}
		}
	}
	if ((oit = scnObj.find("fixedCollisionSets")) != scnObj.end()) {
		json::Object hullObj = oit->second.ToObject();
		std::string lsPath;
		std::vector<int> vIdx;
		for (suboit = hullObj.begin(); suboit != hullObj.end(); ++suboit) {
			lsPath = dataDirectory + suboit->first;
			json::Array polyArr;
			polyArr = suboit->second.ToArray();
			vIdx.reserve(polyArr.size());
			for (int i = 0; i < polyArr.size(); ++i)
				vIdx.push_back( polyArr[i].ToInt());
			_tetCol.addFixedCollisionSet(lsPath, vIdx);
		}
	}
	int maxDimSubdivs = 180;  // 0.5 million tets for cleft model
	if ((oit = scnObj.find("tetrahedralProperties")) != scnObj.end()) {
		json::Object hullObj = oit->second.ToObject();
		float lowTetWeight, highTetWeight, strainMin, strainMax, collisionWeight, fixedWeight, periferalWeight, hookWeight, sutureWeight, autoSutureSpacing, selfCollisionWeight;
		std::string lsFname(dataDirectory), collisionObject;
		for (suboit = hullObj.begin(); suboit != hullObj.end(); ++suboit) {
			if (suboit->first == "minStrain")
				strainMin = suboit->second.ToFloat();
			else if (suboit->first == "maxStrain")
				strainMax = suboit->second.ToFloat();
			else if (suboit->first == "lowTetWeight")
				lowTetWeight = suboit->second.ToFloat();
			else if (suboit->first == "highTetWeight")
				highTetWeight = suboit->second.ToFloat();
			else if (suboit->first == "collisionWeight")
				collisionWeight = suboit->second.ToFloat();
			else if (suboit->first == "selfCollisionWeight")
				selfCollisionWeight = suboit->second.ToFloat();
			else if (suboit->first == "fixedWeight")
				fixedWeight = suboit->second.ToFloat();
			else if (suboit->first == "periferalWeight")
				periferalWeight = suboit->second.ToFloat();
			else if (suboit->first == "sutureWeight")
				sutureWeight = suboit->second.ToFloat();
			else if (suboit->first == "hookWeight")
				hookWeight = suboit->second.ToFloat();
			else if (suboit->first == "autoSutureSpacing")
				autoSutureSpacing = suboit->second.ToFloat();
			else if (suboit->first == "collisionObject")
				collisionObject = suboit->second.ToString();
			else if (suboit->first == "maxDimSubdivisions")
				maxDimSubdivs = suboit->second.ToInt();
			else
				_surgAct->sendUserMessage("Unknown tetrahedral property in scene file-", "File Error Message");
		}
		_ptp.setTetProperties(lowTetWeight, highTetWeight, strainMin, strainMax, collisionWeight, selfCollisionWeight, fixedWeight, periferalWeight);
		_ptp.setHookSutureWeights(hookWeight, sutureWeight, 0.3f);
		_surgAct->getSutures()->setAutoSutureSpacing(autoSutureSpacing);
		lsFname += collisionObject;  //  "collision_proxy_5_4_20.obj";
	}
	struct tetSubset {
		std::string objFile;
		float lowTetWeight;
		float highTetWeight;
		float strainMin;
		float strainMax;
//		std::list<std::string> objFiles;
	};
	std::list<tetSubset> tetSubsets;
	if ((oit = scnObj.find("tetrahedralSubsets")) != scnObj.end()) {
		json::Object tetSubObj = oit->second.ToObject();
		tetSubset ts;
		for (suboit = tetSubObj.begin(); suboit != tetSubObj.end(); ++suboit) {
			path = dataDirectory + suboit->first;
			ts.objFile = path;
			json::Object tetSubData = suboit->second.ToObject();
			for (auto dataoit = tetSubData.begin(); dataoit != tetSubData.end(); ++dataoit) {
				if (dataoit->first == "minStrain")
					ts.strainMin = dataoit->second.ToFloat();
				else if (dataoit->first == "maxStrain")
					ts.strainMax = dataoit->second.ToFloat();
				else if (dataoit->first == "lowTetWeight")
					ts.lowTetWeight = dataoit->second.ToFloat();
				else if (dataoit->first == "highTetWeight")
					ts.highTetWeight = dataoit->second.ToFloat();
				else;
			}
			tetSubsets.push_back(ts);
		}
	}
	else
		;

//	_mt->writeObjFile("objMat.obj");

	// COURT - put back in after cutter debug!
	createNewPhysicsLattice(maxDimSubdivs);  // now creating operable lattice on load
	_surgAct->getDeepCutPtr()->setMaterialTriangles(_mt);
	if (!_surgAct->getDeepCutPtr()->setDeepBed(_mt, deepBedFilepath.c_str(), &_vnTets)){
		_surgAct->sendUserMessage("Undermine layer .bed file could not be found-", "Error Message");
	}
	if (!tetSubsets.empty()) {
		for (auto& ts : tetSubsets)
			_tetSubsets.createSubset(&_vnTets, ts.objFile, ts.lowTetWeight, ts.highTetWeight, ts.strainMin, ts.strainMax);
		_tetSubsets.sendTetSubsets(&_vnTets, _mt, &_ptp);
	}

	_gl3w->frameScene(true);  // computes bounding spheres
	return true;
}

void bccTetScene::updateOldPhysicsLattice()
{  // at present regenerate.  May do incrementally later.
	try {
		remapTetPhysics rtp;
		rtp.getOldPhysicsData(&_vnTets);  // must be done before any new incisions
		if(!_tc.remakeVnTets(_mt))
			throw(std::logic_error("New recut of tet lattice failed.\n"));
		rtp.remapNewPhysicsNodePositions(&_vnTets);

#ifdef NO_PHYSICS
		_firstSpatialCoords.assign(_vnTets.nodeNumber(), Vec3f());
		_vnTets.setNodeSpatialCoordinatePointer(&_firstSpatialCoords[0]);  // for no physics debug
#else
		std::array<float, 3> *nodeSpatialCoords = _ptp.createBccTetStructure(_vnTets.getTetNodeArray(), (float)_vnTets.getTetUnitSize());
		_vnTets.setNodeSpatialCoordinatePointer(nodeSpatialCoords);  // vector created in _ptp
#endif
		rtp.restoreOldNodePositions(&_vnTets);
		_tetCol.initSoftCollisions(_mt, &_vnTets);  // after every topo change
		_tetSubsets.sendTetSubsets(&_vnTets, _mt, &_ptp);
		if (_forcesApplied) {  // _tetsModified not necessary as implied by calling this routine
			initPdPhysics();
			_tetsModified = true;
		}
		_physicsPaused = false;
	}  // end try block
	catch (char *e) {
		_surgAct->sendUserMessage(e, "Exception thrown:", false);
	}
}

void bccTetScene::createNewPhysicsLattice(int maximumDimensionSubdivisions)
{
//	std::chrono::time_point<std::chrono::system_clock> start, end;
//	start = std::chrono::system_clock::now();

	try {
		_tetsModified = false;
#ifdef _DEBUG
		_tc.makeFirstVnTets(_mt, &_vnTets, 30);  // 30
#else
		_tc.makeFirstVnTets(_mt, &_vnTets, maximumDimensionSubdivisions); // 70 for cleft scene with Eigen solver gives 36K tets in 0.36 seconds, 100 for cleft scene with oldMKL solver gives 102K tets, 20 for sphere test.
			// in release mode 180 subdivs gives 503,910 tets in 2.27 seconds without multithreading the tet cutter, 0.57 seconds multithreaded. 
#endif
		_surgAct->getDeepCutPtr()->setVnBccTetrahedra(&_vnTets);

		// COURT - could redo this only using nodes that are outside of boundary.  Revisit after periosteal undermining done.
		// get fixed scene boundary nodes
	//	for (int i = 0; i < (int)_fixedRegions.size(); ++i)
	//		setFixedRegion(_fixedRegions[i].corners);
	//	for (int i = 0; i < (int)_fixedMeshes.size(); ++i)
	//		setFixedMesh(&_fixedMeshes[i]);

#ifdef NO_PHYSICS
		_firstSpatialCoords.assign(_vnTets.nodeNumber(), Vec3f());
		_vnTets.setNodeSpatialCoordinatePointer(&_firstSpatialCoords[0]);  // for no physics debug
#else
		std::array<float, 3> *nodeSpatialCoords = _ptp.createBccTetStructure(_vnTets.getTetNodeArray(), (float)_vnTets.getTetUnitSize());
		_vnTets.setNodeSpatialCoordinatePointer(nodeSpatialCoords);  // vector created in _ptp
#endif

		_vnTets.materialCoordsToSpatialVector();

		//			end = std::chrono::system_clock::now();
//			std::chrono::duration<double> elapsed_seconds = end - start;
//			std::time_t end_time = std::chrono::system_clock::to_time_t(end);
//			std::string message("Physics initial creation took ");
//			message += std::to_string(elapsed_seconds.count());
//			message += " seconds for ";
//			message += std::to_string(_vnTets.tetNumber());
//			message += " tets.";
//			_surgAct->sendUserMessage(message.c_str(), "Timer");

		_tetsModified = false;
		_physicsPaused = false;
	}  // end try block
	catch (char *e) {
		_surgAct->sendUserMessage(e, "Exception thrown", false);
	}
}

void bccTetScene::initPdPhysics()
{
	fixPeriostealPeriferalVertices();
	if (!_tetCol.empty()) {
		_tetCol.updateFixedCollisions(_mt, &_vnTets);
	}
#ifndef NO_PHYSICS
	_surgAct->getHooks()->updateHookPhysics();
	_surgAct->getSutures()->updateSuturePhysics();
	_ptp.initializePhysics();
#endif
}

void bccTetScene::updatePhysics()
{
	if (_vnTets.empty())
		return;
	if (!_tetsModified && _forcesApplied) {
		_tetsModified = true;
		initPdPhysics();
	}

#ifndef NO_PHYSICS
	if (_tetsModified || _forcesApplied) {
		_tetCol.findSoftCollisionPairs();
		_ptp.solve();
	}
#endif

#ifdef WRITE_FOR_RENDER
	RenderHelper<float>::writeMesh(*_mt);
	RenderHelper<float>::frame++;
#endif
}
 
void bccTetScene::setVisability(char surface, char physics)
{  // 0=off, 1=on, 2=don't change
	if (surface < 1)
		_surgAct->getSurgGraphics()->getSceneNode()->visible = false;
	if (surface == 1)
		_surgAct->getSurgGraphics()->getSceneNode()->visible = true;
	if (physics < 1)
		_gl3w->getLines()->setLinesVisible(false);
	if (physics == 1) {
		if (!_gl3w->getLines()->getSceneNode()) {
			createTetLatticeDrawing();
			drawTetLattice();
		}
		else
			_gl3w->getLines()->setLinesVisible(true);
	}
}

void bccTetScene::updateSurfaceDraw()
{
	int nv;
	float *fp = _mt->getPositionArray(nv);
	for (int i = 0; i < nv; ++i) {
		if (_vnTets.getVertexTetrahedron(i) > -1)  // an excision may have occurred leaving an empty vertex
			_vnTets.getBarycentricTetPosition(_vnTets.getVertexTetrahedron(i), *(_vnTets.getVertexWeight(i)), (Vec3f &)*fp);
		fp += 3;
	}
	_surgAct->getSurgGraphics()->updatePositionsNormalsTangents();
	if (_gl3w->getLines()->linesVisible())
		drawTetLattice();
}

 void bccTetScene::fixPeriostealPeriferalVertices()
{  // this routine should be done only once after original lattice constructed
	int numElem, i, k;
	struct anchorPoint {
		bool isPeriferal;
		std::array<float, 3> baryWeight, pos;
	}ap;
	std::unordered_map<int, anchorPoint> fixPoints;  // key is tet index
	// fixed nodes take precedence over periferal nodes
	// Fix all nodes surrounding a periosteal or periferalvertex.
	materialTriangles::matTriangle *tr = _mt->getTriangleArray(numElem);
	ap.isPeriferal = false;
	for (i = 0; i < numElem; ++i) {
		if (tr[i].material == 7) {  // periosteal triangle
			for (k = 0; k < 3; ++k) {
				const Vec3f* vp = _vnTets.getVertexWeight(tr[i].v[k]);
				ap.baryWeight[0] = vp->X;
				ap.baryWeight[1] = vp->Y;
				ap.baryWeight[2] = vp->Z;
				_vnTets.vertexMaterialCoordinate(tr[i].v[k], ap.pos);
				fixPoints.insert(std::make_pair(_vnTets.getVertexTetrahedron(tr[i].v[k]), ap));
			}
		}
	}
	ap.isPeriferal = true;
	for (i = 0; i < numElem; ++i) {
		if (tr[i].material == 1) {  // soft border triangle
			for (k = 0; k < 3; ++k) {
				const Vec3f *vp = _vnTets.getVertexWeight(tr[i].v[k]);
				ap.baryWeight[0] = vp->X;
				ap.baryWeight[1] = vp->Y;
				ap.baryWeight[2] = vp->Z;
				_vnTets.vertexMaterialCoordinate(tr[i].v[k], ap.pos);
				fixPoints.insert(std::make_pair(_vnTets.getVertexTetrahedron(tr[i].v[k]), ap));
			}
		}
	}
	size_t n = fixPoints.size();
	std::vector<int> fixedTets, peripheralTets;
	fixedTets.reserve(n);
	peripheralTets.reserve(n);
	std::vector<std::array<float, 3> > fixedWeights, peripheralWeights, fixedPos, peripheralPos;
	fixedWeights.reserve(n);
	peripheralWeights.reserve(n);
	fixedPos.reserve(n);
	peripheralPos.reserve(n);
	for (auto &fp : fixPoints) {
		if (fp.second.isPeriferal) {
			peripheralTets.push_back(fp.first);
			peripheralWeights.push_back(fp.second.baryWeight);
			peripheralPos.push_back(fp.second.pos);
		}
		else {
			fixedTets.push_back(fp.first);
			fixedWeights.push_back(fp.second.baryWeight);
			fixedPos.push_back(fp.second.pos);
		}
	}
#ifndef NO_PHYSICS
	_ptp.setFixedVertices(fixedTets, fixedWeights, fixedPos, peripheralTets, peripheralWeights, peripheralPos);
#endif
}

void bccTetScene::setFixedRegion(const float(&corners)[6])
{  // AFTER lattice has been created, physics nodes inside this zone are fixed
	boundingBox<float> bb;
	bb.Empty_Box();
	bb.Enlarge_To_Include_Point((const float(&)[3])corners[0]);
	bb.Enlarge_To_Include_Point((const float(&)[3])corners[3]);
	Vec3f v;
	for (int n = (int)_vnTets.nodeNumber(), i = 0; i < n; ++i) {
		if (bb.Outside((const float(&)[3])v._v))
			_vnTets.setNodeFixationState(i, false);
		else
			_vnTets.setNodeFixationState(i, true);
	}
}

void bccTetScene::createTetLatticeDrawing()
{
	_nodeGraphicsPositions.clear();
	_nodeGraphicsPositions.assign(_vnTets.nodeNumber() << 2, 1.0f);
	GLfloat *ngp = &_nodeGraphicsPositions[0];
	for (int n = _vnTets.nodeNumber(), i = 0; i < n; ++i){
		float *fp = _vnTets.nodeSpatialCoordinatePtr(i);
		*(ngp++) = fp[0];
		*(ngp++) = fp[1];
		*(ngp++) = fp[2];
		++ngp;
	}
	std::set<std::pair<int, int> > segs;
	for (int n = _vnTets.tetNumber(), i=0; i<n; ++i){
		std::pair<int, int> ll;
		const int* tetNodes = _vnTets.tetNodes(i);
		for (int j = 0; j < 3; ++j){
			for (int k = j+1; k < 4; ++k){
				ll.first = std::min(tetNodes[j], tetNodes[k]);
				ll.second = std::max(tetNodes[j], tetNodes[k]);
				segs.insert(ll);
			}
		}
	}
	std::vector<GLuint> lines;
	lines.reserve(segs.size() * 3);
	for (auto &s : segs){
		lines.push_back(s.first);
		lines.push_back(s.second);
		lines.push_back(0xffffffff);
	}
	_gl3w->getLines()->setGl3wGraphics(_gl3w);
	float white2[4] = {1.0f, 1.0f, 1.0f, 1.0f};
	_gl3w->getLines()->addLines(_nodeGraphicsPositions, lines);
	_gl3w->getLines()->getSceneNode()->setColor(white2);
}

void bccTetScene::eraseTetLattice()
{
	_nodeGraphicsPositions.clear();
	_gl3w->getLines()->clear();
	_gl3w->getLines()->getSceneNode()->visible = false;
}

void bccTetScene::drawTetLattice()
{
	if (_nodeGraphicsPositions.empty())
		return;
	GLfloat *ngp = &_nodeGraphicsPositions[0];
	for (int n = _vnTets.nodeNumber(), i = 0; i < n; ++i){
		float *fp = _vnTets.nodeSpatialCoordinatePtr(i);
		*(ngp++) = fp[0];
		*(ngp++) = fp[1];
		*(ngp++) = fp[2];
		++ngp;
	}
	_gl3w->getLines()->updatePoints(_nodeGraphicsPositions);
}

bccTetScene::bccTetScene() : _physicsPaused(false), _forcesApplied(false), _tetsModified(false)
{
	_tetCol.setPdTetPhysics(&_ptp); // Qisi:set ptp for tetCol so things of ptp are accessible inside of tetCol
}


bccTetScene::~bccTetScene()
{
}
