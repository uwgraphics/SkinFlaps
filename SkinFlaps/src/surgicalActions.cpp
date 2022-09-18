#include "gl3wGraphics.h"
#include "Vec3f.h"
#include "Mat2x2f.h"
#include <sstream>
#include <fstream>
#include <exception>
#include <chrono>
#include <thread>
#include <assert.h>
#include "insidePolygon.h"
#include "prettyPrintJSON.h"
#include "surgGraphics.h"
#include <tbb/task_arena.h>
#include "FacialFlapsGui.h"
#include "surgicalActions.h"

// ReadyPileType ReadyPile;

surgicalActions::surgicalActions() : _toolState(0), _originalTriangleNumber(0), _sceneDir("0"), _historyDir("0"), _strongHooks(false), physicsDone(true), newTopology(false)
{
	_bts.setSurgicalActions(this);
	_historyArray.Clear();
	_historyIt = _historyArray.begin();
	_undermineTriangles.clear();
	_periostealUndermineTriangles.clear();
	_x=0.0f; _y=0.0f; _z=0.0f; _u=0.0f; _f=0.0f; _r=0.0f;
}

surgicalActions::~surgicalActions()
{
}

bool surgicalActions::saveSurgicalHistory(const char *fullFilePath)
{
	std::ofstream outf(fullFilePath);
	std::string ppStr,hstStr = Serialize(_historyArray);
	prettyPrintJSON pp;
	pp.convert(hstStr.c_str(), ppStr);
	outf.write(ppStr.c_str(), ppStr.size());
    outf.close();
	return true;
}

void surgicalActions::sendUserMessage(const char *message, const char *title, bool closeProgram)
{
	_ffg->sendUserMessage(message, title);
}

bool surgicalActions::rightMouseDown(std::string objectHit, float (&position)[3], int triangle)
{	// returns true if a surgical action is taken, false if this is a simple viewer command
	if((_toolState==0 && objectHit[1]!='_') || (_toolState>0 && (objectHit.substr(0,2)=="H_" || objectHit.substr(0,2)=="S_")))
		return false;
	// staticTriangle objects are only scenery. If user selects one, just ignore it. pick() should ignore it
	int hookNum=-1;
	if (_toolState != 7 && !_periostealUndermineTriangles.empty()){  // forgot to finish periosteal undermining so do it now
		int newTs = _toolState;
		_toolState = 7;
		onKeyDown(73);
		setToolState(newTs);
	}
	if (_toolState > 0) {  // active tool requested by user
		_bts.setPhysicsPause(true);  // stop doing physics updates
		// prevent user from doing a new op until previous one is finished
		while (!physicsDone)  // any previously enqueued physics thread must be complete before doing next op.
			;
	}
	if(_toolState==0)	//viewer
	{
		if(objectHit.substr(0,2)=="H_")	// user picked a hook
		{
			_selectedSurgObject = objectHit;
			hookNum = atoi(_selectedSurgObject.c_str()+2);
			_sutures.selectSuture(-1);
			_hooks.selectHook(hookNum);
		}
		else if(objectHit.substr(0,2)=="S_")	// user picked a suture
		{
			_selectedSurgObject = objectHit;
			hookNum = atoi(_selectedSurgObject.c_str()+2);
			_hooks.selectHook(-1);
			int userNum = _sutures.baseToUserSutureNumber(hookNum);
			if(userNum < 0){
				sendUserMessage("User can't select or delete an automatic suture.\nDelete a user applied suture before or after\nto remove the automatic row.", "Usage error", false);
				_sutures.selectSuture(-1);
			}
			else
				_sutures.selectSuture(hookNum);
		}
		else
			;
	}
	else if (_toolState == 1) {	// create hook mode
		auto sn = _gl3w->getNodePtr(objectHit);
		if (sn->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
			return false;
		materialTriangles* tr = _sg.getMaterialTriangles();
		float uv[2] = { 0.0f, 0.0f };
		tr->getBarycentricProjection(triangle, position, uv);
		int material;
		float hTx[2];
		Vec3f hVec;
		if (!setHistoryAttachPoint(triangle, uv, material, hTx, hVec))
			return false;  // try again
		if (_hooks.getNumberOfHooks() < 1) {	// initialize hooks
			_hooks.setHookSize(sn->getRadius() * 0.02f);
			_hooks.setShapes(_gl3w->getShapes());
			_hooks.setGLmatrices(_gl3w->getGLmatrices());
			_hooks.setPhysicsLattice(_bts.getPdTetPhysics_2());
			_hooks.setVnBccTetrahedra(_bts.getVirtualNodedBccTetrahedra());
			_hooks.setIncisions(&_incisions);
		}

		// COURT visual debug use
		/*		for (int j, i = 0; i < tr->numberOfTriangles(); ++i) {
					int* trp = tr->triangleVertices(i);
					for (j = 0; j < 3; ++j)
						if (trp[j] == 14272)
							break;
					if (j < 3) {
						triangle = i;
						if (j == 1)
							uv[0] = 1.0f;
						if (j == 2)
							uv[1] = 1.0f;
						_toolState = 0;
						break;
					}
				} */
		//		triangle = 19520;
		//		uv[0] = 0.33;
		//		uv[1] = 0.33;

		if ((hookNum = _hooks.addHook(tr, triangle, uv, _strongHooks)) > -1)
		{

			//			return true;  // for above debug use

			if (!_bts.getPdTetPhysics_2()->solverInitialized()) {  // solver must be initialized to add a hook
				_ffg->physicsDrag = true;
				_bts.setForcesAppliedFlag();
				physicsDone = false;
				tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
					_bts.updatePhysics();
					physicsDone = true;
					}
				);
			}

			_sutures.selectSuture(-1);
			_hooks.selectHook(hookNum);
			char s[80];
			sprintf(s, "H_%d", hookNum);
			_selectedSurgObject = s;
			if (_historyIt != _historyArray.end()) {
				json::Array tarr;
				for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
					tarr.push_back(*it);
				_historyArray.Clear();
				_historyArray = tarr;
			}
			json::Object hookObj, hookTitle;
			hookObj["hookNum"] = hookNum;
			hookObj["material"] = material;
			if (_strongHooks)
				hookObj["strongHook"] = true;
			json::Array vArr;
			vArr.push_back(hTx[0]);
			vArr.push_back(hTx[1]);
			hookObj["historyTexture"] = vArr;
			vArr.Clear();
			vArr.push_back(hVec[0]);
			vArr.push_back(hVec[1]);
			vArr.push_back(hVec[2]);
			hookObj["displacement"] = vArr;
			hookTitle["addHook"] = hookObj;
			_historyArray.push_back(hookTitle);
			_historyIt = _historyArray.end();
			// don't _frame->setToolState(0) or will get unnecessary hook move on mouse up or motion.  Fix there.
		}
		_bts.setPhysicsPause(false);
	}
	else if (_toolState == 2)	// incision mode
	{
		auto sn = _gl3w->getNodePtr(objectHit);
		if (sn->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
			return false;
		materialTriangles* tr = _sg.getMaterialTriangles();
		if (!_fence.isInitialized()) {	// initialize fence
			_fence.setFenceSize(sn->getRadius() * 0.02f);
			_fence.setGl3wGraphics(_gl3w);
		}
		bool endConn = false;
		Vec3f vtx(position), nrm;
		float uv[2] = { 0.0f,0.0f };  // pos[3], norm[3],
		auto ToutToFirstPoint = [&]() {
			_fence.getPostPos(0, vtx);
			_fence.getPostNormal(0, nrm);
			_fence.addPost(tr, _fence.getPostTriangle(0), vtx._v, nrm._v, endConn, false, false);
			_hooks.selectHook(-1);
			_sutures.selectSuture(-1);
			onKeyDown(GLFW_KEY_ENTER);	// press enter key for user
		};
		if (_ffg->CtrlOrShiftKeyIsDown()) {
			endConn = true;
			int edg, oldTriangle = triangle;
			float param, closeIncisionDistance = _incisions.closestSkinIncisionPoint(vtx, triangle, edg, param);
			if (closeIncisionDistance < FLT_MAX) {
				if (_fence.numberOfPosts() > 2) {  // possible Tout to the first incision point
					_fence.getPostPos(0, vtx);
					float len = (vtx - Vec3f(position)).length();
					if (len < closeIncisionDistance) {
						ToutToFirstPoint();
						return true;
					}
				}
			}
			else {
				if (_fence.numberOfPosts() < 1) {
					sendUserMessage("There is no existing skin incision edge to T in to.", "Usage error", false);
					return false;
				}
				else if (_fence.numberOfPosts() < 3) {
					sendUserMessage("Can only self T out if there are three incision points already present.  Try again-", "Usage error", false);
					return false;
				}
				else {  // can only Tout to the first point
					ToutToFirstPoint();
					return true;
				}
			}
			position[0] = vtx.X; position[1] = vtx.Y; position[2] = vtx.Z;
			if (edg < 1)
				uv[0] = param;
			else if (edg > 1)
				uv[1] = 1.0f - param;
			else {
				uv[0] = 1.0f - param;
				uv[1] = param;
			}
		}
		else
			tr->getBarycentricProjection(triangle, position, uv);
		if (tr->triangleMaterial(triangle) != 2) {
			sendUserMessage("With this tool you can only incise from top side of skin.", "USER ERROR");
			return true;
		}
		tr->getBarycentricPosition(triangle, uv, vtx._v);
		tr->getBarycentricNormal(triangle, uv, nrm._v);
		_fence.addPost(tr, triangle, vtx._v, nrm._v, endConn, false, false);
		_hooks.selectHook(-1);
		_sutures.selectSuture(-1);
		if (_fence.numberOfPosts() > 1 && endConn)	// this must finish an incision
			onKeyDown(GLFW_KEY_ENTER);	// press enter key for user
	}
	else if (_toolState == 3){	// start undermine tool
		auto sn = _gl3w->getNodePtr(objectHit);
		if (sn->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
			return false;
		materialTriangles* tr = _sg.getMaterialTriangles();
		if (tr->triangleMaterial(triangle) != 2 && tr->triangleMaterial(triangle) != 10) {
			sendUserMessage("With this tool you can only undermine from top side of skin.", "USER ERROR");
			return true;
		}
		undermineTriangle ut;
		ut.triangle = triangle;
		ut.incisionConnect = !_ffg->CtrlOrShiftKeyIsDown();
		_undermineTriangles.push_back(ut);
		_bts.updateSurfaceDraw();
		if (!_incisions.addUndermineTriangle(triangle, 2, ut.incisionConnect)) {
			// ignore false return as it does no harm
			_undermineTriangles.pop_back();
		}
	}
	else if (_toolState == 4){	// create suture mode
		auto sn = _gl3w->getNodePtr(objectHit);
		if (sn->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
			return false;
		materialTriangles* tr = _sg.getMaterialTriangles();
		if (_sutures.getNumberOfSutures() < 1) {	// initialize sutures
			_sutures.setSutureSize(sn->getRadius()*0.003f);
			_sutures.setShapes(_gl3w->getShapes());
			_sutures.setGLmatrices(_gl3w->getGLmatrices());
			_sutures.setPhysicsLattice(_bts.getPdTetPhysics_2());
			_sutures.setVnBccTetrahedra(_bts.getVirtualNodedBccTetrahedra());
			_sutures.setDeepCut(&_incisions);
			_sutures.setSurgicalActions(this);
		}
		int i = 0, edg, triMat = tr->triangleMaterial(triangle);
		int eTri = triangle;
		float param, uv[2];
		tr->getBarycentricProjection(triangle, position, uv);
		if (triMat == 2) {
			_sutures.nearestSkinIncisionEdge(uv, eTri, edg, param);
			if (edg < 1) {
				uv[0] = param;
				uv[1] = 0.0f;
			}
			else if (edg > 1){
				uv[1] = 1.0f - param;
				uv[0] = 0.0f;
			}
			else {
				uv[0] = 1.0f - param;
				uv[1] = param;
			}
		}
		else if (triMat == 3) {
			int aTE = tr->triAdjs(triangle)[0];
			if (tr->triangleMaterial(aTE >> 2) > 3 && tr->triangleMaterial(aTE >> 2) < 7)
				aTE = tr->triAdjs(triangle - 1)[0];  // incision convention
			else
				assert(tr->triangleMaterial(aTE >> 2) == 2);
			eTri = aTE >> 2;
			edg = aTE & 3;
			Vec3f V0, V1, P(position);
			int* tp = tr->triangleVertices(eTri);
			tr->getVertexCoordinate(tp[edg], V0._v);
			tr->getVertexCoordinate(tp[(edg + 1) % 3], V1._v);
			V1 -= V0;
			float denom = (V1 * V1);
			if (denom < 1e-8f)
				param = 0.0;
			else {
				param = (P - V0) * V1 / denom;
				if (param >= 1.0f)
					param = 1.0f;
				if (param < 0.0f)
					param = 0.0f;
			}
			if (edg < 1) {
				uv[0] = param;
				uv[1] = 0.0f;
			}
			else if (edg > 1) {
				uv[1] = 1.0f - param;
				uv[0] = 0.0f;
			}
			else {
				uv[0] = 1.0f - param;
				uv[1] = param;
			}
			triMat = 2;  // corrected
		}
		else if (triMat == 6) {
			sendUserMessage("Currently can't suture within cut muscle or fat. Please suture above or below ths point-", "USER ERROR");
			return true;
		}
		else {  // suturing to a non-skin surface object.
			if (uv[0] + uv[1] > 0.67f) {  // force to an edge
				edg = 1;
				param = uv[1] / (uv[0] + uv[1]);
			}
			else if (uv[0] > uv[1]) {
				edg = 0;
				param = uv[0];
			}
			else {
				edg = 2;
				param = 1.0f - uv[1];
			}
		}
		tr->getBarycentricPosition(eTri, uv, _dragXyz);
		i = _sutures.addUserSuture(tr, eTri, edg, param);
		if (_ffg->CtrlOrShiftKeyIsDown()) {
			int prevMat = _sutures.previousUserSuture(i);
			if (prevMat > -1)
				prevMat = _sutures.firstVertexMaterial(prevMat);
			if (prevMat != 2 || triMat != 2) {
sendUserMessage("Can only create an automatic suture line on a skin/mucosal edges-", "USER ERROR");
_sutures.deleteSuture(i);
return true;
			}
			else
			_sutures.setLinked(i, true);
		}
		else
		_sutures.setLinked(i, false);
		_sutures.setSecondVertexPosition(i, _dragXyz);
		char s[10];
		sprintf(s, "S_%d", i);
		_selectedSurgObject = s;
	}
	else if (_toolState == 5)	// excise mode
	{
	auto sn = _gl3w->getNodePtr(objectHit);
	if (sn->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
		return false;
	materialTriangles* tr = _sg.getMaterialTriangles();
	int mat = tr->triangleMaterial(triangle);
	if (mat == 3 || mat == 6) {
		sendUserMessage("Can't excise from a skin/mucosal edge or a cut muscle belly,  Try again-", "USER ERROR");
		return true;
	}
	float uv[2], hTx[2];
	tr->getBarycentricProjection(triangle, position, uv);
	int material;
	Vec3f hVec;
	if (!setHistoryAttachPoint(triangle, uv, material, hTx, hVec))
		return true;
	if (_historyIt != _historyArray.end()) {
		json::Array tarr;
		for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
			tarr.push_back(*it);
		_historyArray.Clear();
		_historyArray = tarr;
	}
	json::Object exciseObj, exciseTitle;
	exciseObj["material"] = material;
	json::Array vArr;
	vArr.push_back(hTx[0]);
	vArr.push_back(hTx[1]);
	exciseObj["historyTexture"] = vArr;
	vArr.Clear();
	vArr.push_back(hVec[0]);
	vArr.push_back(hVec[1]);
	vArr.push_back(hVec[2]);
	exciseObj["displacement"] = vArr;
	exciseTitle["excise"] = exciseObj;
	_historyArray.push_back(exciseTitle);
	_historyIt = _historyArray.end();
	_incisions.excise(triangle);
	physicsDone = false;
	//		_ffg->physicsDrag = true;
	tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
		_bts.updateOldPhysicsLattice();
		newTopology = true;
		physicsDone = true;
		}
	);
	_bts.setPhysicsPause(false);
	_hooks.selectHook(-1);
	_sutures.selectSuture(-1);
	_selectedSurgObject = "";
	_ffg->setToolState(0);
	setToolState(0);
	}
	else if (_toolState == 6)	// deep cut mode
	{
		if (objectHit.substr(0, 3) == "NP_") {	// user picked a fence handle
			_selectedSurgObject = objectHit;
			hookNum = atoi(_selectedSurgObject.c_str() + 3);
			_sutures.selectSuture(-1);
			_hooks.selectHook(-1);
			_fence.selectPost(hookNum);
			Vec3f pos(position[0], position[1], position[2]);
			_fence.setSpherePos(hookNum, pos);
			return true;
		}
		auto sn = _gl3w->getNodePtr(objectHit);
		if (sn->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
			return false;
		materialTriangles* tr = _sg.getMaterialTriangles();
		if (tr->triangleMaterial(triangle) != 2 && tr->triangleMaterial(triangle) != 5) {
			sendUserMessage("Can only deep cut from unelevated skin top or deep bed.  Try again-", "USER ERROR");
			return true;
		}
		if (!_fence.isInitialized()) {	// initialize fence
			_fence.setFenceSize(sn->getRadius() * 0.02f);
			_fence.setGl3wGraphics(_gl3w);
		}
		bool closedEnd = true;
		if (_ffg->CtrlOrShiftKeyIsDown())
			closedEnd = false;
		Vec3f norm;
		float pos[3], uv[2] = { 0.0f, 0.0f };
		tr->getBarycentricProjection(triangle, position, uv);
		// if triangle selected is material 2 which has been undermined, xRay through it to its corresponding deep bed triangle
		if (tr->triangleMaterial(triangle) == 2 && _incisions.triangleUndermined(triangle)) {
			float tx[2];
			tr->getBarycentricTexture(triangle, uv, tx);
			Vec3f displ(0.0f, 0.0f, 0.0f);
			if (!getHistoryAttachPoint(5, tx, displ, triangle, uv, false)) {
				std::string msg = "Can't Xray through top to a deep bed location.";
				historyAttachFailure(msg);
				return false;
			}
		}
		tr->getBarycentricPosition(triangle, uv, pos);
		if(_fence.numberOfPosts() < 1)
			tr->getTriangleNormal(triangle, norm._v, true);
		else
			_fence.getPostNormal(_fence.numberOfPosts() - 1, norm);  // use previous normal as starting point to minimize potential crossover
		_fence.addPost(tr, triangle, pos, norm._v, false, true, !closedEnd);  // never connect to nearest hard edge
		hookNum = _fence.numberOfPosts() - 1;
		_fence.selectPost(hookNum);
		char s[80];
		sprintf(s, "NP_%d", hookNum);
		_selectedSurgObject = s;
		_hooks.selectHook(-1);
		_sutures.selectSuture(-1);
	}
	else if (_toolState == 7){	// periosteal undermine mode
		auto sn = _gl3w->getNodePtr(objectHit);
		if (sn->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
			return false;
		materialTriangles* tr = _sg.getMaterialTriangles();
		Vec3f cameraPos, dir;
		_gl3w->getTrianglePickLine(cameraPos._v, dir._v);  // this routine only used here as of 3/22/2022
		_bts.updateSurfaceDraw();
		perioTri pt;
		pt.incisionConnect = !_ffg->CtrlOrShiftKeyIsDown();
		pt.periostealTriangle = _incisions.addPeriostealUndermineTriangle(triangle, dir, pt.incisionConnect);
		if (pt.periostealTriangle > 0x7ffffffe){
			sendUserMessage("No periosteal triangle hit.  Try again-", "USER ERROR");
			return true;
		}
		_periostealUndermineTriangles.push_back(pt);
	}
	else if (_toolState == 8)	// creates a collision proxy suture.  COURT NUKE AFTER SOFT-SOFT COLLISIONS ADDED!
	{
		auto sn = _gl3w->getNodePtr(objectHit);
		if (sn->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
			return false;
		materialTriangles* tr = _sg.getMaterialTriangles();
		if (_sutures.getNumberOfSutures() < 1) {	// initialize sutures
			_sutures.setSutureSize(sn->getRadius()*0.003f);
			_sutures.setShapes(_gl3w->getShapes());
			_sutures.setGLmatrices(_gl3w->getGLmatrices());
			_sutures.setPhysicsLattice(_bts.getPdTetPhysics_2());
			_sutures.setVnBccTetrahedra(_bts.getVirtualNodedBccTetrahedra());
			_sutures.setDeepCut(&_incisions);
		}
		int material = tr->triangleMaterial(triangle), edg;
		if (material != 2 && material != 5) {
			sendUserMessage("Can only add a collisionProxy suture to material 4 (via 2) or 5-", "USER ERROR");
			_ffg->setToolState(0);
			setToolState(0);
		}
		float uv4[2], uv5[2];
		int tri4, tri5;
		if (material == 2) {
			tr->closestPoint(position, tri4, uv4, 4);
			tr->closestPoint(position, tri5, uv5, 5);
		}
		else {
			tr->getBarycentricProjection(triangle, position, uv5);
			tri5 = triangle;
			tr->closestPoint(position, tri4, uv4, 4);
		}
		float hTx4[2], hTx5[2], param, xyz[3];
		Vec3f hVec;
		if (!setHistoryAttachPoint(tri5, uv5, material, hTx5, hVec))
			return false;  // try again
		if (!setHistoryAttachPoint(tri4, uv4, material, hTx4, hVec))
			return false;  // try again
		if (material != 4) {
			sendUserMessage("Error in collisionProxy suture to material 4-", "USER ERROR");
			_ffg->setToolState(0);
			setToolState(0);
			_bts.setPhysicsPause(false);
		}
		auto triEdge = [&](float(&uv)[2]) {
			if (uv[0] + uv[1] > 0.67f) {  // force to an edge
				edg = 1;
				param = uv[1] / (uv[0] + uv[1]);
			}
			else if (uv[0] > uv[1]) {
				edg = 0;
				param = uv[0];
			}
			else {
				edg = 2;
				param = 1.0f - uv[1];
			}
		};
		triEdge(uv4);
		int i = _sutures.addUserSuture(tr, tri4, edg, param);
		_sutures.setLinked(i, false);
		triEdge(uv5);
		_sutures.setSecondEdge(i, tr, tri5, edg, param);
		tr->getBarycentricPosition(tri5, uv5, xyz);
		_sutures.setSecondVertexPosition(i, xyz);
		json::Array hArr;
		json::Object pObj, sutureTitle;
		int sNum = _sutures.baseToUserSutureNumber(i);
		_bts.setPhysicsPause(false);
		if (sNum < 0) {
			sendUserMessage("Dr. Cutting error-", "PROGRAM ERROR");
			return false;
		}
		pObj["collisionSutureNum"] = sNum;
		hArr.Clear();
		hArr.push_back(hTx4[0]);
		hArr.push_back(hTx4[1]);
		pObj["mat4Texture"] = hArr;
		hArr.Clear();
		hArr.push_back(hTx5[0]);
		hArr.push_back(hTx5[1]);
		pObj["mat5Texture"] = hArr;
		sutureTitle["collisionProxySuture"] = pObj;
		_historyArray.push_back(sutureTitle);
		_historyIt = _historyArray.end();
		_hooks.selectHook(-1);
		_sutures.selectSuture(i);
		_ffg->setToolState(0);
		setToolState(0);
	}
	else
		;
	return true;
}

bool surgicalActions::rightMouseUp(std::string objectHit, float (&position)[3], int triangle)
{  // this routine only called when terminating a surgical drag operation
	std::string hStr;
	if((_toolState==2 || _toolState==0) && _selectedSurgObject.substr(0,2)=="P_")	// fence post selected in viewer or incision mode
		_selectedSurgObject = "";
	else if (_toolState == 4)	{	// finish applying a suture
		// prevent user from doing a new op until previous one is finished
		assert(physicsDone);  // physics update thread must be complete before doing next op.
		assert(_selectedSurgObject.substr(0,2)=="S_");
		materialTriangles *tr = NULL;
		int i = atoi(_selectedSurgObject.c_str()+2);
		surgGraphics *sg = NULL;
		if (objectHit != "") {
			auto sn = _gl3w->getNodePtr(objectHit);
			if (sn->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
				return false;
			tr = _sg.getMaterialTriangles();
		}
		int eTri = triangle;
		int edge, triMat = tr->triangleMaterial(triangle);
		float param, uv[2];
		auto invalidate = [&]() {
			// prevent user from doing a new op until previous one is finished
			if (!physicsDone)  // physics update thread must be complete before doing next op.
				throw(std::logic_error("Trying to invalidate a suture while a physics thread is active.\n"));
			_sutures.deleteSuture(i);
			_bts.setPhysicsPause(false);
			_selectedSurgObject = "";
			_hooks.selectHook(-1);
			_sutures.selectSuture(-1);
			setToolState(0);
			_ffg->setToolState(0);
		};
		if (tr == NULL){
			invalidate();
			return true;
		}
		if (triMat == 2){
			if (_sutures.firstVertexMaterial(i) != 2){
				sendUserMessage("A skin/mucosal edge can only be sutured to another skin/mucosal edge-", "USER ERROR");
				invalidate();
				return true;
			}
			tr->getBarycentricProjection(triangle, position, uv);
			_sutures.nearestSkinIncisionEdge(uv, eTri, edge, param);
			if (edge < 1) {
				uv[0] = param;
				uv[1] = 0.0f;
			}
			else if (edge > 1) {
				uv[1] = 1.0f - param;
				uv[0] = 0.0f;
			}
			else {
				uv[0] = 1.0f - param;
				uv[1] = param;
			}
			if (eTri < 0) {
				sendUserMessage("Couldn't find a suitable edge point for this suture.  Try again-", "PROGRAM ERROR");
				invalidate();
				return true;
			}
		}
		else if (triMat == 3) {
			int aTE = tr->triAdjs(triangle)[0];
			if (tr->triangleMaterial(aTE >> 2) > 3 && tr->triangleMaterial(aTE >> 2) < 7) 
				aTE = tr->triAdjs(triangle - 1)[0];  // incision convention
			else
				assert(tr->triangleMaterial(aTE >> 2) == 2);
			eTri = aTE >> 2;
			edge = aTE & 3;
			Vec3f V0, V1, P(position);
			int* tp = tr->triangleVertices(eTri);
			tr->getVertexCoordinate(tp[edge], V0._v);
			tr->getVertexCoordinate(tp[(edge + 1) % 3], V1._v);
			V1 -= V0;
			float denom = (V1 * V1);
			if (denom < 1e-8f)
				param = 0.0;
			else {
				param = (P - V0) * V1 / denom;
				if (param >= 1.0f)
					param = 1.0f;
				if (param < 0.0f)
					param = 0.0f;
			}
			if (edge < 1) {
				uv[0] = param;
				uv[1] = 0.0f;
			}
			else if (edge > 1) {
				uv[1] = 1.0f - param;
				uv[0] = 0.0f;
			}
			else {
				uv[0] = 1.0f - param;
				uv[1] = param;
			}
		}
		else if (triMat == 6){
			sendUserMessage("Currently can't suture in the middle of cut muscle belly or fat. Please suture above or below this point-", "USER ERROR");
			invalidate();
			return true;
		}
		else{
			if (_sutures.firstVertexMaterial(i) == 2){
				sendUserMessage("A deep tissue can only be sutured to deep tissue and not a skin/mucosal edge-", "USER ERROR");
				invalidate();
				return true;
			}
			tr->getBarycentricProjection(triangle, position, uv);
			// suturing to a non-skin surface object.
			if (uv[0] + uv[1] > 0.67f){  // force to an edge
				edge = 1;
				param = uv[1] / (uv[0] + uv[1]);
			}
			else if (uv[0] > uv[1]){
				edge = 0;
				param = uv[0];
			}
			else{
				edge = 0;
				param = 1.0f - uv[1];
			}
		}
		int sRet = _sutures.setSecondEdge(i, tr, eTri, edge, param);
		_bts.setForcesAppliedFlag();
		if (sRet < 1){
			float pos[3], uv[2] = {0.0f, 0.0f};
			if (edge < 1)
				uv[0] = param;
			else if (edge > 1)
				uv[1] = 1.0f - param;
			else{
				uv[0] = 1.0f - param;
				uv[1] = param;
			}
			tr->getBarycentricPosition(eTri, uv, pos);
			if (!physicsDone)  // physics update thread must be complete before doing next op.
				throw(std::logic_error("Trying to add a suture while a physics thread is active.\n"));
			_sutures.setSecondVertexPosition(i, pos);
			if (_sutures.isLinked(i)) {
				physicsDone = false;
				_ffg->physicsDrag = true;
				tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {
					_sutures.laySutureLine(i);
					physicsDone = true;
					}
				);
			}
		}
		else if (sRet < 2){
			sendUserMessage("Trying to suture to same side of incision is not allowed-", "USER ERROR");
			invalidate();
			return false;
		}
		else
			assert(false);
		_bts.setPhysicsPause(false);

		if (_historyIt != _historyArray.end()) {
			json::Array tarr;
			for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
				tarr.push_back(*it);
			_historyArray.Clear();
			_historyArray = tarr;
		}
		auto getSutureUv = [&]() {
			param += (param < 0.002f) ? 0.001f : -0.001f;
			if (edge < 1) {
				uv[0] = param;
				uv[1] = 0.001f;
			}
			else if (edge > 1) {
				uv[0] = 0.001f;
				uv[1] = 1.0f - param;
			}
			else {
				uv[1] = param;
				uv[0] = 0.999f - param;
			}
		};
		getSutureUv();
		int material;
		float hTx[2];
		Vec3f hVec;
		if (!setHistoryAttachPoint(eTri, uv, material, hTx, hVec))
			invalidate();
		json::Array hArr;
		json::Object pObj, sutureTitle;
		int sNum = _sutures.baseToUserSutureNumber(i);
		if(sNum < 0) {
			sendUserMessage("Trying to finish a non-user suture manually-", "PROGRAM ERROR");
			invalidate();
			return false;
		}
		pObj["sutureNum"] = sNum;
		pObj["linked"] = _sutures.isLinked(i);
		pObj["material1"] = material;  // should probably allow multiple material suturing
		hArr.Clear();
		hArr.push_back(hTx[0]);
		hArr.push_back(hTx[1]);
		pObj["historyTexture1"] = hArr;
		hArr.Clear();
		hArr.push_back(hVec[0]);
		hArr.push_back(hVec[1]);
		hArr.push_back(hVec[2]);
		pObj["displacement1"] = hArr;
		int tri;
		_sutures.getEdgeAttachment(i, true, tr, tri, edge, param);
		getSutureUv();
		if (!setHistoryAttachPoint(tri, uv, material, hTx, hVec))
			invalidate;
		pObj["material0"] = material;
		hArr.Clear();
		hArr.push_back(hTx[0]);
		hArr.push_back(hTx[1]);
		pObj["historyTexture0"] = hArr;
		hArr.Clear();
		hArr.push_back(hVec[0]);
		hArr.push_back(hVec[1]);
		hArr.push_back(hVec[2]);
		pObj["displacement0"] = hArr;
		sutureTitle["addSuture"] = pObj;
		_historyArray.push_back(sutureTitle);
		_historyIt = _historyArray.end();
		_hooks.selectHook(-1);
		_sutures.selectSuture(i);
		_ffg->setToolState(0);
		_bts.setPhysicsPause(false);
		setToolState(0);
	}
	else if (_selectedSurgObject.substr(0, 2) == "H_")	// hook selected. Can only drag hooks.
	{
		if (_toolState == 1) {
			_bts.setPhysicsPause(false);
			setToolState(0);
			_ffg->setToolState(0);
			return true;
		}
		if (_toolState == 0) {  // Too many spurius hook moves recorded due to zoom releases. Fixed in cleftSimViewer.
			Vec3f xyz, selXyz;
			int hookNum = atoi(_selectedSurgObject.c_str() + 2);
			_hooks.getHookPosition(hookNum, xyz._v);
			_hooks.getSelectPosition(hookNum, selXyz._v);
			selXyz -= xyz;
//			if (selXyz.length2() < 0.01f)  // ignore small movements to unclutter history file
//				return true;
			if (_historyIt != _historyArray.end()) {
				json::Array tarr;
				for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
					tarr.push_back(*it);
				_historyArray.Clear();
				_historyArray = tarr;
			}
			json::Array hArr;
			hArr.push_back(hookNum);
			hArr.push_back((double)xyz._v[0]);
			hArr.push_back((double)xyz._v[1]);
			hArr.push_back((double)xyz._v[2]);
			json::Object mObj;
			mObj["moveHook"] = hArr;
			_historyArray.push_back(mObj);
			_historyIt = _historyArray.end();
			setToolState(0);
		}
	}
	else if (_toolState == 6){
		if (_selectedSurgObject.substr(0, 3) == "NP_") {	// fence post selected in viewer or incision mode
			int postNum = atoi(_selectedSurgObject.c_str() + 3);
			_fence.selectPost(postNum);
//			if (postNum > 0 && !_incisions.inputCorrectFence(&_fence)) {  / changed to correcting fence problems on <enter> key
//				_ffg->sendUserMessage("This post does not connect to the previous one-", "Try again");
//				_fence.deleteLastPost();
//				_incisions.popLastDeepPost();
//				_selectedSurgObject = "";
//			}
		}
	}
	else
		return false;
	return true;
}

bool surgicalActions::mouseMotion(float dScreenX, float dScreenY)
{
	Vec3f xyz, dv;
	if(_toolState==6 && _selectedSurgObject.substr(0,3)=="NP_")
	{
		int postNum = atoi(_selectedSurgObject.c_str()+3);
		_fence.getSpherePos(postNum, xyz);
		_gl3w->getGLmatrices()->getDragVector(dScreenX, dScreenY, xyz._v, dv._v);
		xyz += dv;
		_fence.setSpherePos(postNum, xyz);
	}
	else if(_toolState==4)	{
		assert(_selectedSurgObject.substr(0,2)=="S_");
		int sutNum = atoi(_selectedSurgObject.c_str()+2);
		_gl3w->getGLmatrices()->getDragVector(dScreenX,dScreenY,_dragXyz,dv._v);
		_dragXyz[0]+=dv._v[0]; _dragXyz[1]+=dv._v[1]; _dragXyz[2]+=dv._v[2];
		const float *mm=_gl3w->getGLmatrices()->getFrameAndRotationMatrix();
		transformVector3(_dragXyz,mm,xyz._v);
		xyz *= 0.7f;
		xyz._v[0]-=mm[12]; xyz._v[1]-=mm[13]; xyz._v[2]-=mm[14];
		dv._v[0] = mm[0]*xyz._v[0] + mm[1]*xyz._v[1] + mm[2]*xyz._v[2];
		dv._v[1] = mm[4]*xyz._v[0] + mm[5]*xyz._v[1] + mm[6]*xyz._v[2];
		dv._v[2] = mm[8]*xyz._v[0] + mm[9]*xyz._v[1] + mm[10]*xyz._v[2];
		_sutures.setSecondVertexPosition(sutNum, dv._v);
	}
	else if (_toolState != 1 && _selectedSurgObject.substr(0, 2) == "H_")	// hook selected.
	{
		int hookNum = atoi(_selectedSurgObject.c_str() + 2);
		_hooks.getHookPosition(hookNum, xyz._v);
		_gl3w->getGLmatrices()->getDragVector(dScreenX, dScreenY, xyz._v, dv._v);
		xyz += dv;
		_bts.setForcesAppliedFlag();  // this is a hook move so forces are applied
		if (!_bts.isPhysicsPaused() || !physicsDone) {
			_bts.setPhysicsPause(true);  // stop doing physics updates
			// prevent user from doing a new op until previous one is finished
			while (!physicsDone)  // any previously enqueued physics thread must be complete before doing next op.
				;
		}
		_hooks.setHookPosition(hookNum, xyz._v);
	}
	else
		;
	return true;
}

void surgicalActions::onKeyDown(int key)
{
	std::string hStr;
	// ctrl and shift keys now handled by frame calls
	if(key == GLFW_KEY_DELETE)	// delete key
	{
///		// can't delete periosteal undermines (toolState 7) already done.
//		if (_toolState == 7) {
//			sendUserMessage("Sorry. Periosteal undermining can't be undone-", "USER ERROR");
//			return;
//		}
//		else 
		if (_toolState == 2) {
			if (_fence.numberOfPosts() > 0)
				_fence.deleteLastPost();
			return;  // don't reset toolstate
		}
		else if (_toolState == 6) {
			if (_fence.numberOfPosts() > 0) {
				_fence.deleteLastPost();
				_incisions.popLastDeepPost();
				_selectedSurgObject = "";
			}
			return;  // don't reset toolstate
		}
		else if (_toolState == 3){
			_incisions.clearCurrentUndermine(2);
			_undermineTriangles.clear();
			return;
		}
		else if (_selectedSurgObject.substr(0, 2) == "H_")
		{
			int hookNum = atoi(_selectedSurgObject.c_str()+2);
			// prevent user from doing a new op until previous one is finished
			_bts.setPhysicsPause(true);  // don't spawn another physics update till complete
			while (!physicsDone)  // physics update thread must be complete before doing next op.
				;
			_hooks.deleteHook(hookNum);
			if (_historyIt != _historyArray.end()) {
				json::Array tarr;
				for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
					tarr.push_back(*it);
				_historyArray.Clear();
				_historyArray = tarr;
			}
			json::Object dObj;
			dObj["deleteHook"] = hookNum;
			_historyArray.push_back(dObj);
			_historyIt = _historyArray.end();
		}
		else if(_selectedSurgObject.substr(0,2)=="S_")
		{
			if (_historyIt != _historyArray.end()) {
				json::Array tarr;
				for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
					tarr.push_back(*it);
				_historyArray.Clear();
				_historyArray = tarr;
			}
			json::Object sObj;
			int sutNum = atoi(_selectedSurgObject.c_str() + 2);
			int userNum = _sutures.baseToUserSutureNumber(sutNum);
			_bts.setPhysicsPause(true);  // don't spawn another physics update till complete
			// prevent user from doing a new op until previous one is finished
			while (!physicsDone)  // physics update thread must be complete before doing next op.
				;
			int linkNum = _sutures.deleteSuture(sutNum);
			if (userNum < 0) {
				json::Object lObj;
				lObj["autoSuturesFor"] = _sutures.baseToUserSutureNumber(linkNum);
				sObj["deleteSuture"] = lObj;
			}
			else
				sObj["deleteSuture"] = userNum;
			_historyArray.push_back(sObj);
			_historyIt = _historyArray.end();
		}
		else
			;
		_ffg->setToolState(0);
		setToolState(0);
		_bts.setPhysicsPause(false);
	}
	else if (key == GLFW_KEY_ENTER)	// <enter> key
	{
		// prevent user from doing a new op until previous one is finished
		assert(physicsDone);  // physics update thread must be complete before doing next op.
		if (_toolState == 7){	//periosteal undermine mode
			_bts.setPhysicsPause(true);  // should already be done
			while (!physicsDone)
				;
			materialTriangles *mt = _sg.getMaterialTriangles();
			for (int n = mt->numberOfTriangles(), i = 0; i < n; ++i){
				if (mt->triangleMaterial(i) == 10)
					mt->setTriangleMaterial(i, 8);  // 8 is a periosteal triangle that has been undermined
			}
			_bts.updateSurfaceDraw();
			while (!physicsDone)  // physics update thread must be complete before doing next op.
				;
			physicsDone = false;
			_ffg->physicsDrag = true;
			tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
				_bts.fixPeriostealPeriferalVertices();
				_bts.nonTetPhysicsUpdate();
				newTopology = true;
				physicsDone = true;
				}
			);
			_bts.setPhysicsPause(false);
			if (_historyIt != _historyArray.end()) {
				json::Array tarr;
				for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
					tarr.push_back(*it);
				_historyArray.Clear();
				_historyArray = tarr;
			}
			float hTx[2], uv[2] = { 0.333f, 0.333f };
			int material;
			Vec3f hVec;
			json::Array uArr;
			json::Object uObj, pObj;
			auto ptit = _periostealUndermineTriangles.begin();
			while (ptit != _periostealUndermineTriangles.end()) {
				uObj.Clear();
				setHistoryAttachPoint(ptit->periostealTriangle, uv, material, hTx, hVec);
				uObj["material"] = material;
				json::Array sArr;
				sArr.push_back(hTx[0]);
				sArr.push_back(hTx[1]);
				uObj["historyTexture"] = sArr;
				sArr.Clear();
				sArr.push_back(hVec[0]);
				sArr.push_back(hVec[1]);
				sArr.push_back(hVec[2]);
				uObj["displacement"] = sArr;
				uObj["incisionConnect"] = (bool)ptit->incisionConnect;
				pObj["periostealTriangle"] = uObj;
				uArr.push_back(pObj);
				++ptit;
			}
			uObj.Clear();
			uObj["periostealUndermine"] = uArr;
			_historyArray.push_back(uObj);
			_incisions.clearCurrentUndermine(8);  // set all periosteal undermined triangles to material 8 and reset.
			_periostealUndermineTriangles.clear();
			_historyIt = _historyArray.end();
			_hooks.selectHook(-1);
			_sutures.selectSuture(-1);
			_selectedSurgObject = "";
		}
		else if (_toolState == 2)	//incision mode
		{
			std::vector<Vec3f> positions, normals;
			std::vector<float> postUvs;
			std::vector<int> postTriangles;
			bool edgeStart = false, edgeEnd = false, Tout = false, nukeThis = false, sOpen, eOpen;
			int n = _fence.getPostData(positions, normals, postTriangles, postUvs, edgeStart, edgeEnd, sOpen, eOpen);
			if (_historyIt != _historyArray.end()) {
				json::Array tarr;
				for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
					tarr.push_back(*it);
				_historyArray.Clear();
				_historyArray = tarr;
			}
			json::Object iObj;
			iObj["incisedObject"] = 0;	// for now only one object incisable
			iObj["Tin"] = edgeStart;
			iObj["Tout"] = edgeEnd;
			iObj["pointNumber"] = n;
			json::Array iArr;
			iArr.push_back(iObj);
			materialTriangles *tri=_sg.getMaterialTriangles();
			float uv[2];  //  , minParam = 1.0e15f;
			int material;
			float hTx[2];
			Vec3f hVec;
			json::Array hArr;
			for (int i = 0; i<n; ++i)	{
				uv[0] = postUvs[i << 1];
				uv[1] = postUvs[(i << 1) + 1];
				setHistoryAttachPoint(postTriangles[i], uv, material, hTx, hVec);
				json::Object pObj;
				pObj["material"] = material;
				hArr.Clear();
				hArr.push_back(hTx[0]);
				hArr.push_back(hTx[1]);
				pObj["historyTexture"] = hArr;
				hArr.Clear();
				hArr.push_back(hVec[0]);
				hArr.push_back(hVec[1]);
				hArr.push_back(hVec[2]);
				pObj["displacement"] = hArr;
				iObj.Clear();
				iObj["incisionPoint"] = pObj;
				iArr.push_back(iObj);
			}
			iObj.Clear();
			iObj["makeIncision"] = iArr;
			_historyArray.push_back(iObj);
			_historyIt = _historyArray.end();
			if(!nukeThis)	{
				if (!_incisions.skinCut(positions, normals, edgeStart, edgeEnd)) {
						sendUserMessage("Incision tool error.  Please save history file for debugging-", "Error Message");
				}
				else {
					if (_incisions.physicsRecutRequired()){
						_bts.setPhysicsPause(true);  // don't spawn another physics update till complete
						while (!physicsDone)  // physics update thread must be complete before doing next op.
							;
						physicsDone = false;
						_ffg->physicsDrag = true;
						tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
							_bts.updateOldPhysicsLattice();
							newTopology = true;
							physicsDone = true;
							}
						);
					}
					else {
						newTopology = true;
					}
				}
			}
			_bts.setPhysicsPause(false);
			_fence.clear();
		}
		else if (_toolState == 3) {	// undermine mode
			if (_historyIt != _historyArray.end()) {
				json::Array tarr;
				for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
					tarr.push_back(*it);
				_historyArray.Clear();
				_historyArray = tarr;
			}
			float hTx[2], uv[2] = {0.333f, 0.333f};
			int material;
			Vec3f hVec;
			json::Array uArr;
			json::Object uObj, pObj;
			auto uit = _undermineTriangles.begin();
			while (uit != _undermineTriangles.end()) {
				uObj.Clear();
				setHistoryAttachPoint(uit->triangle, uv, material, hTx, hVec);
				uObj["material"] = 2;  // at time executed all set to 10, but they came in as 2
				uObj["incisionConnect"] = (bool)uit->incisionConnect;
				json::Array sArr;
				sArr.push_back(hTx[0]);
				sArr.push_back(hTx[1]);
				uObj["historyTexture"] = sArr;
				sArr.Clear();
				sArr.push_back(hVec[0]);
				sArr.push_back(hVec[1]);
				sArr.push_back(hVec[2]);
				uObj["displacement"] = sArr;
				pObj["underminePoint"] = uObj;
				uArr.push_back(pObj);
				++uit;
			}
			uObj.Clear();
			uObj["undermine"] = uArr;
			_historyArray.push_back(uObj);
			_historyIt = _historyArray.end();
			_bts.setPhysicsPause(true);  // should already be done
			while (!physicsDone)
				;
			_bts.updateSurfaceDraw();
			_incisions.undermineSkin();
			_undermineTriangles.clear();
			physicsDone = false;
			_ffg->physicsDrag = true;
			tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
				_bts.updateOldPhysicsLattice();
				newTopology = true;
				physicsDone = true;
				}
			);
			_bts.setPhysicsPause(false);
		}
		else if (_toolState == 6)	// deep cut mode
		{
			if (!_incisions.inputCorrectFence(&_fence, _ffg)) {
				return;
			}
			std::vector<Vec3f> positions, rays;
			std::vector<float> postUvs;
			std::vector<int> postTriangles;
			bool edgeStart, edgeEnd, startOpen, endOpen;  //  , Tout = false, nukeThis = false;
			if (_historyArray.size()>0 && _historyIt != _historyArray.end()) {
				json::Array tarr;
				for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
					tarr.push_back(*it);
				_historyArray.Clear();
				_historyArray = tarr;
			}
			int n = _fence.getPostData(positions, rays, postTriangles, postUvs, edgeStart, edgeEnd, startOpen, endOpen); // bools not relevant
			materialTriangles *tri = _sg.getMaterialTriangles();
			float hTx[2], uv[2];
			int material;
			Vec3f hVec;
			json::Array iArr;
			json::Object iObj, dObj;
			iObj["deepCutObject"] = 0;	// for now only one object incisable
			iObj["openIn"] = startOpen;
			iObj["openOut"] = endOpen;
			iObj["pointNumber"] = n;
			iArr.push_back(iObj);
			json::Array pArr;
			for (int i = 0; i<n; ++i)	{
				iObj.Clear();
				dObj.Clear();
				uv[0] = postUvs[i << 1];
				uv[1] = postUvs[(i << 1) + 1];
				setHistoryAttachPoint(postTriangles[i], uv, material, hTx, hVec);
				dObj["material"] = material;
				json::Array sArr;
				sArr.push_back(hTx[0]);
				sArr.push_back(hTx[1]);
				dObj["historyTexture"] = sArr;
				sArr.Clear();
				sArr.push_back(hVec[0]);
				sArr.push_back(hVec[1]);
				sArr.push_back(hVec[2]);
				dObj["displacement"] = sArr;
				sArr.Clear();
				sArr.push_back(rays[i].x());
				sArr.push_back(rays[i].y());
				sArr.push_back(rays[i].z());
				dObj["postNormal"] = sArr;
				iObj["deepCutPoint"] = dObj;
				iArr.push_back(iObj);
			}
			iObj.Clear();
			iObj["makeDeepCut"] = iArr;
			_historyArray.push_back(iObj);
			_historyIt = _historyArray.end();
			if (!_bts.isPhysicsPaused())
				throw(std::logic_error("Physics must be paused before deep cut."));
			while (!physicsDone)
				;
			_bts.updateSurfaceDraw();
			if (!_incisions.cutDeep()) {
				sendUserMessage("Attempted deepCut failed. Save history to debug.", "PROGRAM ERROR");
				return;
			}
			physicsDone = false;
			_ffg->physicsDrag = true;
			_ffg->user_message_flag = false;
//			tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
				_bts.updateOldPhysicsLattice();
				newTopology = true;
				physicsDone = true;
//				}
//			);
			_fence.clear();
			_incisions.clearDeepCutter();
			_bts.setPhysicsPause(false);
		}
		else
			;
		_ffg->setToolState(0);
		setToolState(0);
	}
	else
		;
}

void surgicalActions::onKeyUp(int key)
{  // ctrl and shift key now handled by frame call
}

bool surgicalActions::loadScene(const char *modelDirectory, const char *sceneFilename)
{
	bool ret = _bts.loadScene(modelDirectory, sceneFilename);  // computes bounding spheres
	_sceneDir.assign(modelDirectory);
	_originalTriangleNumber = _sg.getMaterialTriangles()->numberOfTriangles();
	if(ret && _historyArray.size() < 1) {
		std::string dstr(modelDirectory),fstr(sceneFilename);
		_historyArray.Clear();
		std::size_t n;
		while ((n = dstr.find("\\")) < dstr.npos)
			dstr.replace(n, 1, "/");
		json::Object loadObj;
		loadObj["loadSceneFile"] = fstr;
		_historyArray.push_back(loadObj);
		_historyIt = _historyArray.end();
	}
	_gl3w->zeroViewRotations();
	return ret;
}

bool surgicalActions::setHistoryAttachPoint(const int triangle, const float(&uv)[2], int &material, float(&historyTexture)[2], Vec3f &historyVec)
{  // Input an attach point in current environment. Outputs a historyTriangle, historyUv, and historyVec for storage in a history file.
	// This attachment point is created by program with variable physics state at time of incisions.  For this reason move away a safe distance to an original triangle and use
	// historyVec to find closest original location.  historyVec is in material coords so less sensitive to physics state.
	materialTriangles *mtp = _sg.getMaterialTriangles();
	material = mtp->triangleMaterial(triangle);
	if (material == 3 || material == 6) {
		sendUserMessage("Can't attach to side of skin incision or middle of cut muscle. Try again-", "USER ERROR", false);
		return false;
	}
	// attachments on an edge must be moved inside triangle
	float tp[2];
	if (uv[0] < 1e-5f) {
		if (uv[1] < 1e-5f) {
			tp[0] = 0.001f;
			tp[1] = 0.001f;
		}
		else {
			tp[0] = 0.001f;
			tp[1] = uv[1] * 0.996f;
		}
	}
	else if (uv[1] < 1e-5f) {
		tp[0] = uv[0] * 0.996f;
		tp[1] = 0.001f;
	}
	else if (uv[0] + uv[1] > 0.998f) {
		tp[0] = uv[0] *0.996f;
		tp[1] = uv[1] * 0.996f;
	}
	else {
		tp[0] = uv[0];
		tp[1] = uv[1];
	}
	auto isBorderTriangle = [mtp, material](int tri) ->bool {  // on incision edge?
		unsigned int *adjs;
		adjs = mtp->triAdjs(tri);
		for (int i = 0; i < 3; ++i) {
			int mat = mtp->triangleMaterial(adjs[i] >> 2);
			if ((mat > 2 && mat < 4) || mat == 6)  // hard intermaterial cut edge to move away from
				return true;
		}
		return false;
	};
	historyVec.set(0.0f, 0.0f, 0.0f);
	if (!isBorderTriangle(triangle)) {
		Vec2f tx;
		int *tr = mtp->triangleTextures(triangle);
		float *fp = mtp->getTexture(tr[0]);
		tx.set(fp[0], fp[1]);
		tx *= 1.0f - tp[0] - tp[1];
		fp = mtp->getTexture(tr[1]);
		tx += Vec2f(fp[0], fp[1]) * tp[0];
		fp = mtp->getTexture(tr[2]);
		tx += Vec2f(fp[0], fp[1]) * tp[1];
		historyTexture[0] = tx[0];
		historyTexture[1] = tx[1];
		return true;
	}
	auto vbt = _bts.getVirtualNodedBccTetrahedra();
	auto triangleNormal = [mtp, vbt](const int tri) -> const Vec3f {
		Vec3f N, vM[3];
		int *tr = mtp->triangleVertices(tri);
		for (int j = 0; j < 3; ++j)
			vbt->vertexGridLocus(tr[j], vM[j]);
		vM[1] -= vM[0];
		vM[2] -= vM[0];
		N = vM[1] ^ vM[2];
		N.normalize();
		return N;
	};
	// get material coord displacement direction
	Vec3f planeN, edgeN;
	unsigned int *adjs = mtp->triAdjs(triangle);
	int eNum = 0;
	edgeN.set(0.0f, 0.0f, 0.0f);
	for (int i = 0; i < 3; ++i) {
		int mat = mtp->triangleMaterial(adjs[i] >> 2);
		if ((mat > 2 && mat < 4) || mat == 6) {  // hard intermaterial cut edge to move away from
			edgeN += triangleNormal(adjs[i] >> 2);
			++eNum;
		}
	}
	assert(eNum < 3 && eNum>0);
	if(eNum > 1)
		edgeN.normalize();
	planeN = edgeN ^ triangleNormal(triangle);
	int *tr = mtp->triangleVertices(triangle);
	vbt->vertexGridLocus(tr[0], historyVec);
	historyVec *= 1.0f - tp[0] - tp[1];
	for (int i = 0; i<2; ++i){
		Vec3f v;
		vbt->vertexGridLocus(tr[i+1], v);
		historyVec += v*tp[i];
	}
	float d = planeN * historyVec, tetSizeSq = (float)vbt->getTetUnitSize();
	tetSizeSq *= tetSizeSq;
	int nextTri = triangle;
	int lastEdge = -1;
	do {
		tr = mtp->triangleVertices(nextTri);
		Vec3f now, nextV;  // last, 
		vbt->vertexGridLocus(tr[0], nextV);
		float dNext = planeN * nextV - d;
		for (int i = 2; i > -1; --i) {
			vbt->vertexGridLocus(tr[i], now);
			float dNow = planeN * now - d;
			if (i != lastEdge && std::signbit(dNext) != std::signbit(dNow)) {
				Vec3f vI;
				vI = nextV * dNow + now * -dNext;
				vI /= dNow - dNext;
				vI -= historyVec;
				if (vI*edgeN < 0.0) {
					lastEdge = (i + 2) % 3;
					// COURT - what size to use? Could be a physics/model specific value based on variability in number of iterations to stability.
					// creator of a history file will usually have lots of physics iterations before applying a suture, but someone playing back history quickly may have very few.
					if (vI.length2()*tetSizeSq > 0.0001f && !isBorderTriangle(nextTri)) {
						// pull inside this triangle to ensure an insideTest() texture find on retrieval
//						float edgeParam = -dLast / (dNow - dLast);
						float edgeParam = -dNext / (dNow - dNext);
						if (lastEdge < 1) {
							if (edgeParam < 0.005f) {
								tp[0] = 0.001f;
								tp[1] = 0.001f;
							}
							else if (edgeParam > 0.995f) {
								tp[0] = 0.998f;
								tp[1] = 0.001f;
							}
							else {
								tp[0] = edgeParam * 0.996f;
								tp[1] = 0.001f;
							}
						}
						 else if (lastEdge < 2) {
							if (edgeParam < 0.005f) {
								tp[0] = 0.998f;
								tp[1] = 0.001f;
							}
							else if (edgeParam > 0.995f) {
								tp[1] = 0.998f;
								tp[0] = 0.001f;
							}
							else {
								tp[1] = edgeParam * 0.996f;
								tp[0] = (1.0f - edgeParam) * 0.996f;
							}
						}
						 else {
							if (edgeParam > 0.998f) {
								tp[0] = 0.001f;
								tp[1] = 0.001f;
							}
							else if (edgeParam < .005f) {
								tp[1] = 0.998f;
								tp[0] = 0.001f;
							}
							else {
								tp[1] = (1.0f - edgeParam) * 0.996f;
								tp[0] = 0.001f;
							}

						}
						int* triTx = mtp->triangleTextures(nextTri);
						Vec2f tx;
						float *fp = mtp->getTexture(triTx[0]);
						tx.set(fp[0], fp[1]);
						tx *= 1.0f - tp[0] - tp[1];
						fp = mtp->getTexture(triTx[1]);
						tx += Vec2f(fp[0], fp[1]) * tp[0];
						fp = mtp->getTexture(triTx[2]);
						tx += Vec2f(fp[0], fp[1]) * tp[1];
						historyVec = -vI;
						historyTexture[0] = tx[0];
						historyTexture[1] = tx[1];
						return true;
					}
					else {
						unsigned int adj = mtp->triAdjs(nextTri)[lastEdge];
						nextTri = adj >> 2;
//						lastEdge = ((adj & 3) + 1) %3;
						lastEdge = adj & 3;
						break;
					}
				}
			}
			if (i > 1)
				material = 40;  // send error
			nextV = now;
			dNext = dNow;
		}
	} while (mtp->triangleMaterial(nextTri) == material);
	// have just failed to find a suitable attach point using a search perpendicular to the wound edge.  This is usually due to trying to suture
	// to a flap corner point where that corner is an acute angle.  Instead will look for the nearest non-border triangle.
	material = mtp->triangleMaterial(triangle);
	struct nearTri {
		bool borderTriangle;
		Vec3f centroid;
		int triangle;
	} nt;
	nt.borderTriangle = false;
	nt.triangle = triangle;
	std::multimap<float, nearTri> nearTris;
	nearTris.insert(std::make_pair(0.0f, nt));
	auto ntit = nearTris.begin();
	std::set<int> trisDone;
	while (ntit != nearTris.end()) {
		ntit->second.borderTriangle = false;
		for (int i = 0; i < 3; ++i) {
			unsigned int adj = mtp->triAdjs(ntit->second.triangle)[i];
			int mat = mtp->triangleMaterial(adj >> 2);
			if ((mat > 2 && mat < 4) || mat == 6) {  // hard intermaterial cut edge to move away from
				ntit->second.borderTriangle = true;
				continue;
			}
			else if (!trisDone.insert(adj >> 2).second)
				continue;
			else {
				nt.triangle = adj >> 2;
				int *tr = mtp->triangleVertices(adj >> 2);
				Vec3f v;
				nt.centroid = { 0.0f, 0.0f, 0.0f };
				for (int j = 0; j < 3; ++j) {
					vbt->vertexGridLocus(tr[j], v);
					nt.centroid += v * 0.33333f;
				}
				nearTris.insert(std::make_pair((historyVec - nt.centroid).length2(), nt));
			}
		}
		if (!ntit->second.borderTriangle) {
			historyVec -= ntit->second.centroid;
			int *ttx = mtp->triangleTextures(ntit->second.triangle);
			Vec2f tx, vt;
			tx.set(0.0f, 0.0f);
			for (int i = 0; i < 3; ++i) {
				mtp->getTexture(ttx[i], vt._v);
				tx += vt * 0.33333f;
			}
			historyTexture[0] = tx[0];
			historyTexture[1] = tx[1];
			return true;
		}
		else {
			nearTris.erase(ntit);
			ntit = nearTris.begin();
		}
	}
	sendUserMessage("No suitable attachment triangle for this point. Try again-", "PROGRAM LIMITATION", false);
	return false;
}

void surgicalActions::historyAttachFailure(std::string& errorDescription) {
	json::Array tarr;
	for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
		tarr.push_back(*it);
	_historyArray.Clear();
	_historyArray = tarr;
	std::string msg = errorDescription;
	msg.append("\nSetting history back one step and truncating further forward.");
	sendUserMessage(msg.c_str(), "Program error");
}

bool surgicalActions::getHistoryAttachPoint(const int material, const float(&historyTexture)[2], const Vec3f &displacement, int &triangle, float(&uv)[2], bool findEdge)
{  // Input a history attach point from history file. Outputs a triangle, and parametric uv coord in current environment.
	materialTriangles *mtp = _sg.getMaterialTriangles();
	std::vector<Vec2f> triTex;
	triTex.assign(3, Vec2f());
	Vec2f txIn(historyTexture[0], historyTexture[1]);
	int k, n = mtp->numberOfTriangles(), matIn = material;
	if (matIn == 8)  // make all periosteal materials 7
		matIn = 7;
	insidePolygon ip;
	for (k = 0; k < n; ++k) {
		if (material > 6) {
			if (mtp->triangleMaterial(k) < 7)  // in an undermine periosteum may have already been labelled as 7, 8, or 10.
				continue;
		}
		else {
			if (mtp->triangleMaterial(k) != material && mtp->triangleMaterial(k) != 10)  // in an undermine may already have been labelled as 10
				continue;
		}
		int *tr = mtp->triangleTextures(k);
		float *fp;
		for (int j = 0; j < 3; ++j) {
			fp = mtp->getTexture(tr[j]);
			triTex[j].set(fp[0], fp[1]);
		}
		if (ip.insidePolygon2f(txIn, triTex)) {  // COURT texture seams may screw this up. See following backup strategy
			Mat2x2f M;
			M.Initialize_With_Column_Vectors(triTex[1] - triTex[0], triTex[2] - triTex[0]);
			Vec2f R = M.Robust_Solve_Linear_System(txIn - triTex[0]);
			uv[0] = R.X;
			uv[1] = R.Y;
			break;
		}
	}
	if (k >= n) {  // this section to handle texture seam case
		triangle = -1;
		float dsq, minDsq = FLT_MAX;
		for (k = 0; k < n; ++k) {
			if (material > 6) {
				if (mtp->triangleMaterial(k) < 7)  // in an undermine periosteum may have already been labelled as 7, 8, or 10.
					continue;
			}
			else {
				if (mtp->triangleMaterial(k) != material && mtp->triangleMaterial(k) != 10)  // in an undermine may already have been labelled as 10
					continue;
			}
			int* tr = mtp->triangleTextures(k);
			float* fp;
			for (int j = 0; j < 3; ++j) {
				fp = mtp->getTexture(tr[j]);
				triTex[0].set(fp[0], fp[1]);
				dsq = (triTex[0] - txIn).length2();
				if (minDsq > dsq) {
					minDsq = dsq;
					triangle = k;
					uv[0] = j==1 ? 1.0f : 0.0f;
					uv[1] = j > 1 ? 1.0f : 0.0f;
				}
			}
		}
		k = triangle;
		if (triangle < 0)
			throw(std::logic_error("Program error in finding history attach point."));
		return true;
	}
	if (displacement[0] == 0.0f && displacement[1] == 0.0f &&displacement[2] == 0.0f) {
		triangle = k;
		return true;
	}
	// attachments on an edge must be moved inside
	if (uv[0] < 1e-5f) {
		if (uv[1] < 1e-5f) {
			uv[0] = 0.001f;
			uv[1] = 0.001f;
		}
		else
			uv[0] = 0.001f;
	}
	else if (uv[1] < 1e-5f)
		uv[1] = 0.001f;
	else if (uv[0] + uv[1] > 0.998f) {
		uv[0] = uv[0] * 0.996f;
		uv[1] = uv[1] * 0.996f;
	}
	else
		;
	Vec3f triV[3], N, V, startV;
	// switch to material coords to minimize effect of varying physics state
	auto vbt = _bts.getVirtualNodedBccTetrahedra();
	int *tr = mtp->triangleVertices(k);
	for (int j = 0; j < 3; ++j)
		vbt->vertexGridLocus(tr[j], triV[j]);
	startV = triV[0] * (1.0f - uv[0] - uv[1]) + triV[1] * uv[0] + triV[2] * uv[1];
	V = (triV[1] - triV[0])^(triV[2] - triV[0]);
	N = displacement ^ V;
	float d = N * startV, dsqFinal = displacement.length2();
	int nextTri = k;
	int lastEdge = -1;
	do {
		tr = mtp->triangleVertices(nextTri);
		Vec3f now, nextV;  // last, 
		vbt->vertexGridLocus(tr[0], nextV);
		float dNext = N * nextV - d;
		int i;
		for (i = 2; i > -1; --i) {
			vbt->vertexGridLocus(tr[i], now);
			float dNow = N * now - d;
			if (i != lastEdge && std::signbit(dNext) != std::signbit(dNow)) {
				Vec3f vI;
				vI = nextV * dNow + now * -dNext;
				vI /= dNow - dNext;
				vI -= startV;
				if (vI*displacement > 0.0) {
					lastEdge = i;
					if (!findEdge && vI.length2() > dsqFinal) {  // overshoot displacement
						// pull inside this triangle to ensure an insideTest() texture find on retrieval
						for (int j = 0; j < 3; ++j)
							vbt->vertexGridLocus(tr[j], triV[j]);
						triV[1] -= triV[0];
						triV[2] -= triV[0];
						V = startV + displacement - triV[0];
						Mat2x2f M2;
						M2.Initialize_With_Column_Vectors(Vec2f(triV[1]*triV[1], triV[1]*triV[2]), Vec2f(triV[2]*triV[1], triV[2]*triV[2]));
						Vec2f R, B = {V*triV[1], V*triV[2]};
						R = M2.Robust_Solve_Linear_System(B);
						uv[0] = R[0];
						uv[1] = R[1];
						triangle = nextTri;
						return true;
					}
					else {
						unsigned int adj = mtp->triAdjs(nextTri)[lastEdge];
						int nextMaterial = mtp->triangleMaterial(adj>>2);
						if (nextMaterial == 8)
							nextMaterial = 7;
						if (nextMaterial != matIn && nextMaterial != 10) {  // crossed incision edge
							triangle = nextTri;
							float edgeParam = dNow / (dNow - dNext);
							if (lastEdge < 1) {
								uv[0] = edgeParam;
								uv[1] = 0.0f;
							}
							else if (lastEdge < 2) {
								uv[1] = edgeParam;
								uv[0] = 1.0f - edgeParam;
							}
							else {
								uv[1] = 1.0f - edgeParam;
								uv[0] = 0.0f;
							}
							return true;
						}
						nextTri = adj >> 2;
						lastEdge = adj & 3;
						break;
					}
				}
			}
			nextV = now;
			dNext = dNow;
		}
		if (i < 0) {  // no way out
			triangle = -1;
			return false;
		}
	} while (true);
	return false;
}

bool surgicalActions::loadHistory(const char *historyDir, const char *historyFile)
{
	if (_historyArray.begin() != _historyArray.end())
		return false;
	_historyDir.assign(historyDir);
	// set scene dir elsewhere. Don't lock to history location.
	std::size_t found = _historyDir.rfind("History");
	if (found == _historyDir.size())
		sendUserMessage("History directory specified incorrectly.", "Program error", false);
	_historyArray.Clear();
	std::string hPath(_historyDir);
	hPath.append(historyFile);
	std::ifstream is(hPath.c_str());
	if(!is.is_open())
		return false;
	std::stringstream buffer;
	buffer << is.rdbuf();
	std::string str = buffer.str();
	json::Value hstData = json::Deserialize(str);
	if(hstData.GetType() != json::ArrayVal)
		return false;
	_historyArray = hstData.ToArray();
	_historyIt = _historyArray.begin();
	nextHistoryAction();  // loads scene in history file
	return true;
}

void surgicalActions::promoteFakeSutures()
{
	if (_historyIt != _historyArray.end()) {
		json::Array tarr;
		for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
			tarr.push_back(*it);
		_historyArray.Clear();
		_historyArray = tarr;
	}
	json::Object title;
	title["promoteSutureApproximations"] = 0;
	_historyArray.push_back(title);
	_historyIt = _historyArray.end();
	_bts.promoteSutures();
}

void surgicalActions::pausePhysics()
{
	if (_historyIt != _historyArray.end()) {
		json::Array tarr;
		for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
			tarr.push_back(*it);
		_historyArray.Clear();
		_historyArray = tarr;
	}
	json::Object title;
	title["pausePhysics"] = 0;
	_historyArray.push_back(title);
	_historyIt = _historyArray.end();
	_bts.setPhysicsPause(true);
}

void surgicalActions::nextHistoryAction()
{
	try {
		if (_historyIt == _historyArray.end())
		{
			sendUserMessage("There are no more actions found in this history file-", "SURGICAL HISTORY INFORMATION", false);
			return;
		}
		_bts.setPhysicsPause(true);  // don't spawn another physics update till complete
		// prevent user from doing a new op until previous one is finished
		while (!physicsDone)  // physics update thread must be complete before doing next op.
			;
		_gl3w->drawAll();
		if (_historyIt->HasKey("loadSceneFile"))
		{
			const json::Object& fObj = _historyIt->ToObject();
			if (!loadScene(_sceneDir.c_str(), fObj.begin()->second.ToString().c_str())) {
				sendUserMessage("The scene file in the history file can't be loaded-", "SURGICAL HISTORY INFORMATION", false);
				_historyArray.Clear();
			}
			else {
				_ffg->setModelFile(fObj.begin()->second.ToString());
				++_historyIt;
			}
		}
		else if (_historyIt->HasKey("addHook"))
		{
			materialTriangles *tr = _sg.getMaterialTriangles();
			if (tr == NULL)
				return;
			json::Object hookObj = (*_historyIt)["addHook"].ToObject();
			int hookNum, triangle;
			int material;
			float uv[2], historyTx[2];
			assert(hookObj.HasKey("material"));
			material = hookObj["material"].ToInt();
			hookNum = hookObj["hookNum"].ToInt();
			assert(hookObj.HasKey("historyTexture"));
			json::Array vArr = hookObj["historyTexture"];
			historyTx[0] = vArr[0];
			historyTx[1] = vArr[1];
			Vec3f V;
			vArr.Clear();
			vArr = hookObj["displacement"];
			V[0] = vArr[0].ToFloat();
			V[1] = vArr[1].ToFloat();
			V[2] = vArr[2].ToFloat();
			if (!getHistoryAttachPoint(material, historyTx, V, triangle, uv, false)) {
				std::string msg = "History file attachment failure at hook number ";
				msg.append(std::to_string(hookNum));
				historyAttachFailure(msg);
				return;
			}
			if (_hooks.getNumberOfHooks() < 1) {	// initialize hooks
				_hooks.setHookSize(_sg.getSceneNode()->getRadius()*0.02f);
				_hooks.setShapes(_gl3w->getShapes());
				_hooks.setGLmatrices(_gl3w->getGLmatrices());
				_hooks.setPhysicsLattice(_bts.getPdTetPhysics_2());
				_hooks.setVnBccTetrahedra(_bts.getVirtualNodedBccTetrahedra());
				_hooks.setIncisions(&_incisions);
			}
			hookNum = -1;
			bool strongHook = false;
			if(hookObj.HasKey("strongHook"))
				strongHook = true;
			int newHookNum;
			if ((newHookNum = _hooks.addHook(tr, triangle, uv, strongHook)) > -1)
			{
				if (!_bts.getPdTetPhysics_2()->solverInitialized()) {  // solver must be initialized to add a hook. Done once.
					_bts.setForcesAppliedFlag();
					physicsDone = false;
					_ffg->physicsDrag = true;
					tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
						_bts.updatePhysics();
						physicsDone = true;
						}
					);
				}
				_sutures.selectSuture(-1);
				_hooks.selectHook(hookNum);
				char s[20];
#ifdef _WINDOWS
				sprintf_s(s, 19, "H_%d", hookNum);
#else
				sprintf(s, "H_%d", hookNum);
#endif
				_selectedSurgObject = s;
			}
			++_historyIt;
		}
		else if (_historyIt->HasKey("moveHook"))
		{
			Vec3f xyz;
			int hookNum;
			const json::Array& hArr = (*_historyIt)["moveHook"];
			hookNum = hArr[0].ToInt();
			xyz._v[0] = hArr[1].ToFloat();
			xyz._v[1] = hArr[2].ToFloat();
			xyz._v[2] = hArr[3].ToFloat();
			_hooks.setHookPosition(hookNum, xyz._v);
			_bts.setForcesAppliedFlag();
			char s[20];
			sprintf(s, "H_%d", hookNum);
			_sutures.selectSuture(-1);
			_hooks.selectHook(hookNum);	// deselect hooks
			_selectedSurgObject = s;
			++_historyIt;
		}
		else if (_historyIt->HasKey("deleteHook"))
		{
			int hookNum = (*_historyIt)["deleteHook"].ToInt();
			_hooks.deleteHook(hookNum);
			_sutures.selectSuture(-1);
			_hooks.selectHook(-1);
			++_historyIt;
		}
		else if (_historyIt->HasKey("makeIncision"))
		{
			json::Array iArr = (*_historyIt)["makeIncision"].ToArray();
			json::Object iObj = iArr[0].ToObject();
			int i, incisPointNum;
			bool startIncis, endIncis;
			// for now ignore "incisedObject"] as there is only one incisable object.  May change this later.
			startIncis = iObj["Tin"].ToBool();
			endIncis = iObj["Tout"].ToBool();
			incisPointNum = iObj["pointNumber"].ToInt();
			std::vector<Vec3f> positions, normals;
			positions.assign(incisPointNum, Vec3f());
			normals.assign(incisPointNum, Vec3f());
			materialTriangles *mtp = _sg.getMaterialTriangles();
			for (i = 0; i < incisPointNum; ++i)
			{
				iObj = iArr[i + 1].ToObject();
				if (!iObj.HasKey("incisionPoint")) {
					sendUserMessage("There is an error in this history file.  Truncating from this point forward-", "", false);
					if (_historyIt != _historyArray.end()) {
						json::Array tarr;
						for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
							tarr.push_back(*it);
						_historyArray.Clear();
						_historyArray = tarr;
						_historyIt = _historyArray.end();
					}
					_bts.setPhysicsPause(false);
					return;
				}
				json::Object ipObj = iObj["incisionPoint"].ToObject();
				int tri;
				int material;
				float hTx[2], uv[2];
				material = ipObj["material"].ToInt();
				json::Array sArr = ipObj["historyTexture"].ToArray();
				hTx[0] = sArr[0].ToFloat();
				hTx[1] = sArr[1].ToFloat();
				Vec3f hV;
				sArr.Clear();
				sArr = ipObj["displacement"];
				hV[0] = sArr[0].ToFloat();
				hV[1] = sArr[1].ToFloat();
				hV[2] = sArr[2].ToFloat();
				if (!getHistoryAttachPoint(material, hTx, hV, tri, uv, false)) {
					std::string msg = "History file incision point location failure.";
					historyAttachFailure(msg);
					return;
				}
				assert(mtp->triangleMaterial(tri) == 2);
				mtp->getBarycentricPosition(tri, uv, positions[i]._v);
				mtp->getBarycentricNormal(tri, uv, normals[i]._v);
			}
			if (!_incisions.skinCut(positions, normals, startIncis, endIncis)) {
				sendUserMessage("Incision in history file failed.", "Program error", false);
			}
			else {
				if (_incisions.physicsRecutRequired()) {
					physicsDone = false;
					_ffg->physicsDrag = true;
					tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
						_bts.updateOldPhysicsLattice();
						newTopology = true;
						physicsDone = true;
						}
					);
				}
				else
					newTopology = true;
			}
			++_historyIt;
		}
		else if (_historyIt->HasKey("undermine"))
		{
			_bts.updateSurfaceDraw();
			json::Array pArr, uArr = (*_historyIt)["undermine"].ToArray();
			float hTx[2], uv[2];
			int tri;
			int material;
			Vec3f hVec;
			json::Object uObj, pObj;
			for (int n = (int)uArr.size(), i = 0; i < n; ++i) {
				uObj = uArr[i].ToObject();
				pObj = uObj["underminePoint"].ToObject();
				material = pObj["material"].ToInt();
				bool ic = true;  // compatibility with old history files
				if(pObj.HasKey("incisionConnect"))
					ic = pObj["incisionConnect"].ToBool();
				pArr = pObj["historyTexture"].ToArray();
				hTx[0] = pArr[0].ToFloat();
				hTx[1] = pArr[1].ToFloat();
				pArr = pObj["displacement"].ToArray();
				hVec[0] = pArr[0].ToFloat();
				hVec[1] = pArr[1].ToFloat();
				hVec[2] = pArr[2].ToFloat();
				if (!getHistoryAttachPoint(material, hTx, hVec, tri, uv, false)) {
					std::string msg = "History file undermine point location failure.";
					historyAttachFailure(msg);
					return;
				}
				_incisions.addUndermineTriangle(tri, 2, ic);
			}
			_gl3w->drawAll();
			glfwSwapBuffers(_ffg->FFwindow);
			std::this_thread::sleep_for(std::chrono::milliseconds(800));
			_incisions.undermineSkin();
			_undermineTriangles.clear();
			physicsDone = false;
			_ffg->physicsDrag = true;
			tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
				_bts.updateOldPhysicsLattice();
				newTopology = true;
				physicsDone = true;
				}
			);
			++_historyIt;
		}
		else if (_historyIt->HasKey("excise"))
		{
			json::Object exciseObj = (*_historyIt)["excise"].ToObject();
			json::Array pArr;
			float hTx[2], uv[2];
			int tri;
			int material;
			Vec3f hVec;
			material = exciseObj["material"].ToInt();
			pArr = exciseObj["historyTexture"].ToArray();
			hTx[0] = pArr[0];
			hTx[1] = pArr[1];
			pArr.Clear();
			pArr = exciseObj["displacement"].ToArray();
			hVec[0] = pArr[0];
			hVec[1] = pArr[1];
			hVec[2] = pArr[2];
			if (!getHistoryAttachPoint(material, hTx, hVec, tri, uv, false)) {
				std::string msg = "History file excise point location failure.";
				historyAttachFailure(msg);
				return;
			}
			_incisions.excise(tri);
			physicsDone = false;
			_ffg->physicsDrag = true;
			tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
				_bts.updateOldPhysicsLattice();
				newTopology = true;
				physicsDone = true;
				}
			);
			++_historyIt;
		}
		else if (_historyIt->HasKey("collisionProxySuture")) {
			json::Object sutureObj = (*_historyIt)["collisionProxySuture"].ToObject();
			int edge, sutNum = sutureObj["collisionSutureNum"].ToInt();
			float param, uv[2];
			if (_sutures.getNumberOfSutures() < 1) {	// initialize sutures
				_sutures.setSutureSize(_sg.getSceneNode()->getRadius()*0.003f);
				_sutures.setShapes(_gl3w->getShapes());
				_sutures.setGLmatrices(_gl3w->getGLmatrices());
				_sutures.setPhysicsLattice(_bts.getPdTetPhysics_2());
				_sutures.setVnBccTetrahedra(_bts.getVirtualNodedBccTetrahedra());
				_sutures.setDeepCut(&_incisions);
			}
			materialTriangles *tr = _sg.getMaterialTriangles();
			float hTx[2];
			Vec3f hVec(0.0f, 0.0f, 0.0f);
			json::Array pArr = sutureObj["mat4Texture"].ToArray();
			hTx[0] = pArr[0].ToFloat();
			hTx[1] = pArr[1].ToFloat();
			pArr.Clear();
			int eTri;
			if (!getHistoryAttachPoint(4, hTx, hVec, eTri, uv, false)) {
				std::string msg = "History file collision suture attachment failure.";
				historyAttachFailure(msg);
				return;
			}
			auto triEdge = [&](float(&uv)[2]) {
				if (uv[0] + uv[1] > 0.67f) {  // force to an edge
					edge = 1;
					param = uv[1] / (uv[0] + uv[1]);
				}
				else if (uv[0] > uv[1]) {
					edge = 0;
					param = uv[0];
				}
				else {
					edge = 2;
					param = 1.0f - uv[1];
				}
			};
			triEdge(uv);
			int k = _sutures.addUserSuture(tr, eTri, edge, param);
			_sutures.setLinked(k, false);
			pArr.Clear();
			pArr = sutureObj["mat5Texture"].ToArray();
			hTx[0] = pArr[0].ToFloat();
			hTx[1] = pArr[1].ToFloat();
			if (!getHistoryAttachPoint(5, hTx, hVec, eTri, uv, false)) {
				std::string msg = "History file collision suture attachment failure.";
				historyAttachFailure(msg);
				return;
			}
			triEdge(uv);
			_sutures.setSecondEdge(k, tr, eTri, edge, param);
			tr->getBarycentricPosition(eTri, uv, hVec._v);
			_sutures.setSecondVertexPosition(k, hVec._v);
			_bts.setForcesAppliedFlag();
			_hooks.selectHook(-1);	// deselect hooks
			_sutures.selectSuture(-1);
			++_historyIt;
		}
		else if (_historyIt->HasKey("addSuture"))
		{
			json::Object sutureObj = (*_historyIt)["addSuture"].ToObject();
			int edge, sutNum = sutureObj["sutureNum"].ToInt();
			float param, uv[2], xyz[3];
			if (_sutures.getNumberOfSutures() < 1) {	// initialize sutures
				_sutures.setSutureSize(_sg.getSceneNode()->getRadius()*0.003f);
				_sutures.setShapes(_gl3w->getShapes());
				_sutures.setGLmatrices(_gl3w->getGLmatrices());
				_sutures.setPhysicsLattice(_bts.getPdTetPhysics_2());
				_sutures.setVnBccTetrahedra(_bts.getVirtualNodedBccTetrahedra());
				_sutures.setDeepCut(&_incisions);
				_sutures.setSurgicalActions(this);
			}
			materialTriangles *tr = _sg.getMaterialTriangles();
			int material;
			float hTx[2];
			Vec3f hVec;
			json::Array pArr = sutureObj["historyTexture0"].ToArray();
			hTx[0] = pArr[0].ToFloat();
			hTx[1] = pArr[1].ToFloat();
			pArr.Clear();
			pArr = sutureObj["displacement0"].ToArray();
			hVec[0] = pArr[0].ToFloat();
			hVec[1] = pArr[1].ToFloat();
			hVec[2] = pArr[2].ToFloat();
			material = sutureObj["material0"].ToInt();
			int eTri;
			if (!getHistoryAttachPoint(material, hTx, hVec, eTri, uv, material == 2 ? true : false)) {
				std::string msg = "Attempted attachment of suture number ";
				msg.append(std::to_string(sutNum));
				msg.append(" in history file failed.");
				historyAttachFailure(msg);
				return;
			}
			assert(material == tr->triangleMaterial(eTri));
			if (material == 2) {
				if (uv[1] == 0.0f) {
					edge = 0;
					param = uv[0];
				}
				else if (uv[0] == 0.0f) {
					edge = 2;
					param = 1.0f - uv[1];
				}
				else if (uv[0] + uv[1] > 0.998f) {
					edge = 1;
					param = uv[1];
				}
				else
					assert(false);
			}
			else {
				if (uv[0] + uv[1] > 0.67f) {  // force to an edge
					edge = 1;
					param = uv[1] / (uv[0] + uv[1]);
				}
				else if (uv[0] > uv[1]) {
					edge = 0;
					param = uv[0];
				}
				else {
					edge = 2;
					param = 1.0f - uv[1];
				}
			}
			int sn = _sutures.addUserSuture(tr, eTri, edge, param);
			assert(_sutures.baseToUserSutureNumber(sn) == sutNum);
			pArr.Clear();
			pArr = sutureObj["historyTexture1"].ToArray();
			hTx[0] = pArr[0].ToFloat();
			hTx[1] = pArr[1].ToFloat();
			pArr.Clear();
			pArr = sutureObj["displacement1"].ToArray();
			hVec[0] = pArr[0].ToFloat();
			hVec[1] = pArr[1].ToFloat();
			hVec[2] = pArr[2].ToFloat();
			material = sutureObj["material1"].ToInt();
			if (!getHistoryAttachPoint(material, hTx, hVec, eTri, uv, material == 2 ? true : false)) {
				std::string msg = "Attempted attachment of suture number ";
				msg.append(std::to_string(sutNum));
				msg.append(" in history file failed.");
				historyAttachFailure(msg);
				return;
			}
			tr->getBarycentricPosition(eTri, uv, xyz);
			assert(material == tr->triangleMaterial(eTri));
			if (material == 2) {
				if (uv[1] == 0.0f) {
					edge = 0;
					param = uv[0];
				}
				else if (uv[0] == 0.0f) {
					edge = 2;
					param = 1.0f - uv[1];
				}
				else if (uv[0] + uv[1] > 0.998f) {
					edge = 1;
					param = uv[1];
				}
				else
					assert(false);
			}
			else {
				if (uv[0] + uv[1] > 0.67f) {  // force to an edge
					edge = 1;
					param = uv[1] / (uv[0] + uv[1]);
				}
				else if (uv[0] > uv[1]) {
					edge = 0;
					param = uv[0];
				}
				else {
					edge = 0;
					param = 1.0f - uv[1];
				}
			}
			int sRet = _sutures.setSecondEdge(sn, tr, eTri, edge, param);
			_bts.setForcesAppliedFlag();
			if (sRet < 1)
				_sutures.setSecondVertexPosition(sn, xyz);
			else if (sRet < 2) {
				sendUserMessage("Trying to suture to same side of incision is not allowed-", "USER ERROR", false);
				_selectedSurgObject = "";
				setToolState(0);
				++_historyIt;
				return;
			}
			else
				assert(false);
			if (sutureObj["linked"].ToBool()) {
				_sutures.setLinked(sn, true);
				physicsDone = false;
 				_ffg->physicsDrag = true;
//				tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {
					_sutures.laySutureLine(sn);
					physicsDone = true;
//					}
//				);
			}
			else
				_sutures.setLinked(sn, false);
			_hooks.selectHook(-1);	// deselect hooks
			_sutures.selectSuture(sn);
			char s[20];
#ifdef _WINDOWS
			sprintf_s(s, 19, "S_%d", sn);
#else
			sprintf(s, "S_%d", sn);
#endif
			_selectedSurgObject = s;
			++_historyIt;
			if (_historyIt == _historyArray.end()) {  // automatically promote any fake sutures if this is the last one
				while (!physicsDone)
					;
				_bts.promoteSutures();
			}
		}
		else if (_historyIt->HasKey("deleteSuture"))
		{
			int sutNum;
			if ((*_historyIt)["deleteSuture"].GetType() == json::ObjectVal) {
				sutNum = (*_historyIt)["deleteSuture"]["autoSuturesFor"].ToInt();
				sutNum = _sutures.userToBaseSutureNumber(sutNum);
				++sutNum;
				if (_sutures.baseToUserSutureNumber(sutNum) < 0)
					_sutures.deleteSuture(sutNum);
			}
			else {
				sutNum = (*_historyIt)["deleteSuture"].ToInt();
				sutNum = _sutures.userToBaseSutureNumber(sutNum);
				_sutures.deleteSuture(sutNum);
			}
			_sutures.selectSuture(-1);
			_hooks.selectHook(-1);
			++_historyIt;
		}
		else if (_historyIt->HasKey("makeDeepCut"))
		{
			// COURT - this is somewhat flawed. If physics state on creation is different than that on execution, the deep side of the normal may have very different outcomes.
			// consider assuring a certain number of physics iterations before execution.  Could also put position of deep post point in history and compute N.  Intermediate
			// incision intersections would still be different.
			json::Array pArr, iArr = (*_historyIt)["makeDeepCut"].ToArray();
			json::Object iObj = iArr[0].ToObject();
			assert(iObj["deepCutObject"].ToInt() == 0);	// for now only one object incisable
			bool startOpen, endOpen;
			startOpen = iObj["openIn"].ToBool();
			endOpen = iObj["openOut"].ToBool();
			int pointNum = iObj["pointNumber"].ToInt();
			materialTriangles* tr = _sg.getMaterialTriangles();
			_incisions.clearDeepCutter();
			if (!_fence.isInitialized()) {	// initialize fence
				_fence.setFenceSize(tr->getDiameter() * 0.01f);
				_fence.setGl3wGraphics(_gl3w);
			}
			_fence.clear();
			float hTx[2], uv[2];
			int tri;
			int material;
			Vec3f hVec, postN, xyz;
			json::Object uObj, pObj;
			for (int i = 0; i < pointNum; ++i) {
				uObj = iArr[i + 1].ToObject();
				if (!uObj.HasKey("deepCutPoint")) {
					sendUserMessage("There is an error in this history file.  Truncating from this point forward-", "", false);
					if (_historyIt != _historyArray.end()) {
						json::Array tarr;
						for (json::Array::ValueVector::iterator it = _historyArray.begin(); it != _historyIt; ++it)
							tarr.push_back(*it);
						_historyArray.Clear();
						_historyArray = tarr;
						_historyIt = _historyArray.end();
					}
					return;
				}
				pObj = uObj["deepCutPoint"].ToObject();
				material = pObj["material"].ToInt();
				pArr = pObj["historyTexture"].ToArray();
				hTx[0] = pArr[0].ToFloat();
				hTx[1] = pArr[1].ToFloat();
				pArr = pObj["displacement"].ToArray();
				hVec[0] = pArr[0].ToFloat();
				hVec[1] = pArr[1].ToFloat();
				hVec[2] = pArr[2].ToFloat();
				if (!getHistoryAttachPoint(material, hTx, hVec, tri, uv, false)) {
					std::string msg = "Can't retrieve deep cut point from history file.";
					historyAttachFailure(msg);
					return;
				}
				tr->getBarycentricPosition(tri, uv, xyz._v);
				pArr = pObj["postNormal"].ToArray();
				postN.X = pArr[0].ToFloat();
				postN.Y = pArr[1].ToFloat();
				postN.Z = pArr[2].ToFloat();
				if (i == pointNum - 1)
					startOpen = endOpen;
				else if (i > 0)
					startOpen = false;
				else
					;
				_fence.addPost(tr, tri, xyz._v, postN._v, false, true, startOpen);
			}
			if (!_incisions.inputCorrectFence(&_fence, _ffg)) {
				sendUserMessage("The deepCut in this history file failed.", "PROGRAM ERROR");
				_fence.clear();
				_incisions.clearDeepCutter();
				_ffg->setToolState(0);
				setToolState(0);
				_bts.setPhysicsPause(false);
				physicsDone = true;
				_ffg->physicsDrag = false;
				return;
			}
			_bts.updateSurfaceDraw();
			if (!_incisions.cutDeep()) {
				sendUserMessage("Attempted deepCut failed. Save history to debug.", "PROGRAM ERROR");
				_fence.clear();
				_incisions.clearDeepCutter();
				_ffg->setToolState(0);
				setToolState(0);
				_bts.setPhysicsPause(false);
				physicsDone = true;
				_ffg->physicsDrag = false;
				return;
			}
			physicsDone = false;
			_ffg->physicsDrag = true;
			tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
				_bts.updateOldPhysicsLattice();
				newTopology = true;
				physicsDone = true;
				}
			);
			_fence.clear();
			++_historyIt;
		}
		else if (_historyIt->HasKey("periostealUndermine")) {
			_bts.updateSurfaceDraw();
			json::Array pArr, uArr = (*_historyIt)["periostealUndermine"].ToArray();
			float hTx[2], uv[2];
			int tri;
			int material;
			Vec3f hVec;
			json::Object uObj, pObj;
			for (int n = (int)uArr.size(), i = 0; i < n; ++i) {
				uObj = uArr[i].ToObject();
				pObj = uObj["periostealTriangle"].ToObject();
				material = pObj["material"].ToInt();
				pArr = pObj["historyTexture"].ToArray();
				bool ic = true;
				if (pObj.HasKey("incisionConnect"))
					ic = pObj["incisionConnect"].ToBool();
				hTx[0] = pArr[0].ToFloat();
				hTx[1] = pArr[1].ToFloat();
				pArr = pObj["displacement"].ToArray();
				hVec[0] = pArr[0].ToFloat();
				hVec[1] = pArr[1].ToFloat();
				hVec[2] = pArr[2].ToFloat();
				if (!getHistoryAttachPoint(material, hTx, hVec, tri, uv, false)) {
					std::string msg = "Can't retrieve periosteal undermine point from history file.";
					historyAttachFailure(msg);
					return;
				}
				_incisions.addPeriostealUndermineTriangle(tri, hVec, ic);
			}
			// all periosteal undermine triangles now marked as material 10
			_incisions.clearCurrentUndermine(8);  // set all periosteal undermined triangles to material 8 and reset.
			_bts.fixPeriostealPeriferalVertices();
			physicsDone = false;
			_ffg->physicsDrag = true;
			tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
				_bts.nonTetPhysicsUpdate();
				newTopology = true;
				physicsDone = true;
				}
			);
			++_historyIt;
		}
		else if (_historyIt->HasKey("promoteSutureApproximations"))
		{
			_bts.promoteSutures();
			++_historyIt;
		}
		else if (_historyIt->HasKey("pausePhysics"))
		{
			_bts.setPhysicsPause(true);
			++_historyIt;
			return;  // don't setToolState(0) as will unpause physics
		}
		else
			++_historyIt;
		_ffg->setToolState(0);
		setToolState(0);
		_bts.setPhysicsPause(false);
	}  // end try block
	catch (std::exception &e) {
		sendUserMessage(e.what(), "Program error", false);
//		throw(std::logic_error("Program error in executing next history action"));
	}
}

bool surgicalActions::texturePickCode(const int triangle, const float (&uv)[2], float (&txUv)[2], float &triangleDuv, int &material)
{ // texture seam triangles will have a large deltaUV in cylindrical or spherical texture mapping
	// return true=user selected a top or bottom triangle, false if an edge triangle was selected
	float *tx[3],mm[4]={1e30f,-1e30f,1e30f,-1e30f},p=1.0f-uv[0]-uv[1];

	// COURT
	materialTriangles *tri = NULL;  //  _bts.getElasticSkinGraphics()->getMaterialTriangles();

	int *tr = tri->triangleTextures(triangle);
	if(tri->triangleMaterial(triangle)<0)
		return false;
	for(int j,i=0; i<3; ++i) {
		tx[i] = tri->getTexture(tr[i]);
		for(j=0; j<2; ++j) {
			if(mm[j<<1]>tx[i][j])
				mm[j<<1]=tx[i][j];
			if(mm[(j<<1)+1]<tx[i][j])
				mm[(j<<1)+1]=tx[i][j];
		}
	}
	for(int i=0; i<2; ++i)
		txUv[i] = p*tx[0][i] + uv[0] * tx[1][i] + uv[1] * tx[2][i];
	triangleDuv = mm[1]-mm[0] + mm[3]-mm[2];
	material = tri->triangleMaterial(triangle);
	return true;
}

bool surgicalActions::closestTexturePick(const float(&txUv)[2], const float triangleDuv, int &material, int &triangle, float(&uv)[2])
{ // at present search limited to top and bottom triangles. May change later.
	triangle = -1;

	// COURT
	materialTriangles *tri = NULL;  //  _bts.getElasticSkinGraphics()->getMaterialTriangles();

	std::vector<materialTriangles::matTriangle> *tArr = tri->getTriangleArray();
	int i,j,n=(int)tArr->size();
	int *tr;
	float *tx[3],minErr=1e30f,minDuv=1000.0f;
	float minimax[4];
	for(i=0; i<n; ++i)	{
		if(tri->triangleMaterial(i)!=material)
			continue;
		minimax[0] = 1e30f; minimax[1] = -1e30f; minimax[2] = 1e30f; minimax[3] = -1e30f;
		tr = tri->triangleTextures(i);  // COURT - if we switch to allowing edge picks this will be wrong. Would use lipGraphicsMt vertices instead of materialTriangles verts.
		for (j = 0; j<3; ++j) {
			tx[j] = tri->getTexture(tr[j]);
			if(minimax[0]>tx[j][0])
				minimax[0]=tx[j][0];
			if(minimax[1]<tx[j][0])
				minimax[1]=tx[j][0];
			if(minimax[2]>tx[j][1])
				minimax[2]=tx[j][1];
			if(minimax[3]<tx[j][1])
				minimax[3]=tx[j][1];
		}
		if (txUv[0] + 1e-5f<minimax[0] || txUv[0] - 1e-5f>minimax[1])
			continue;
		if (txUv[1] + 1e-5f<minimax[2] || txUv[1] - 1e-5f>minimax[3])
			continue;
		float err=0.0f,u,v,det=(tx[1][0]-tx[0][0])*(tx[2][1]-tx[0][1]) - (tx[1][1]-tx[0][1])*(tx[2][0]-tx[0][0]);
		if(fabs(det)<1e-16f)
			continue;
		u = (txUv[0] - tx[0][0])*(tx[2][1] - tx[0][1]) - (txUv[1] - tx[0][1])*(tx[2][0] - tx[0][0]);
		u /= det;
		v = (tx[1][0] - tx[0][0])*(txUv[1] - tx[0][1]) - (tx[1][1] - tx[0][1])*(txUv[0] - tx[0][0]);
		v /= det;
		if(u<-1e-4f)
			err += u*u;
		else if(u>1.0001f)
			err += (u-1.0f)*(u-1.0f);
		else ;
		if(v<-1e-4f)
			err += v*v;
		else if(v>1.0001f)
			err += (v-1.0f)*(v-1.0f);
		else ;
		if(err==0.0f) {
			if(u+v<1.0001f) {
				det = fabs(minimax[1]-minimax[0]+minimax[3]-minimax[2]-triangleDuv);
				if(det<minDuv) {
					minDuv = det;
					triangle = i;
					uv[0] = u;
					uv[1] = v;
				}
			}
			else {
				err = 1.0001f - u - v;
				err *= err;
			}
		}
		if(err<minErr && minDuv>100.0f) {
			minErr = err;
			triangle = i;
			uv[0] = u;
			uv[1] = v;
		}
	}
	if(minDuv>100.0f)
		return false;
	else
		return true;
}

bool surgicalActions::saveCurrentObj(const char* fullFilePath, const char* fileNamePrefix) {
	materialTriangles* tr = _sg.getMaterialTriangles();
	if (tr == nullptr)
		return false;
	return tr->writeObjFile(fullFilePath, fileNamePrefix);
}
