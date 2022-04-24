#include <assert.h>
#include <set>
#include "remapTetPhysics.h"

void remapTetPhysics::getOldPhysicsData(vnBccTetrahedra *oldVnbt)
{
	_oldNodePositions.clear();
	_oldNodePositions.reserve(oldVnbt->nodeNumber());
	for (int n = oldVnbt->nodeNumber(), i = 0; i < n; ++i)
		_oldNodePositions.push_back(oldVnbt->_nodeSpatialCoords[i]);
//	_oldNodePositions = std::move(oldVnbt->_nodeSpatialCoords);
	_oldFixedNodes.clear();
	_oldFixedNodes.assign(oldVnbt->_fixedNodes.begin(), oldVnbt->_fixedNodes.end());
	_oldTets.clear();
	_oldTets.assign(oldVnbt->_tetNodes.begin(), oldVnbt->_tetNodes.end());  // just copy these as remakeNewTets() needs them.
	_oldCentroids.clear();
	_oldCentroids = std::move(oldVnbt->_tetCentroids);
	_oldTetHash.clear();
	_oldTetHash = std::move(oldVnbt->_tetHash);
	_oldVertexTets.clear();
	_oldVertexTets.assign(oldVnbt->_vertexTets.begin(), oldVnbt->_vertexTets.end());  // just copy these as remakeNewTets() needs them.
}

void remapTetPhysics::remapNewPhysicsNodePositions(vnBccTetrahedra *newVnbt)
{  // all new physics nodes are clones of old ones.  After a pure topo change their spatial positions will remain identical as no new forces have yet been applied.
	std::vector<long> newToOldTets;
	newToOldTets.assign(newVnbt->tetNumber(), -1);
	std::vector<char> oldTetsUsed;
	oldTetsUsed.assign(_oldTets.size(), 0);
	for (int n = _oldVertexTets.size(), j, i = 0; i < n; ++i){  // new verts won't have an identical old correspondence
		// quickest first by using old vertex correspondence
		j = _oldVertexTets[i];
		if (j < 0)  		// handles deleted vertices after an excise
			continue;
		oldTetsUsed[j] = 1;
		if (newToOldTets[newVnbt->_vertexTets[i]] < 0)
			newToOldTets[newVnbt->_vertexTets[i]] = j;
		else{
			if (newToOldTets[newVnbt->_vertexTets[i]] != j)
				int junk = 0;
	//		assert(newToOldTets[newVnbt->_vertexTets[i]] == j);
		}
	}
	// next find singleton unused old tets and assign to new matching centroid tets
	for (int n = oldTetsUsed.size(), i = 0; i < n; ++i){
		if (oldTetsUsed[i])
			continue;
		auto npair = _oldTetHash.equal_range(_oldCentroids[i].ll);
		if (std::distance(npair.first, npair.second) != 1){
			continue;
		}
		oldTetsUsed[i] = 1;  // ? nuke as not used anymore
		npair = newVnbt->_tetHash.equal_range(_oldCentroids[i].ll);
		// deleted cube from an excise operation can leave npair.first == npair.second
		while (npair.first != npair.second){
			if (newToOldTets[npair.first->second] < 0)
				newToOldTets[npair.first->second] = i;
			else {
				if (newToOldTets[npair.first->second] != i)
					int junk = 0;
				//		assert(newToOldTets[newVnbt->_vertexTets[i]] == j);
			}
			++npair.first;
		}
	}
	// get correspondence between new and old physics nodes found so far
	_newToOldNodes.assign(newVnbt->nodeNumber(), -1);
	for (int n = newToOldTets.size(), j, i = 0; i < n; ++i){
		if (newToOldTets[i] < 0)
			continue;
		long *oldNodes = _oldTets[newToOldTets[i]].data(), *newNodes = newVnbt->_tetNodes[i].data();
		for(j=0; j<4; ++j){
			if (_newToOldNodes[newNodes[j]] < 0)
				_newToOldNodes[newNodes[j]] = oldNodes[j];
			else{  // above two tests pass, but this fails frequently
				int junk = 0;
				//				assert(oldNodes[j] == _newToOldNodes[newNodes[j]]);
			}
		}
	}
	// remaining original tets are intersected by the trianglulated surface, but don't contain a vertex
	for (int n = newToOldTets.size(), i = 0; i < n; ++i){
		if (newToOldTets[i] > -1)
			continue;
		long *newNodes = newVnbt->_tetNodes[i].data();
		auto npair = _oldTetHash.equal_range(newVnbt->_tetCentroids[i].ll);
		assert(npair.first != npair.second);
		long *oldN;
		if (std::distance(npair.first, npair.second) == 1){  // this singleton match is certain. Always replace.
			oldN = _oldTets[npair.first->second].data();  // this singleton old tet was used for the _vertexTet match in the beginning
			for (int j = 0; j < 4; ++j){
#ifdef _DEBUG
//				if (_newToOldNodes[newNodes[j]] > -1)
//					assert(_newToOldNodes[newNodes[j]] == oldN[j]);
#endif
				_newToOldNodes[newNodes[j]] = oldN[j];
			}
		}
		else{
			// find best matched old cube. Last resort only replace node correspondence if desperate.
			long matched, bestOldCube, nMatched = -1;
			while (npair.first != npair.second){
				oldN = _oldTets[npair.first->second].data();
				matched = 0;
				for (int j = 0; j < 4; ++j){
					if (oldN[j] == _newToOldNodes[newNodes[j]])
						++matched;
				}
				if (matched > nMatched){
					nMatched = matched;
					bestOldCube = npair.first->second;
				}
				++npair.first;
			}
			assert(nMatched > 0);
			oldN = _oldTets[bestOldCube].data();
			for (int j = 0; j < 4; ++j){
				if (_newToOldNodes[newNodes[j]] < 0)  // use only if no better one has been found
					_newToOldNodes[newNodes[j]] = oldN[j];
			}
		}
	}
	_oldTets.clear();
	_oldCentroids.clear();
	_oldTetHash.clear();
	_oldVertexTets.clear();
}

void remapTetPhysics::restoreOldNodePositions(vnBccTetrahedra *newVnbt)
{
	// can't set new physics node fixation state.  Only old and new node positions match up after incisions, no real node corespondence.
	// move new physics nodes to their last spatial position
	for (int n = _newToOldNodes.size(), i = 0; i < n; ++i){
		assert(_newToOldNodes[i] > -1);
		newVnbt->_nodeSpatialCoords[i].set(_oldNodePositions[_newToOldNodes[i]]);
	}
	_newToOldNodes.clear();
	_oldNodePositions.clear();
}

remapTetPhysics::remapTetPhysics()
{
}


remapTetPhysics::~remapTetPhysics()
{
}
