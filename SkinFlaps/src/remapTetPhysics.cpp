#include <assert.h>
#include <set>
#include "remapTetPhysics.h"

void remapTetPhysics::getOldPhysicsData(vnBccTetrahedra *oldVnbt)
{
	_oldNodePositions.clear();
	_oldNodePositions.reserve(oldVnbt->nodeNumber());
	for (int n = oldVnbt->nodeNumber(), i = 0; i < n; ++i)
		_oldNodePositions.push_back(oldVnbt->nodeSpatialCoordinate(i));
	_oldTets.clear();
	const auto otna = oldVnbt->getTetNodeArray();
	_oldTets.assign(otna.begin(), otna.end());  // just copy these as remakeNewTets() needs them.
	_oldCentroids.clear();
	_oldCentroids = std::move(oldVnbt->_tetCentroids);
	_oldTetHash.clear();
	_oldTetHash.reserve(_oldCentroids.size());
	_oldTetHash.insert(oldVnbt->_tetHash.begin(), oldVnbt->_tetHash.end());
	_oldVertexTets.clear();
	_oldVertexTets.assign(oldVnbt->_vertexTets.begin(), oldVnbt->_vertexTets.end());  // just copy these as remakeNewTets() needs them.
}

void remapTetPhysics::remapNewPhysicsNodePositions(vnBccTetrahedra *newVnbt)
{  // with new multires tet formulation all new physics nodes are no longer clones of old ones following decimation.
	// new low end nodes may be subnodes of larger tets.
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
