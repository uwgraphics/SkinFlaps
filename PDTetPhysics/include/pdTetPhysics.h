// Proposed interface for pdTetPhysics library:
//  Wang, Tao, Cutting, Sifakis
// Date: 1/30/2020

#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include "PDTetSolver.h"
#include "Utilities.h"


class pdTetPhysics {
	// a very unsafe wrapper class
private:
	using T = float;
	static constexpr int d = 3;
	PDTetSolver<T, d> m_solver;
	T m_hookWeight{ 0 };
	T m_sutureWeight{ 0 };
	T m_fixedWeight{ 0 };
	T m_peripheralWeight{ 0 };

	T m_stressLimit{ 1 };

	std::vector<int> fixedTetConstraints;

public:
	/* loaded with model file in history as static variables applied to all tet constraints. Later
	 * could have different properties for different tissues (e.g. cartilage versus skin
	 * versus breast) with different closed manifold surfaces around different types of tissue
	 * in the model.  If done that way global collision weight would be handled separately.
	 */

	inline void addCollisionObject(const std::string& collisionObjPath) {
		m_solver.addLevelSet(collisionObjPath);
	}

	void addSoftCollisionTets(const std::vector<int> &tets) {
		m_solver.addSelfCollisionElements(&tets[0], tets.size());
	}

	void addFixedCollisionSet(const std::string &levelSetFile, const std::vector<int> &tets, const std::vector<std::array<float, 3> > &weights) {
		if (!m_levelsetInited)
			addCollisionObject(levelSetFile);
		inputCollisionProxies(tets, weights);
			// called after every topological change as tets and weights will change.  Check to see if levelSetFile has already been loaded.  This only needs to be done on initial load and not repeated.
	}

	void currentSoftCollisionPairs(const std::vector<int> &topTets, const std::vector<std::array<float, 3> > &topBarys,
		const std::vector<int> &bottomTets, const std::vector<std::array<float, 3> > &bottomBarys, const std::vector<std::array<float, 3> > &collisionNormals) {
		assert(topTets.size() == topBarys.size() && topTets.size() == bottomTets.size() && topTets.size() == bottomBarys.size() && topTets.size() == collisionNormals.size());

		if (topTets.size())
			m_solver.updateCollisionSutures(topTets.size(), topTets.data(), bottomTets.data(), topBarys[0].data(), bottomBarys[0].data(), collisionNormals[0].data());
	}

	inline void setTetProperties(const float lowTetWeight, const float highTetWeight, const float strainMin, const float strainMax, const float collisionWeight, const float selfCollisionWeight, const float fixedWeight, const float peripheralWeight) {
		// guard against reset tetProperties
		if (m_tetPropsSet)
			throw std::logic_error("tet properties can only be set once");
		else
			m_tetPropsSet = true;
		m_fixedWeight = fixedWeight;
		m_peripheralWeight = peripheralWeight;
		m_solver.setParameters(highTetWeight / 2, lowTetWeight / highTetWeight, strainMin, strainMax, collisionWeight, selfCollisionWeight);
	}

	inline void tetSubset(const float lowTetWeight, const float highTetWeight, const float strainMin, const float strainMax, const std::vector<int>& tets) {
		m_solver.addSubset(highTetWeight / 2, lowTetWeight / highTetWeight, strainMin, strainMax, tets);

		// QISI - write me.  First 4 arguments are the tet properties of this subset of tets.  Last argument contains the indices of the tets these
		// properties should be assigned to.  This will be called after every topo change.  Currently there will be only one of these for cartilage.
		// In the there may be several of these for different tissue type within the body.
		std::cout << "calling tetSubset with set size : "<<tets.size() << std::endl;
	}

	inline bool solverInitialized() { return m_solverInited; }

	inline std::array<float, 3>* createBccTetStructure(const std::vector< std::array<int, 4> > &tetIndices, float tetScale) {
		m_solver.initializeDeformer(reinterpret_cast<const int(*)[4]>(&tetIndices[0][0]), tetIndices.size(), tetScale * 2);
		m_deformerInited = true;
		m_solverInited = false;
		fixedTetConstraints.clear();
		return reinterpret_cast<std::array<T, d>(*)>(m_solver.getPositionPtr());
	}

	inline std::array<float, 3>* createBccTetStructure_multires(const std::vector< std::array<int, 4> >& tetIndices, const std::vector<uint8_t>& tetSizeMultiples, float tetScale) {
		m_solver.initializeDeformer_multires(reinterpret_cast<const int(*)[4]>(&tetIndices[0][0]), reinterpret_cast<const uint8_t*>(&tetSizeMultiples[0]), tetIndices.size(), tetScale * 2);
		m_deformerInited = true;
		m_solverInited = false;
		fixedTetConstraints.clear();
		return reinterpret_cast<std::array<T, d>(*)>(m_solver.getPositionPtr());
	}

	// Next routine for inputting nodes on the face separating a large tet from possible multiple smaller ones. The subnodes input are present on a smaller tet, but not on the larger one.
	// These are constrained by internodeWeight to be barycentrically located on the larger face by faceNodes.  With multiple levels of physics resolution faceNodes and their
	// corresponding barycentric multipliers can number more than three for a single subtet.
	void addInterNodeConstraints(const std::vector<int>& subNodes, const std::vector<std::vector<int> >& faceNodes, const std::vector<std::vector<float> >& faceBarycentrics, const float internodeWeight) {
		int sns = subNodes.size();
		assert(sns == faceNodes.size() && sns == faceBarycentrics.size());
		for (int i = 0; i < sns; i++) {
			if (faceNodes[i].size() > 3) {
				// std::cout << i << "th internode needs " << faceNodes[i].size() << " macro nodes" << std::endl;
				m_solver.addInterNodeConstraint(subNodes[i], faceNodes[i].size(), &faceNodes[i][0], &faceBarycentrics[i][0], 0);  // QISI - currently not processing these?
			} else {
				int l = faceNodes[i].size() < 3 ? faceNodes[i].size() : 3;
				int fN[3]{}; for (int v = 0; v < l; ++v) fN[v] = faceNodes[i][v];
				float bC[3]{}; for (int v = 0; v < l; ++v) bC[v] = faceBarycentrics[i][v];
				int handle = m_solver.addInterNodeConstraint(subNodes[i], fN, bC, internodeWeight);
			}
		}
	}

	/* Doesn’t always follow createNewTetTopology() and may happen without changing
	 * tets. For example periosteal undermining releases some of these removing some
	 * Dirichlet constraints but all the tet constraints stay the same.	 */

	 // Court's new fixed vertex constraint version to replace setFixedNodes(). Should have done it this way in the first place. Sorry Qisi.
	// Also with a periosteal undermine there is no topo change, but previous fixed constraints are released so a pd reinit is required?
	inline void setFixedVertices(const std::vector<int> &fixedTets, const std::vector<std::array<float, 3> > &fixedWeights, const std::vector<std::array<float, 3> > &fixedPositions,
		const std::vector<int> &peripheralTets, const std::vector<std::array<float, 3> > &peripheralWeights, const std::vector<std::array<float, 3> > &peripheralPositions) {
		if (!m_deformerInited)
			throw std::logic_error("need to init tet topology before setFixedNodes");
		const T weight[d] = { 0,0,0 };

		// QISI - I see your m_grid_deformer handles tet constraints for hooks and fixed vertices.  Will hooks always be the last?  This is somewhat ugly.  Please advise.
		if (m_solver.numberOfTetConstraints() >= fixedTetConstraints.size()) {
			for (auto& tc : fixedTetConstraints)
				m_solver.deleteConstraint(tc);
		}
		fixedTetConstraints.clear();

		assert(fixedTets.size() == fixedWeights.size() && fixedWeights.size() == fixedPositions.size());
		assert(peripheralTets.size() == peripheralWeights.size() && peripheralWeights.size() == peripheralPositions.size());
		fixedTetConstraints.reserve(fixedTets.size() + peripheralTets.size());
		int fts = fixedTets.size();
		for (int i = 0; i < fts; i++) {
			int handle = m_solver.addConstraint(reinterpret_cast<const int(&)[4]>(m_solver.getTetIndices(fixedTets[i])), reinterpret_cast<const T(&)[3]>(fixedWeights[i]), reinterpret_cast<const T(&)[3]>(fixedPositions[i]), m_fixedWeight); // change weight
			fixedTetConstraints.push_back(handle);
		}
		for (int i = 0; i < peripheralTets.size(); i++) {
			int handle = m_solver.addConstraint(reinterpret_cast<const int(&)[4]>(m_solver.getTetIndices(peripheralTets[i])), reinterpret_cast<const T(&)[3]>(peripheralWeights[i]), reinterpret_cast<const T(&)[3]>(peripheralPositions[i]), m_peripheralWeight); // change weight
			fixedTetConstraints.push_back(handle);
		}
	}

	/* returns constraint index */
	inline int addHook(const int tet, const std::array<float, 3> &barycentricWeight, const std::array<float, 3> &hookPosition, bool strong = false) {
		if (!m_deformerInited)
			throw std::logic_error("need to init tet topology before addHook");
		int number;
		if(strong)
			number = m_solver.addConstraint(tet, reinterpret_cast<const T(&)[d]>(barycentricWeight), reinterpret_cast<const T(&)[d]>(hookPosition), m_hookWeight*200.0f, m_stressLimit*2000.0f);
		else
			number = m_solver.addConstraint(tet, reinterpret_cast<const T(&)[d]>(barycentricWeight), reinterpret_cast<const T(&)[d]>(hookPosition), m_hookWeight, m_stressLimit);
		initializePhysics();
		return number;
	}

	inline void moveHook(const int hookHandle, const std::array<float, 3> &newPosition) {
		if (!m_deformerInited)
			throw std::logic_error("need to init tet topology before moveHook");
		m_solver.moveConstraint(hookHandle, reinterpret_cast<const T(&)[d]>(newPosition));
	}

	/* Could also just nullify it. */
	inline void deleteHook(const int hookHandle) {
		if (!m_deformerInited)
			throw std::logic_error("need to init tet topology before deleteHook");
		m_solver.deleteConstraint(hookHandle);
	}

	/*Sets static variables for these parameters. */
	inline void setHookSutureWeights(const float hookWeight, const float sutureWeight, const float stressLimit = FLT_MAX) {
		m_hookWeight = hookWeight;
		m_sutureWeight = sutureWeight;
		m_stressLimit = stressLimit;
	}

	/* returns constraint index */
	inline int addSuture(const int(&tets)[2], const std::array<float, 3>(&barycentricWeights)[2]) {
		if (!m_deformerInited)
			throw std::logic_error("need to init tet topology before addSuture");
		// m_solverInited = false;
		return m_solver.addSuture(tets, reinterpret_cast<const T(&)[2][d]>(barycentricWeights[0]), sqrt(m_sutureWeight));
	}

	/* Perhaps there should be a single delete(or nullify?)Constraint(int Handle); call for both hooks and sutures.*/
	inline void deleteSuture(const int sutureHandle) {
		if (!m_deformerInited)
			throw std::logic_error("need to init tet topology before deleteSuture");
		m_solver.deleteSuture(sutureHandle);
	}

	// After constraints have changed computes ATA and does its LDLT() if needed
	inline void initializePhysics() {
		if (m_solverInited) {
			reInitializePhysics();
		}
		else {
			// the tetProperties need to be set before calling this
			if (!m_deformerInited)
				throw std::logic_error("need to init tet topology before init solver");
			// guard against init with no tet properties
			if (!m_tetPropsSet)
				throw std::logic_error("need to set tetProperties before initializePhysics");
			initializeCollisionObject(0.03f);
			promoteAllSutures();
			m_solver.initializeSolver();
			m_solverInited = true;
		}
	}
	// Adding and deleting hooks can have an initializePhysics() call at the end of them as they are always unique. This should not be done with sutures since often a whole line of sutures
	// will be added after which only a single initializePhysics() call is necessary.
	private:
	inline void reInitializePhysics() {
		if (!m_solverInited)
			throw std::logic_error("need to init solver before reinit");
		m_solver.reInitializeSolver();
	}
	public:

	inline void inputCollisionProxies(const std::vector<int> &tets, const std::vector<std::array<float, 3> > &weights) {
		if (!m_deformerInited)
			throw std::logic_error("need to init tet topology before add proxies");
		m_solver.addCollisionProxies(&tets[0], reinterpret_cast<const T(*)[d]>(&weights[0]), tets.size());
	}

	// do least squares solve and process collisions
	inline void solve() {
		if (!m_solverInited)
			throw std::logic_error("need to init solver before solve");
		m_solver.solve();
	}

	pdTetPhysics() : m_tetPropsSet(false), m_solverInited(false), m_deformerInited(false), m_levelsetInited(false) {}

	~pdTetPhysics() {
		m_solver.releaseSolver();
		m_solver.releaseDeformer();
	}

	inline void promoteAllSutures() { m_solver.premoteSutures(); m_solverInited = false;}

	inline void initializeCollisionObject(const T levelSetDx) { if (!m_levelsetInited) { m_solver.initializeLevelSet(levelSetDx); m_levelsetInited = true; } }

private:
	// static variables as properties

	// cannot reinit if initted
	bool m_deformerInited;
	bool m_solverInited; // only init if deformer inited, only numfact if inited

	bool m_tetPropsSet;
	bool m_levelsetInited;
};
