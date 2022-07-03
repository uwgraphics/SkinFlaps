#include <assert.h>
#include <limits>
#include <algorithm>
#include <cmath>
#include <functional>
#include "boundingBox.h"
#include "Mat2x2d.h"
#include "Mat3x3f.h"

#include <omp.h>
#include <chrono>  // for openMP timing
#include <ctime>  // nuke after openMP debug
#include <fstream>

#include "tbb/parallel_for_each.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"

#include "surgicalActions.h"  // for debug.  May nuke later

#include "vnBccTetCutterTbb.h"

class PlanePolygonMaker {
	const int pSet;
	materialTriangles *mt;
	const Vec3d *matCoords;
	vnBccTetCutterTbb::bccPlane* planes;
public:
	void operator() (const tbb::blocked_range<std::size_t>& r) const {
		for (std::size_t plane = r.begin(); plane != r.end(); ++plane) {
			vnBccTetCutterTbb::bccPlane *pp = &planes[plane];
			std::vector<double> planeDist;
			planeDist.assign(mt->numberOfVertices(), 0.0);
			int c2, c1, c0 = pSet >> 1;
			c1 = (c0 + 1) % 3;
			c2 = (c0 + 2) % 3;
			bool odd = (pSet & 1), somePos = false, someNeg = false;
			for (int n = mt->numberOfVertices(), i = 0; i < n; ++i) {
				if (odd)
					planeDist[i] = -matCoords[i]._v[c0] + matCoords[i]._v[c1] + pp->D;
				else
					planeDist[i] = -matCoords[i]._v[c0] - matCoords[i]._v[c1] + pp->D;
				std::signbit(planeDist[i]) ? somePos = true : someNeg = true;
			}
			if (!(somePos & someNeg))
				return;
			int *tr;
			auto getIntersect = [&](int edge, Vec2d &intrsct) {
				double denom = abs(planeDist[tr[edge]] - planeDist[tr[(edge + 1) % 3]]);
				if (denom < 1e-16f) {
					intrsct.X = matCoords[tr[edge]][c2];
					intrsct.Y = matCoords[tr[edge]][c1];
				}
				else {
					denom = 1.0 / denom;
					const double *dp = matCoords[tr[(edge + 1) % 3]]._v;
					Vec2d I1(dp[c2], dp[c1]);
					dp = matCoords[tr[edge]]._v;
					intrsct.set(dp[c2], dp[c1]);
					intrsct *= abs(planeDist[tr[(edge + 1) % 3]]) * denom;
					intrsct += I1 * abs(planeDist[tr[edge]]) * denom;
				}
			};
			std::vector<char> trisDone;
			trisDone.assign(mt->numberOfTriangles(), 0);
			// for the cleft lip entry data set max number of polygon vertices generated for any plane was 596
			for (int n = mt->numberOfTriangles(), j, i = 0; i < n; ++i) {
				if (trisDone[i] || mt->triangleMaterial(i) < 0)
					continue;
				tr = mt->triangleVertices(i);
				for (j = 0; j < 2; ++j) {
					if (std::signbit(planeDist[tr[j]]) != std::signbit(planeDist[tr[(j + 1) % 3]]))  // polygon start
						break;
				}
				if (j > 1)
					continue;
				// list polygon counterclockwise in the plane as the triangle edge polygons will be listed that way as well
				// counterclockwise going edge of triangle will always have signbit going from false to true
				// triSegment2 elements in polygon list triangle and its EXITING uv.  Entering uv is not listed in a ts2.
				if (j > 0) {
					if (!std::signbit(planeDist[tr[1]]))
						j = 2;
				}
				else {
					if (!std::signbit(planeDist[tr[0]])) {  // not edge 0. Will be remaining edge crossing
						if (std::signbit(planeDist[tr[2]]) != std::signbit(planeDist[tr[0]]))
							j = 2;
						else
							j = 1;
					}
				}
				// this polygon starts at triangle i edge j
				pp->polygons2.push_back(std::list<vnBccTetCutterTbb::triSegment2>());
				vnBccTetCutterTbb::triSegment2 ts;
				ts.triangle = i;
				do {
					getIntersect(j, ts.uv);
					pp->polygons2.back().push_back(ts);
					unsigned int adj = mt->triAdjs(ts.triangle)[j];
					if (adj == 3)
						return;
					ts.triangle = adj >> 2;
					tr = mt->triangleVertices(ts.triangle);
					trisDone[ts.triangle] = 1;
					j = adj & 3;
					for (int k = 1; k < 3; ++k) {
						if (std::signbit(planeDist[tr[(j + k) % 3]]) != std::signbit(planeDist[tr[(j + k + 1) % 3]])) {  // next edge
							j = (j + k) % 3;
							break;
						}
					}
				} while (ts.triangle != i);
			}
		}
	}
	PlanePolygonMaker(materialTriangles *matTri, const int planeSet, const Vec3d *vertexMatCoords, vnBccTetCutterTbb::bccPlane *planeVec) : mt(matTri), pSet(planeSet), matCoords(vertexMatCoords), planes(planeVec)
	{}
};

class IntersectPlanePolygon {
	const int pSet;
	materialTriangles *mt;
	vnBccTetrahedra *vbt;
	vnBccTetCutterTbb::bccPlane *planes;
	vnBccTetCutterTbb::LocusVec *nodes;
public:
	void operator() (const tbb::blocked_range<std::size_t>& r) const {
		for (std::size_t plane = r.begin(); plane != r.end(); ++plane) {
			auto uvToPlaneFace = [](const Vec2d &v) ->std::pair<short, short> {
				// plane face indexing has V component in the horizontal span minimum
				// The U component alternates with even permutation U+V has even U as upward pointing triangles, odd U down pointing.
				// Similarly odd permutation U+V has odd U as upward pointing triangles, even U down pointing.
				std::pair<short, short> sp;
				sp.first = (short)std::floor(v.X);
				sp.second = (short)std::floor(v.Y);
				if ((sp.first + sp.second) & 1) {
					if (v.Y - sp.second < sp.first + 1.0 - v.X)
						--sp.first;
				}
				else {
					if (v.Y - sp.second > v.X - sp.first)
						--sp.first;
				}
				return sp;
			};
			vnBccTetCutterTbb::bccPlane *pp = &planes[plane];
			std::vector<vnBccTetCutterTbb::PLANE_LINE> planeHorizLines, planeDiagonals[2];
			auto getFaceEdgePath = [&planeHorizLines, &planeDiagonals](vnBccTetrahedra *vBT, const Vec2d *startUv, std::pair<short, short> startFace, const Vec2d *endUv, std::pair<short, short> endFace, int triangle, int secondAxis,
				std::vector<std::pair<short, short> > &triFaces) {  // , std::vector<vnBccTetCutterTbb::PLANE_LINE> &horizLines, std::vector<vnBccTetCutterTbb::PLANE_LINE>(&diagonals)[2]
				triFaces.clear();
				triFaces.reserve(4);
				int previousEdge = -1;
				vnBccTetCutterTbb::planeLineCrossing plc;
				plc.solidRight = 0;
				plc.triangle = triangle;
				Vec2d N;
				N = *endUv - *startUv;  // due to counterclockwiseness around solid, its solid surface normal is y, -x
				double s, t;  // s is parameter along line from startUv to endUv. t is parameter along a candidate triangle edge.
				std::pair<short, short> faceNow = startFace;
				while (faceNow != endFace) {
					if (faceNow != startFace)
						triFaces.push_back(faceNow);
					bool oddTri = ((faceNow.first + faceNow.second) & 1);
					int i;
					for (i = 0; i < 3; ++i) {
						if (i == previousEdge)
							continue;
						if (i < 1) {  // horizontal axis crossing?
							if (faceNow.second == endFace.second)
								continue;
							int horizIndx = faceNow.second + (oddTri ? 1 : 0);
							s = (horizIndx - startUv->Y) / N[1];
							if (s < 0.0 || s > 1.0)
								continue;
							if (!(s > 0.0) && N.Y >= 0.0) {  // On horizLine but can't be a starting edge in this direction
								assert(previousEdge < 0);
								continue;
							}
							t = s * N[0] + startUv->X;
							if (t < (double)faceNow.first || t >= (double)faceNow.first + 2.0)
								continue;
							plc.solidRight = N[1] < 0.0 ? 1 : 0;
							planeHorizLines[horizIndx].insert(std::make_pair(t, plc));
							if (oddTri)
								++faceNow.second;
							else
								--faceNow.second;
							previousEdge = i;
							break;
						}
						else if (i < 2) {  // diag0 axis crossing - increasing direction in uv is [1, 1].
							if (N.X == N.Y) // crossing not possible
								continue;
							double dy = startUv->Y - faceNow.second + (oddTri ? 1.0 : 0.0), dx = startUv->X - faceNow.first;
							if (dy == dx && ((oddTri && -N.X + N.Y >= 0.0) || (!oddTri && N.X - N.Y >= 0.0))) {    // use same test for being exactly on diagonal line as uvToPlaneFace() or possible roundoff error
								assert(previousEdge < 0);
								continue;
							}
							s = dx - dy;
							s /= (N.Y - N.X);
							t = startUv->X - faceNow.first + N.X * s + (oddTri ? -1.0 : 0.0);
							if (s < 0.0 || s > 1.000001 || t < 0.0 || t > 1.0)  // allow a little destination roundoff error
								continue;
							int diagIndx = vBT->_gridSize[secondAxis] + ((faceNow.first - faceNow.second + 1) >> 1);
							plc.solidRight = (N[1] - N[0]) < 0.0 ? 1 : 0;
							planeDiagonals[0][diagIndx].insert(std::make_pair(t + (double)faceNow.second, plc));
							if (oddTri)
								++faceNow.first;
							else
								--faceNow.first;
							previousEdge = i;
							break;
						}
						else {	// diag1 axis crossing - increasing direction in uv is [-1, 1].
							if (N.X == -N.Y) // crossing not possible
								continue;
							double dy = startUv->Y - faceNow.second, dx = faceNow.first + (oddTri ? 1.0 : 2.0) - startUv->X;
							if (dx == dy && ((oddTri && N.X + N.Y >= 0.0) || (!oddTri && -N.X + -N.Y >= 0.0))) {    // use same test for being exactly on diagonal line as uvToPlaneFace() or possible roundoff error
								assert(previousEdge < 0);
								continue;
							}
							s = dx - dy;
							s /= N.Y + N.X;
							t = startUv->Y - faceNow.second + N.Y * s;
							if (s < 0.0 || s > 1.000001 || t < 0.0 || t > 1.0)  // allow a little destination roundoff error
								continue;
							int diagIndx = (faceNow.second + faceNow.first + 2) >> 1;
							plc.solidRight = (-N[1] - N[0]) < 0.0 ? 1 : 0;
							planeDiagonals[1][diagIndx].insert(std::make_pair(t + (double)faceNow.second, plc));
							if (oddTri)
								--faceNow.first;
							else
								++faceNow.first;
							previousEdge = i;
							break;
						}
					}
					assert(i < 3);
				}
			};

			auto makeConnectedComponentTetFaces = [&planeHorizLines, &planeDiagonals, &pp](int planeSet, vnBccTetrahedra *vBT, std::pair <vnBccTetCutterTbb::PFIT, vnBccTetCutterTbb::PFIT> &face)
			{
				std::pair<short, short> fac = face.first->first;
				bool oddFace = ((fac.first + fac.second) & 1);
				vnBccTetCutterTbb::PLANE_LINE *h, *d0, *d1;
				if (oddFace)
					h = &planeHorizLines[fac.second + 1];
				else
					h = &planeHorizLines[fac.second];
				d0 = &planeDiagonals[0][vBT->_gridSize[((planeSet >> 1) + 1) % 3] + ((fac.first - fac.second + 1) >> 1)];
				d1 = &planeDiagonals[1][(fac.second + fac.first + 2) >> 1];
				std::list<std::list<vnBccTetCutterTbb::triSegment2> > edgePolygons;
				// first get all open edgePolygons
				vnBccTetCutterTbb::triSegment2 ts;
				auto getEdgeSolids = [&](vnBccTetCutterTbb::PLANE_LINE *pl, float low, float high, bool fillX, bool reverse) {
					auto lb = pl->lower_bound(low);
					auto ub = pl->upper_bound(high);
					if (lb == ub) {
						if (lb != pl->end() && !lb->second.solidRight) {  // entire edge is interior
							std::list<vnBccTetCutterTbb::triSegment2> lt;
							if (fillX)
								ts.uv.X = low;
							else
								ts.uv.Y = low;
							ts.triangle = -1;
							lt.push_back(ts);
							if (fillX)
								ts.uv.X = high;
							else
								ts.uv.Y = high;
							lt.push_back(ts);
							if (reverse) {
								lt.reverse();
								edgePolygons.push_front(lt);
							}
							else
								edgePolygons.push_back(lt);
						}
						return;
					}
					if (!lb->second.solidRight) {
						std::list<vnBccTetCutterTbb::triSegment2> lt;
						if (fillX)
							ts.uv.X = low;
						else
							ts.uv.Y = low;
						ts.triangle = -2;
						lt.push_back(ts);
						if (fillX)
							ts.uv.X = lb->first;
						else
							ts.uv.Y = lb->first;
						ts.triangle = lb->second.triangle;
						lt.push_back(ts);
						if (reverse) {
							lt.reverse();
							edgePolygons.push_front(lt);
						}
						else
							edgePolygons.push_back(lt);
						++lb;
					}
					while (lb != ub) {
						assert(lb->second.solidRight);
						std::list<vnBccTetCutterTbb::triSegment2> lt;
						if (fillX)
							ts.uv.X = lb->first;
						else
							ts.uv.Y = lb->first;
						ts.triangle = lb->second.triangle;
						lt.push_back(ts);
						++lb;
						if (lb == ub) {
							if (fillX)
								ts.uv.X = high;
							else
								ts.uv.Y = high;
							ts.triangle = -2;
						}
						else {
							if (fillX)
								ts.uv.X = lb->first;
							else
								ts.uv.Y = lb->first;
							ts.triangle = lb->second.triangle;
						}
						lt.push_back(ts);
						if (reverse) {
							lt.reverse();
							edgePolygons.push_front(lt);
						}
						else
							edgePolygons.push_back(lt);
						if (lb != ub)
							++lb;
					}
				};
				if (oddFace) {
					ts.uv.X = -3.0f;
					getEdgeSolids(d1, (float)fac.second, (float)fac.second + 1.0f, false, true);
					ts.uv.X = -1.0f;
					getEdgeSolids(d0, (float)fac.second, (float)fac.second + 1.0f, false, false);
					ts.uv.Y = fac.second + 1.0f;
					getEdgeSolids(h, (float)fac.first, (float)fac.first + 2.0f, true, true);
					// compute diagonal x values
					for (auto &tp : edgePolygons) {
						for (auto &p : tp) {
							if (p.uv.X < -2.0f)
								p.uv.X = fac.first + 1.0f - (p.uv.Y - fac.second);
							else if (p.uv.X < -0.5f)
								p.uv.X = fac.first + 1.0f + (p.uv.Y - fac.second);
							else;
						}
					}
				}
				else {
					ts.uv.Y = (float)fac.second;
					getEdgeSolids(h, (float)fac.first, (float)fac.first + 2.0f, true, false);
					ts.uv.X = -1.0f;
					getEdgeSolids(d0, (float)fac.second, (float)fac.second + 1.0f, false, true);
					ts.uv.X = -3.0f;
					getEdgeSolids(d1, (float)fac.second, (float)fac.second + 1.0f, false, false);
					// compute diagonal x values
					for (auto &tp : edgePolygons) {
						for (auto &p : tp) {
							if (p.uv.X < -2.0f)
								p.uv.X = fac.first + 2.0f - (p.uv.Y - fac.second);
							else if (p.uv.X < -0.5f)
								p.uv.X = fac.first + (p.uv.Y - fac.second);
							else;
						}
					}
				}
				// remove any doubled corner points
				auto epit = edgePolygons.begin();
				while (epit != edgePolygons.end()) {
					if (epit->back().triangle < 0) {  // all corners points should be doubled
						auto pNext = epit;
						++pNext;
						if (pNext == edgePolygons.end())
							pNext = edgePolygons.begin();
						assert(pNext->front().triangle < 0 && epit->back().uv == pNext->front().uv);
						epit->pop_back();
						if (epit != pNext) {
							epit->splice(epit->end(), *pNext);
							edgePolygons.erase(pNext);
						}
						else
							break;
					}
					else
						++epit;
				}
				std::list<std::list<vnBccTetCutterTbb::triSegment2> > facePolygons;
				// now splice or add (if closed interior polygon) triangle-edge strings traversing this triangle+
				while (face.first != face.second) {
					if (face.first->second.back().uv.X < FLT_MAX) {  // a closed interior polygon
						facePolygons.push_back(std::list<vnBccTetCutterTbb::triSegment2>());
						facePolygons.back().insert(facePolygons.back().end(), face.first->second.begin(), face.first->second.end());
						face.first->second.clear();
						//						facePolygons.back().splice(facePolygons.back().end(), face.first->second);
					}
					else {
						auto lp = &face.first->second;
						assert(lp->back().uv[0] > 1e32f);  // unclosed chain with this uv not inside face
						epit = edgePolygons.begin();
						while (epit != edgePolygons.end()) {
							if (epit->front().triangle == lp->back().triangle) {
								std::list<vnBccTetCutterTbb::triSegment2> lTs2(lp->begin(), lp->end());
								lTs2.pop_back();
								epit->splice(epit->begin(), lTs2);
								face.first->second.clear();
								if (epit->front().triangle == epit->back().triangle) {  // polygon closed - done
									facePolygons.splice(facePolygons.end(), edgePolygons, epit);
									epit = edgePolygons.end();
								}
								break;
							}
							++epit;
						}
						if (epit == edgePolygons.end()) {
							assert(face.first->second.empty());
							++face.first;
							continue;
						}
						auto epit2 = edgePolygons.begin();
						while (epit2 != edgePolygons.end()) {
							if (epit2 == epit) {
								++epit2;
								continue;
							}
							if (epit2->back().triangle == epit->front().triangle) {
								epit->splice(epit->begin(), *epit2);
								edgePolygons.erase(epit2);
								break;
							}
							++epit2;
						}
					}
					assert(face.first->second.empty());
					++face.first;
				}
				if (!edgePolygons.empty()) {  // single polygon of only interior corner points surrounding an interior polygon(s)
					assert(!facePolygons.empty());
					assert(edgePolygons.size() < 2);
					facePolygons.splice(facePolygons.end(), edgePolygons);
				}
				// All polygons now are closed polygons with inside edge vertices labelled as < 0
				// Next section was createSurfaceTetFaces(fac, facePolygons, pp);
				std::list<std::list<vnBccTetCutterTbb::triSegment2> > exteriorPolygons, interiorPolygons;
				std::list<vnBccTetCutterTbb::vnTetFace> exteriorTetFaces;
				auto counterClockwise = [](std::list<vnBccTetCutterTbb::triSegment2> &poly) ->bool {
					std::list<vnBccTetCutterTbb::triSegment2>::iterator minIt, pit;
					double minX=DBL_MAX;
					for (pit = poly.begin(); pit != poly.end(); ++pit) {
						if (pit->uv[0] < minX) {
							minX = pit->uv[0];
							minIt = pit;
						}
					}
					Vec2d uv0, uv1;
					pit = minIt;
					do {
						if (pit == poly.begin())
							pit = poly.end();
						--pit;
					} while (pit->uv == minIt->uv && pit != minIt);
					if (pit == minIt)  // zero area polygon
						return true;
					uv0 = pit->uv - minIt->uv;
					pit = minIt;
					do {
						++pit;
						if (pit == poly.end())
							pit = poly.begin();
						uv1 = pit->uv - minIt->uv;
						double a = uv0[0] * uv1[1] - uv0[1] * uv1[0];
						if (a > 0.0)
							return false;
						else if (a < 0.0)
							return true;
						else;
					} while (pit != minIt);
					return true;  // zero area polygon
				};
				auto windingNumber = [](Vec2d &uv, std::list<vnBccTetCutterTbb::triSegment2> &poly) ->int {
					// Winding number test for a point in a polygon.  Assumes closed polygon without repeat of first point.
					auto lit = poly.end();
					--lit;
					int wn = 0;
					for (auto pit = poly.begin(); pit != poly.end(); ++pit) {
						if (lit->uv.Y <= uv.Y) {  // start
							if (pit->uv.Y > uv.Y)  // an upward crossing
								if (((pit->uv.X - lit->uv.X) * (uv.Y - lit->uv.Y) - (uv.X - lit->uv.X) * (pit->uv.Y - lit->uv.Y)) > 0.0)  // uv left of  edge
									++wn;            // have  a valid up intersect
						}
						else {  // lit->uv.Y > uv.y
							if (pit->uv.Y <= uv.Y)     // a downward crossing
								if (((pit->uv.X - lit->uv.X) * (uv.Y - lit->uv.Y) - (uv.X - lit->uv.X) * (pit->uv.Y - lit->uv.Y)) < 0.0)  // uv right of  edge
									--wn;            // have  a valid down intersect
						}
						lit = pit;
					}
					return wn;  // winding number.  0 only when P is outside. Sign gives clockwiseness.
				};
				while (!facePolygons.empty()) {
					auto fpit = facePolygons.begin();
					vnBccTetCutterTbb::vnTetFace tf;
					tf.interiorNodes = 0;
					auto lastIt = fpit->begin();
					bool exterior = false;
					for (auto tsit = fpit->end(); tsit != fpit->begin();) {
						--tsit;
						if (tsit->triangle < 0) {  // interior corner point
							exterior = true;
							// COURT - decided to have plane horizontal be vertex0 to vertex1 for both odd and even faces. This means even is counterclockwise and odd clockwise
							if (tsit->uv[0] == fac.first) {
								tf.interiorNodes |= 1;
								assert(tsit->uv[1] == fac.second + ((fac.first + fac.second) & 1) ? 1 : 0);
							}
							else if (tsit->uv[0] == fac.first + 2.0) {
								tf.interiorNodes |= 2;
								assert(tsit->uv[1] == fac.second + ((fac.first + fac.second) & 1) ? 1 : 0);
							}
							else if (tsit->uv[0] == fac.first + 1.0) {  // doesn't matter if up or down
								tf.interiorNodes |= 4;
								assert(tsit->uv[1] == fac.second + ((fac.first + fac.second) & 1) ? 0 : 1);
							}
							else
								assert(false);
							if (tsit != fpit->begin()) {
								--tsit;
								if (tsit->triangle > -1)
									tf.edgeTriangles.push_back(tsit->triangle);
								else
									++tsit;
							}
						}
						else {
							if (tsit->triangle == lastIt->triangle) {  // contains a border segment so is an exterior solid polygon
								tf.edgeTriangles.push_back(tsit->triangle);
								--tsit;
								if (tsit->triangle < 0) {
									++tsit;
									continue;
								}
								assert(tsit->triangle != tf.edgeTriangles.back());
								tf.edgeTriangles.push_back(tsit->triangle);
								exterior = true;
							}
							else {
								if (tf.edgeTriangles.empty() || tsit->triangle != tf.edgeTriangles.front())
									tf.interiorTriangles.push_back(tsit->triangle);
							}

						}
						lastIt = tsit;
					}
					if (exterior) {
						exteriorTetFaces.push_back(std::move(tf));
						assert(counterClockwise(*fpit));
						exteriorPolygons.splice(exteriorPolygons.end(), facePolygons, fpit);
					}
					else {  // test this interior polygon for clockwiseness
						if (counterClockwise(*fpit)) {
							exteriorTetFaces.push_back(std::move(tf));
							exteriorPolygons.splice(exteriorPolygons.end(), facePolygons, fpit);
						}
						else
							interiorPolygons.splice(interiorPolygons.end(), facePolygons, fpit);
					}
				}
				if (!interiorPolygons.empty()) {  // any interior must be inside an exterior one, so add its triangles to its interior tf
					for (auto &ip : interiorPolygons) {
						auto efit = exteriorTetFaces.begin();
						auto epit = exteriorPolygons.begin();
						while (epit != exteriorPolygons.end()) {
							if (windingNumber(ip.front().uv, *epit) != 0) {
								for (auto &ts2 : ip)
									efit->interiorTriangles.push_back(ts2.triangle);
								break;
							}
							++epit;
							++efit;
						}
						assert(epit != exteriorPolygons.end());
					}
				}
				// create a new tet face for each exteriorPolygon
				for (auto &etf : exteriorTetFaces)
					pp->vnTetFaces.insert(std::make_pair(fac, std::move(etf)));
			};

			int first = pSet >> 1;
			int second = (first + 1) % 3;
			int third = (first + 2) % 3;
			int diagDim = vbt->_gridSize[second] + vbt->_gridSize[third] + 1;
			planeDiagonals[0].resize(diagDim);
			planeDiagonals[1].resize(diagDim);
			int phlDim = (vbt->_gridSize[second] << 1) + 1;
			planeHorizLines.resize(phlDim);
			pp->vnTetFaces.clear();
			std::list<std::list<vnBccTetCutterTbb::triSegment2> > *L = &(pp->polygons2);
			if (L->empty())
				return;
			std::multimap<std::pair<short, short>, std::vector<vnBccTetCutterTbb::triSegment2> > planeFaces;
			for (auto lit = L->begin(); lit != L->end(); ++lit) {
				Vec2d lastUv = lit->back().uv;
				std::pair<short, short> lp, lpLast = uvToPlaneFace(lastUv);
				std::vector<vnBccTetCutterTbb::triSegment2> segmentLast, segmentNow;
				segmentLast.reserve(lit->size());
				segmentNow.reserve(lit->size());
				std::vector<std::pair<short, short> > triFaces;
				bool segmentLastDone = false;
				for (auto pit = lit->begin(); pit != lit->end(); ++pit) {  // pit contains triangle and its exiting uv.  Entering uv not listed.
					if (segmentLastDone)
						segmentNow.push_back(*pit);
					else
						segmentLast.push_back(*pit);
					lp = uvToPlaneFace(pit->uv);
					if (lp != lpLast) {
						if (segmentLastDone) {
							segmentNow.back().uv[0] = DBL_MAX;  // last uv not inside previous face.  Mark as an unclosed polygon segment
							auto pfit = planeFaces.insert(std::make_pair(lpLast, segmentNow));
							pfit->second.shrink_to_fit();
						}
						else {
							segmentLast.back().uv[0] = DBL_MAX;  // same logic as above
							segmentLastDone = true;
						}
						// next routine collects face path of this polygon segment and logs any triangle edge intersections
						getFaceEdgePath(vbt, &lastUv, lpLast, &pit->uv, lp, pit->triangle, second, triFaces);  // , planeHorizLines, planeDiagonals
						segmentNow.clear();
						// this uv is in the new face
						segmentNow.push_back(*pit);
						if (!triFaces.empty()) {  // complete triangle pass through these faces with no uv inside
							double firstU = pit->uv[0];
							segmentNow.back().uv[0] = DBL_MAX;
							for (auto &tf : triFaces) {
								auto pfit = planeFaces.insert(std::make_pair(tf, segmentNow));
								pfit->second.shrink_to_fit();
							}
							segmentNow.back().uv[0] = firstU;
						}
						lpLast = lp;
					}
					lastUv = pit->uv;
				}
				assert(segmentLast.size() + segmentNow.size() <= lit->size());
				if (segmentLastDone) {
					segmentNow.insert(segmentNow.end(), segmentLast.begin(), segmentLast.end());
					auto pfit = planeFaces.insert(std::make_pair(lpLast, segmentNow));
					pfit->second.shrink_to_fit();
				}
				else {  // closed polygon inside one triangular face
					auto pfit = planeFaces.insert(std::make_pair(lpLast, segmentLast));
					pfit->second.shrink_to_fit();
				}
			}
			L->clear();
			vnBccTetCutterTbb::PLANE_LINE::iterator hit;
			for (int k, j = 0; j < diagDim; ++j) {
				if (!planeDiagonals[0][j].empty()) {
					k = 1;
					for (hit = planeDiagonals[0][j].begin(); hit != planeDiagonals[0][j].end(); ++hit) {
						if (hit->second.solidRight != k) {  // set from polygon in getFaceEdgePath()
							// may be from floating point roundoff of coincident surface
							vnBccTetCutterTbb::planeLineCrossing plc = hit->second;
							auto hit2 = hit;
							++hit2;
							assert(hit2 != planeDiagonals[0][j].end() && hit2->first - hit->first < 1e-6);
							hit->second = hit2->second;
							hit2->second = plc;
						}
						k ^= 1;
					}
					assert(k > 0);
				}
				if (planeDiagonals[1][j].empty())	continue;
				k = 1;
				for (hit = planeDiagonals[1][j].begin(); hit != planeDiagonals[1][j].end(); ++hit) {
					if (hit->second.solidRight != k) {  // set from polygon in getFaceEdgePath()
						// may be from floating point roundoff of coincident surface
						vnBccTetCutterTbb::planeLineCrossing plc = hit->second;
						auto hit2 = hit;
						++hit2;
						assert(hit2 != planeDiagonals[1][j].end() && hit2->first - hit->first < 1e-6);
						hit->second = hit2->second;
						hit2->second = plc;
					}
					k ^= 1;
				}
				assert(k > 0);
			}
			for (int k, j = 0; j < phlDim; ++j) {
				if (planeHorizLines[j].empty())	continue;
				k = 1;
				for (hit = planeHorizLines[j].begin(); hit != planeHorizLines[j].end(); ++hit) {
					if (hit->second.solidRight != k) {  // set from polygon in getFaceEdgePath()
						// may be from floating point roundoff of coincident surface
						vnBccTetCutterTbb::planeLineCrossing plc = hit->second;
						auto hit2 = hit;
						++hit2;
						assert(hit2 != planeHorizLines[j].end() && hit2->first - hit->first < 3e-6);
						hit->second = hit2->second;
						hit2->second = plc;
					}
					k ^= 1;
				}
				assert(k > 0);
			}
			if (pSet < 1) {	// save interior tet nodes.  Careful here with openMP. Would need to create mutex for interiorNodeLatticeLoci. ?not worth the trouble.
				std::vector< std::array<short, 3> > nVec;
				assert(planeHorizLines[0].empty());
				assert(planeHorizLines[phlDim - 1].empty());
				for (int j = 1; j < phlDim - 1; ++j) {
					if (planeHorizLines[j].empty())	continue;
					std::array<short, 3> ll;
					ll[0] = pp->D - j;
					ll[1] = j;
					for (auto hit = planeHorizLines[j].begin(); hit != planeHorizLines[j].end(); ++hit) {
						short z = (short)std::floor(hit->first);
						assert(hit->second.solidRight);
						++z;	++hit;
						if (hit->second.solidRight) {  // this is floating point indecision at collision of coincident surfaces
							assert(false);
							// write me
						}
						// COURT - currently making them, putting into hash table later.
						while (z < hit->first) {
							if ((j & 1) == (z & 1)) {  // all nodes in lattice have either all odd or all even coordinates
								ll[2] = z;
								nVec.push_back(ll);
							}
							++z;
						}
					}
				}
				nodes->grow_by(nVec.begin(), nVec.end());
			}
			auto fit = planeFaces.begin();
			while (fit != planeFaces.end()) {
				auto pr = planeFaces.equal_range(fit->first);
				makeConnectedComponentTetFaces(pSet, vbt, pr);
				fit = pr.second;
			}
		}
	}

	IntersectPlanePolygon(materialTriangles *matTri, vnBccTetrahedra *vBT, const int planeSet, vnBccTetCutterTbb::bccPlane *planeVec, vnBccTetCutterTbb::LocusVec *nodeLoci)
		: mt(matTri), vbt(vBT), pSet(planeSet), planes(planeVec), nodes(nodeLoci)
	{}
};

vnBccTetCutterTbb::vnBccTetCutterTbb(void)
{
}

vnBccTetCutterTbb::~vnBccTetCutterTbb(void)
{
}

bool vnBccTetCutterTbb::makeFirstVnTets(materialTriangles *mt, vnBccTetrahedra *vbt, int maximumGridDimension)
{  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	// WARNING - no complete tests are done to check for non-self-intersecting closed manifold triangulated surface input!!
	// This is essential. findAdjacentTriangles() is the closest test I provide.
	_mt = mt;
	_vbt = vbt;
	_vbt->_mt = mt;
	_vbt->_vertexTets.clear();
	_vbt->_vertexTets.assign(mt->numberOfVertices(), -1);
	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(mt->numberOfVertices(), Vec3f());
	_surfaceTetFaceNumber = 0;
	_vbt->_fixedNodes.clear();
	_vbt->_tetHash.clear();
	_vbt->_tetNodes.clear();
	_vbt->_nodeSpatialCoords = nullptr;
	if (_mt->findAdjacentTriangles(true, false))	return false;
	if (!setupBccIntersectionStructures(maximumGridDimension))
		return false;
	for (int i = 0; i < 6; ++i) {
		int nS = _planeSetsTbb[i].size();
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nS), PlanePolygonMaker(_mt, i, &_vMatCoords[0], &_planeSetsTbb[i][0]));  // best to leave default grain size
	}
	// first plane set creates all the interior tetrahedral nodes.
	_nodeLoci.clear();
	_nodeLoci.reserve((_planeSetsTbb[0].size()*_vbt->_gridSize[1] * _vbt->_gridSize[2]) >> 2);  // reasonable guess
	for (int i = 0; i < 6; ++i) {
		int nS = _planeSetsTbb[i].size();
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nS), IntersectPlanePolygon(_mt, _vbt, i, &_planeSetsTbb[i][0], &_nodeLoci));
	}
	_interiorNodes.clear();  // only nodes created thus far are interior nodes
	_interiorNodes.reserve(_nodeLoci.size());
	for (int n = _nodeLoci.size(), j = 0; j < n; ++j)
		_interiorNodes.insert(std::make_pair(_nodeLoci[j], j));
	collectSurfaceTetCentersFaces();
	createVirtualNodedSurfaceTets();
	createSurfaceTetNodes();
	fillNonVnTetCenter();
	_vbt->_tetCentroids.clear();
	_vbt->_tetCentroids.reserve(_tets.size());
	_vbt->_tetNodes.clear();
	_vbt->_tetNodes.reserve(_tets.size());
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_tets.size());
	for (int n = _tets.size(), i = 0; i < n; ++i) {  // since vnBccTetrahedra not tbb aware can't multithread, but this is quite fast.
		_vbt->_tetNodes.push_back(_tets[i].tetNodes);
		_vbt->_tetCentroids.push_back(_tets[i].centroid);
		_vbt->_tetHash.insert(std::make_pair(_tets[i].centroid.ll, i));
	}
	_vbt->_nodeGridLoci.clear();
	_vbt->_nodeGridLoci.assign(_nodeLoci.begin(), _nodeLoci.end());
	_nodeLoci.clear();

	// all vertices in duplicated, virtual noded tets have their tet index already assigned.  Vertices in unique tets still unassigned.
	for (int n = _vbt->_vertexTets.size(), i = 0; i < n; ++i){
		if (_vbt->_vertexTets[i] < 0) {
			if (_vertexTetLoci[i].ll < LLONG_MAX) {  // checks for excisions
				auto vf = _vbt->_tetHash.find(_vertexTetLoci[i].ll);
				_vbt->_vertexTets[i] = vf->second;
			}
		}
	}
	_vertexTetLoci.clear();
	return true;
}

void vnBccTetCutterTbb::createSurfaceTetNodes()
{
	struct sharedTetNode {
		bool internal;
		std::set<unsigned int> nSet;
	}stn;
	stn.internal = false;
	stn.nSet.clear();
	typedef tbb::concurrent_unordered_map<std::array<short, 3>, std::list<sharedTetNode>, arrayShort3Hasher> NODESETS;
	NODESETS surfaceTetNodes((size_t)(_tets.size() * 0.37f));  // generous guess to avoid rehashing
	auto addTetNodePair = [&](int triVert, std::array<short, 3> &nLoc, vnTetFace *tfp) {
		unsigned int v0 = tfp->tetNodes[0][triVert], v1 = tfp->tetNodes[1][triVert];
		assert(v0 < 0xffffffff && v1 < 0xffffffff);  // Should never happen.
		bool internalNode;
		if (tfp->interiorNodes & (1 << triVert))
			internalNode = true;
		else
			internalNode = false;
		auto tnSets = &surfaceTetNodes.insert(std::make_pair(nLoc, std::list<sharedTetNode>())).first->second;
		auto tsit = tnSets->begin();
		while (tsit != tnSets->end()) {
			if (tsit->nSet.find(v0) != tsit->nSet.end()) {
				tsit->internal |= internalNode;  // if any face says it is an internal node, they all are
				if (tsit->nSet.insert(v1).second) {
					auto tsit2 = tnSets->begin();
					while (tsit2 != tnSets->end()) {
						if (tsit == tsit2) {
							++tsit2;
							continue;
						}
						assert(tsit2->nSet.find(v0) == tsit2->nSet.end());
						if (tsit2->nSet.find(v1) != tsit2->nSet.end()) {
							tsit->internal |= tsit2->internal;
							tsit->nSet.insert(tsit2->nSet.begin(), tsit2->nSet.end());
							tnSets->erase(tsit2);
							tsit2 = tnSets->end();
							break;
						}
						++tsit2;
					}
				}
				break;
			}
			if (tsit->nSet.find(v1) != tsit->nSet.end()) {
				tsit->internal |= internalNode;  // if any face says it is an internal node, they all are
				tsit->nSet.insert(v0);  // must happen
				auto tsit2 = tnSets->begin();
				while (tsit2 != tnSets->end()) {
					if (tsit == tsit2) {
						++tsit2;
						continue;
					}
					assert(tsit2->nSet.find(v1) == tsit2->nSet.end());
					if (tsit2->nSet.find(v0) != tsit2->nSet.end()) {
						tsit->internal |= tsit2->internal;
						tsit->nSet.insert(tsit2->nSet.begin(), tsit2->nSet.end());
						tnSets->erase(tsit2);
						tsit2 = tnSets->end();
						break;
					}
					++tsit2;
				}
				break;
			}
			++tsit;
		}
		if (tsit == tnSets->end()) {
			tnSets->push_back(stn);
			tnSets->back().internal = internalNode;
			tnSets->back().nSet.insert(v0);
			tnSets->back().nSet.insert(v1);
		}
	};
	for (int i = 0; i < 6; ++i) {
		int nS = _planeSetsTbb[i].size();
		tbb::parallel_for(tbb::blocked_range<size_t>(0, nS), [&](const tbb::blocked_range<size_t>& r) {
			for (size_t j = r.begin(); j != r.end(); ++j) {
				int c0 = i >> 1, oddPlane = i & 1;
				std::array<short, 3> nodeLoc;
				for (auto &tf : _planeSetsTbb[i][j].vnTetFaces) {
					unsigned short faceOdd = (tf.first.first + tf.first.second) & 1;
					for (int k = 0; k < 3; ++k) {
						nodeLoc[c0] = oddPlane ? _planeSetsTbb[i][j].D + tf.first.second + faceOdd : _planeSetsTbb[i][j].D - tf.first.second - faceOdd;
						nodeLoc[(c0 + 1) % 3] = tf.first.second + faceOdd;
						nodeLoc[(c0 + 2) % 3] = tf.first.first;
						if (k == 1)
							nodeLoc[(c0 + 2) % 3] += 2;
						if (k == 2) {
							++nodeLoc[(c0 + 2) % 3];
							int vOffset = faceOdd ? -1 : 1;
							nodeLoc[(c0 + 1) % 3] += vOffset;
							nodeLoc[c0] += oddPlane ? vOffset : -vOffset;
						}
						addTetNodePair(k, nodeLoc, &tf.second);
					}
				}
				_planeSetsTbb[i][j].vnTetFaces.clear();
			}
		});
	}
	// create and assign all shared surface tet nodes just found
	// serial version
	tbb::parallel_for_each(surfaceTetNodes.begin(), surfaceTetNodes.end(), [&](NODESETS::value_type &stn) {
		for (auto &sharedNode : stn.second) {
			int node;
			if (sharedNode.internal) {
				auto in = _interiorNodes.find(stn.first);
				assert(in != _interiorNodes.end());
				node = in->second;
			}
			else {
				auto nit = _nodeLoci.push_back(stn.first);
				node = (int)(nit - _nodeLoci.begin());
			}
			for (auto &tetNode : sharedNode.nSet) {
				_tets[tetNode >> 2].tetNodes[tetNode & 3] = node;
			}
		}
	});
	// In some complex solids createVirtualNodedSurfaceTets() will create multiple tets for the same tet locus as they will have no shared connected components,
	// but they can have the same nodes since the tets surrounding them may have connected components.  While these are not strictly an error, they do the user no good
	// since the tets are incapable of independent movement.  Further these are tets that are sparsely filled with solid.  Combining them provides somewhat more accurate
	// physics behavior. The next section removes these identical multiplicities.  Can't multithread this section.
	std::vector<int> tetMap;
	tetMap.assign(_tets.size(), -1);
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_tets.size());  // Used to remove surface tet duplicates and prevent dups in making center tets.
	for (int n = _tets.size(), i = 0; i < n; ++i) {
		for (int j = 0; j < 4; ++j) {
			int *tn = &_tets[i].tetNodes[j];
			// single interior polygon face tet can have 3 faces unpenetrated leaving one vertex isolated.
			if (*tn < 0) {
				std::array<short, 3> locus = _vbt->nodeGridLocation(_tets[i].centroid, j);
				auto nit = _nodeLoci.push_back(locus);
				*tn = (int)(nit - _nodeLoci.begin());
			}
		}
		auto pr = _vbt->_tetHash.equal_range(_tets[i].centroid.ll);
		while (pr.first != pr.second) {
			if (_tets[pr.first->second].tetNodes == _tets[i].tetNodes) {
				tetMap[i] = pr.first->second;
				break;
			}
			++pr.first;
		}
		if (pr.first == pr.second)
			_vbt->_tetHash.insert(std::make_pair(_tets[i].centroid.ll, i));
	}
	// remove surface tet duplicates
	size_t m = 0;
	for (int n = _tets.size(), i = 0; i < n; ++i) {
		if (tetMap[i] > -1) {
			assert(tetMap[i] < i && tetMap[i] > -1);
			tetMap[i] = tetMap[tetMap[i]];
		}
		else {
			if (m < i)
				_tets[m] = _tets[i];
			tetMap[i] = m++;
		}
	}
	_tets.resize(m);
	// now fix vertexTets
	for (int n = _vbt->_vertexTets.size(), i = 0; i < n; ++i) {
		if (_vbt->_vertexTets[i] > -1)
			_vbt->_vertexTets[i] = tetMap[_vbt->_vertexTets[i]];
	}
}

bool vnBccTetCutterTbb::setupBccIntersectionStructures(int maximumGridDimension)
{  // this first one must be in unmoved material coords.
	boundingBox<float> bbf;
	bbf.Empty_Box();
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i)
		bbf.Enlarge_To_Include_Point((const float(&)[3])(*_mt->vertexCoordinate(i)));
	bbf.Minimum_Corner(_vbt->_minCorner._v);
	bbf.Maximum_Corner(_vbt->_maxCorner._v);
	// In this model all grid coords are 1, all diagonal distances are sqrt(3) (~1.732), and all odd or even Cartesian distances are 2.
	_vbt->_unitSpacing = -1.0f;
	int bigDim;
	for (int i = 0; i < 3; ++i) {
		float dimSize = bbf.val[(i << 1) + 1] - bbf.val[i << 1];
		if (dimSize > _vbt->_unitSpacing) {
			_vbt->_unitSpacing = dimSize;
			bigDim = i;
		}
}
	float offset = FLT_EPSILON * (float)_vbt->_unitSpacing;
	_vbt->_minCorner -= Vec3f(offset, offset, offset);
	_vbt->_maxCorner += Vec3f(offset, offset, offset);
	double cs = _vbt->_maxCorner._v[bigDim] - _vbt->_minCorner._v[bigDim];
	_vbt->_unitSpacing = cs / maximumGridDimension;
	_vbt->_unitSpacingInv = 1.0 / _vbt->_unitSpacing;
	{
		Vec3f box = (_vbt->_maxCorner - _vbt->_minCorner)*(float)_vbt->_unitSpacingInv * 0.5f;
		for (int i = 0; i < 3; ++i) {
			_vbt->_gridSize[i] = 1 + (int)std::floor(box._v[i]);
		}
	}
	_vMatCoords.clear();
	_vMatCoords.assign(_mt->numberOfVertices(), Vec3d());
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i) {
		Vec3d Vf;
		Vf.set((const float(&)[3])*_mt->vertexCoordinate(i));
		_vMatCoords[i].set((Vf - _vbt->_minCorner)* _vbt->_unitSpacingInv);
	}
	for (int j, i = 0; i < 6; ++i) {
		// while there will be physics nodes at the boundaries of the grid, there will be no intersections there as the object is inside the boundary
		int nf = _vbt->_gridSize[i >> 1], ns = _vbt->_gridSize[((i >> 1) + 1) % 3];
		_planeSetsTbb[i].assign(nf + ns - 1, bccPlane());
		if (i & 1) {
			for (j = 0; j < ns - 1; ++j)
				_planeSetsTbb[i][j].D = -((ns - j - 1) << 1);
			for (j = 0; j < nf; ++j)
				_planeSetsTbb[i][ns - 1 + j].D = j << 1;
		}
		else {
			for (j = 0; j < nf - 1; ++j)
				_planeSetsTbb[i][j].D = (j + 1) << 1;
			for (j = 0; j < ns; ++j)
				_planeSetsTbb[i][j + nf - 1].D = (j + nf) << 1;
		}
	}
	_vertexTetLoci.clear();
	_vertexTetLoci.assign(_mt->numberOfVertices(), bccTetCentroid());
	int n = _mt->numberOfVertices();
	for (int i = 0; i < n; ++i) {
		Vec3f B;
		B.set(_vMatCoords[i]._v);
		_vbt->gridLocusToTetCentroid(B, _vertexTetLoci[i]);
		// set barycentric coordinate within that tet
		_vbt->gridLocusToBarycentricWeight(B, _vertexTetLoci[i], _vbt->_barycentricWeights[i]);
	}
	return true;
}

bool vnBccTetCutterTbb::remakeVnTets(materialTriangles *mt)
{  // uses old vbt to get material coords of mt vertices and grid data, from which it makes new vbt.
	_mt = mt;
	_vbt->_mt = mt;
	_surfaceTetFaceNumber = 0;
	int newVertexStart = _vMatCoords.size();
	_vMatCoords.resize(_mt->numberOfVertices());
	// all oldVertices should have an unchanged grid locus
	for (int n = _mt->numberOfVertices(), i = newVertexStart; i < n; ++i){  // after incisions all new vertices should have been assigned a grid locus associated with the old _vbt
		if (_vbt->getVertexTetrahedron(i) < 0)
			continue;
		Vec3f Vf;
		_vbt->vertexGridLocus(i, Vf);
		_vMatCoords[i].set(Vf);
	}
	_surfaceTetFaceNumber = 0;
	_vbt->_fixedNodes.clear();
	_vbt->_tetHash.clear();
	_vbt->_tetNodes.clear();
	_vbt->_nodeSpatialCoords = nullptr;
	if (_mt->findAdjacentTriangles(true, false))	return false;
	for (int i = 0; i < 6; ++i)	{
		// while there will be physics nodes at the boundaries of the grid, there will be no intersections there as the object is inside the boundary
		for (int n = _planeSetsTbb[i].size(), j = 0; j < n; ++j){
			_planeSetsTbb[i][j].polygons2.clear();
			_planeSetsTbb[i][j].vnTetFaces.clear();
		}
	}
	_vertexTetLoci.clear();
	bccTetCentroid emptyTC;
	emptyTC.ll = LLONG_MAX;
	_vertexTetLoci.assign(_mt->numberOfVertices(), emptyTC);
	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(_mt->numberOfVertices(), Vec3f());
	// set tet centroids and barycentric coordinates within that tet
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		if (_vbt->_vertexTets[i] < 0)
			continue;
		Vec3f Vf;
		Vf.set((double(&)[3])_vMatCoords[i]._v);
		_vbt->gridLocusToTetCentroid(Vf, _vertexTetLoci[i]);
		_vbt->gridLocusToBarycentricWeight(Vf, _vertexTetLoci[i], _vbt->_barycentricWeights[i]);
	}
	_vbt->_vertexTets.clear();
	_vbt->_vertexTets.assign(mt->numberOfVertices(), -1);
	for (int i = 0; i < 6; ++i) {
		int nS = _planeSetsTbb[i].size();
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nS), PlanePolygonMaker(_mt, i, &_vMatCoords[0], &_planeSetsTbb[i][0]));  // best to leave default grain size
	}

	// first plane set gets all the interior tetrahedral nodes.
	_nodeLoci.clear();
	_nodeLoci.reserve((_planeSetsTbb[0].size()*_vbt->_gridSize[1] * _vbt->_gridSize[2]) >> 2);  // reasonable guess
	for (int i = 0; i < 6; ++i) {
		int nS = _planeSetsTbb[i].size();
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nS), IntersectPlanePolygon(_mt, _vbt, i, &_planeSetsTbb[i][0], &_nodeLoci));
	}
	_interiorNodes.clear();  // only nodes created thus far are interior nodes
	_interiorNodes.reserve(_nodeLoci.size());
	for (int n = _nodeLoci.size(), j = 0; j < n; ++j)
		_interiorNodes.insert(std::make_pair(_nodeLoci[j], j));
	collectSurfaceTetCentersFaces();
	createVirtualNodedSurfaceTets();
	createSurfaceTetNodes();
	fillNonVnTetCenter();
	_vbt->_tetCentroids.clear();
	_vbt->_tetCentroids.reserve(_tets.size());
	_vbt->_tetNodes.clear();
	_vbt->_tetNodes.reserve(_tets.size());
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_tets.size());
	for (int n = _tets.size(), i = 0; i < n; ++i) {
		_vbt->_tetNodes.push_back(_tets[i].tetNodes);
		_vbt->_tetCentroids.push_back(_tets[i].centroid);
		_vbt->_tetHash.insert(std::make_pair(_tets[i].centroid.ll, i));
	}
	_vbt->_nodeGridLoci.clear();
	_vbt->_nodeGridLoci.assign(_nodeLoci.begin(), _nodeLoci.end());
	_nodeLoci.clear();
	// all vertices in duplicated, virtual noded tets have their tet index already assigned.  Vertices in unique tets still unassigned.
	for (int n = _vbt->_vertexTets.size(), i = 0; i < n; ++i) {
		if (_vbt->_vertexTets[i] < 0) {
			if (_vertexTetLoci[i].ll < LLONG_MAX) {  // checks for excisions
				auto vf = _vbt->_tetHash.find(_vertexTetLoci[i].ll);
				_vbt->_vertexTets[i] = vf->second;
			}
		}
	}
	_vertexTetLoci.clear();
	return true;
}

void vnBccTetCutterTbb::collectSurfaceTetCentersFaces()
{
	// Every tet face intersected by the surface will have one surface tet on either side.
	// This routine gets tet centers on either side of each intersected face, hashes them, and adds the face to that hash.
	// Also fills in any interior nodes in the face.
	_stfVec.clear();
	_surfaceTetFaceNumber = 0;
	int bigDim, bdSize = -1;
	for (int i = 0; i < 3; ++i) {
		if (bdSize < _vbt->_gridSize[i]) {
			bdSize = _vbt->_gridSize[i];
			bigDim = i;
		}
	}
	_stfVec.assign(bdSize, UMSTF());
	for (int i = 0; i < 6; ++i) {
		int nS = _planeSetsTbb[i].size();
		tbb::parallel_for(tbb::blocked_range<size_t>(0, nS), [&](const tbb::blocked_range<size_t>& r) {
			for (size_t j = r.begin(); j != r.end(); ++j) {
				int planeD = _planeSetsTbb[i][j].D;
				auto tfit = _planeSetsTbb[i][j].vnTetFaces.begin();
				auto tfEnd = _planeSetsTbb[i][j].vnTetFaces.end();
				_surfaceTetFaceNumber += (int)_planeSetsTbb[i][j].vnTetFaces.size();
				while (tfit != tfEnd) {
					// get tet centers on either side of this face
					bccTetCentroid tc0, tc1;
					int axis = i >> 1, oddFace = (tfit->first.first + tfit->first.second) & 1;
					// midpoint of face horizontal
					tc0.xyz[axis] = (i & 1) ? planeD + tfit->first.second + oddFace : planeD - tfit->first.second - oddFace;
					tc0.xyz[(axis + 1) % 3] = tfit->first.second + oddFace;
					tc0.xyz[(axis + 2) % 3] = tfit->first.first + 1;
					tc0.halfCoordAxis = axis;
					assert((tc0.xyz[axis] & 1) == (tc0.xyz[(axis + 1) % 3] & 1) && (tc0.xyz[axis] & 1) != (tc0.xyz[(axis + 2) % 3] & 1));
					tc1.xyz = tc0.xyz;
					tc1.halfCoordAxis = ((axis + 1) % 3);
					if (oddFace) {
						--tc1.xyz[(axis + 1) % 3];
						if (i & 1)
							--tc0.xyz[axis];
					}
					else {
						if ((i & 1) < 1)
							--tc0.xyz[axis];
					}
					auto fl0 = &_stfVec[tc0.xyz[bigDim] >> 1].insert(std::make_pair(tc0.ll, std::list<vnTetFace*>())).first->second;
					auto fl1 = &_stfVec[tc1.xyz[bigDim] >> 1].insert(std::make_pair(tc1.ll, std::list<vnTetFace*>())).first->second;
					auto tfLast = tfit;
					while (tfit->first == tfLast->first) {
						tfit->second.set = i;
						fl0->push_back(&tfit->second);
						fl1->push_back(&tfit->second);
						++tfit;
						if (tfit == tfEnd)
							break;
					}
				}
			}
		});
	}
}


void vnBccTetCutterTbb::tetConnectedSurface(bccTetCentroid tc, std::set<int> &triangles, std::vector<int> &vertices)
{  // Inputs tc and surface triangles seed.  Returns all the surface triangles and vertices inside the tet.
	std::set<int> verts;
	std::forward_list<int> edgeTris;
	edgeTris.assign(triangles.begin(), triangles.end());
	triangles.clear();
	// input triangles should already contain triangles intersecting tet faces. Must find all new interior triangles using their interior vertices.
	std::function<void(int)> recurseInteriorTriangles = [&](int triangle){
		if (!triangles.insert(triangle).second)
			return;
		int *tr = _mt->triangleVertices(triangle);
		for (int i = 0; i < 3; ++i){
			if (_vertexTetLoci[tr[i]] == tc){
				if (!verts.insert(tr[i]).second)
					continue;
				unsigned int adj = _mt->triAdjs(triangle)[i];
				while (adj >> 2 != triangle){
					recurseInteriorTriangles(adj >> 2);
					assert(_mt->triangleVertices(adj >> 2)[((adj & 3) + 1) % 3] == tr[i]);
					adj = _mt->triAdjs(adj >> 2)[((adj&3)+ 1) % 3];
				}
			}
		}
	};
	for (auto &et : edgeTris)
		recurseInteriorTriangles(et);
	vertices.assign(verts.begin(), verts.end());
}

void vnBccTetCutterTbb::createVirtualNodedSurfaceTets()
{
	auto createTetAssignNodes = [&](std::list<vnTetFace*> &faces, bccTetCentroid &tc, int &firstAxis, bool firstAxisUp) ->int {
		tetType tt;
		tt.centroid = tc;
		tt.tetNodes.fill(-1);
		auto tit = _tets.push_back(tt);
		int newTet = (int)(tit - _tets.begin());
		newTet <<= 2;
		for (auto &f : faces) {
			assert(f->set >> 1 != ((tc.halfCoordAxis + 1) % 3));
			if ((f->set >> 1) == firstAxis) {
				int faceSide = f->set & 1 ? 0 : 1;
				f->tetNodes[faceSide][0] = newTet;
				f->tetNodes[faceSide][1] = newTet | 1;
				if (f->set & 1)
					f->tetNodes[faceSide][2] = newTet | 2;
				else
					f->tetNodes[faceSide][2] = newTet | 3;
			}
			else {  // second axis face
				int faceSide = f->set & 1;
				if (firstAxisUp) {
					f->tetNodes[faceSide][0] = newTet | 2;
					f->tetNodes[faceSide][1] = newTet | 3;
					if (f->set & 1)
						f->tetNodes[faceSide][2] = newTet | 1;
					else
						f->tetNodes[faceSide][2] = newTet;
				}
				else {
					f->tetNodes[faceSide][0] = newTet | 3;
					f->tetNodes[faceSide][1] = newTet | 2;
					if (f->set & 1)
						f->tetNodes[faceSide][2] = newTet;
					else
						f->tetNodes[faceSide][2] = newTet | 1;
				}
			}
		}
		return (int)(newTet >> 2);
	};
	_tets.clear();
	_tets.reserve((size_t)((float)_surfaceTetFaceNumber * 2.4f));  // generous guess
// serial version
//	for (int q = 0; q < _stfVec.size(); ++q) {
//		for (auto &st : _stfVec[q]) {
	tbb::parallel_for(tbb::blocked_range<size_t>(0, _stfVec.size()), [&](const tbb::blocked_range<size_t>& r) {
		for (size_t j = r.begin(); j != r.end(); ++j) {
			for (auto &st : _stfVec[j]) {
				bccTetCentroid tc;
				int firstAxis;  // secondAxis will be halfCoordAxis
				bool firstAxisUp;
				tc.ll = st.first;
				//				tc.ll = svit->first;
				firstAxis = (tc.halfCoordAxis + 2) % 3;
				firstAxisUp = (tc.xyz[tc.halfCoordAxis] + tc.xyz[firstAxis]) & 1;
				// get groups of faces sharing common edge triangle intersects - most common surface tets
				struct edgeIntersectGroup {
					std::set<int> edgeTriangles;
					std::list<vnTetFace*> faces;
				};
				std::list<edgeIntersectGroup> eig;
				auto tfit = st.second.begin();
				while (tfit != st.second.end()) {
					if ((*tfit)->edgeTriangles.empty()) {
						++tfit;
						continue;
					}
					edgeIntersectGroup *prevFused = NULL;
					auto egit = eig.begin();
					while (egit != eig.end()) {
						if (prevFused == NULL) {
							for (auto &fet : (*tfit)->edgeTriangles) {
								if (egit->edgeTriangles.find(fet) != egit->edgeTriangles.end()) {
									egit->edgeTriangles.insert((*tfit)->edgeTriangles.begin(), (*tfit)->edgeTriangles.end());
									egit->faces.push_back(*tfit);
									prevFused = &(*egit);
									break;
								}
							}
							++egit;
						}
						else {
							auto pfit = (*tfit)->edgeTriangles.begin();
							while (pfit != (*tfit)->edgeTriangles.end()) {
								if (egit->edgeTriangles.find(*pfit) != egit->edgeTriangles.end()) {  // fuse these 2 groups
									prevFused->edgeTriangles.insert(egit->edgeTriangles.begin(), egit->edgeTriangles.end());
									prevFused->faces.splice(prevFused->faces.end(), egit->faces);
									egit = eig.erase(egit);
									break;
								}
								++pfit;
							}
							if (pfit == (*tfit)->edgeTriangles.end())
								++egit;
						}
					}
					if (prevFused == NULL) {
						eig.push_back(edgeIntersectGroup());
						eig.back().edgeTriangles.insert((*tfit)->edgeTriangles.begin(), (*tfit)->edgeTriangles.end());
						eig.back().faces.push_back(*tfit);
					}
					tfit = st.second.erase(tfit);
				}
				// have all independent edge groups and interior faces. Do any interior tests for all except single tet centroid locations.
				if (eig.size() < 2 && st.second.empty())
					createTetAssignNodes(eig.front().faces, tc, firstAxis, firstAxisUp);
				else if (eig.empty() && st.second.size() < 2) {
					std::list<vnTetFace*> oneFace;
					oneFace.splice(oneFace.end(), st.second, st.second.begin());
					createTetAssignNodes(oneFace, tc, firstAxis, firstAxisUp);
				}
				else {  // have gone as far as possible using simple common tet edge intersections.  Now must use all triangles and do
					// interior tests for possible further agglomeration or multiple tets at this location
					for (auto &eg : eig) {
						for (auto &fac : eg.faces)
							eg.edgeTriangles.insert(fac->interiorTriangles.begin(), fac->interiorTriangles.end());
					}
					for (auto &fac : st.second) {
						eig.push_back(edgeIntersectGroup());
						assert(fac->edgeTriangles.empty());
						eig.back().edgeTriangles.insert(fac->interiorTriangles.begin(), fac->interiorTriangles.end());
						eig.back().faces.push_back(fac);
					}
					st.second.clear();
					auto eit = eig.begin();
					while (eit != eig.end()) {
						std::vector<int> vertices;
						tetConnectedSurface(tc, eit->edgeTriangles, vertices);
						auto eit2 = eit;
						++eit2;
						bool recollectInterior = false;  // "dents" into a solid can lead to omissions in gathering interior
						while (eit2 != eig.end()) {
							bool notFound = true;
							for (auto tri : eit2->edgeTriangles) {
								if (eit->edgeTriangles.find(tri) != eit->edgeTriangles.end()) {
									for (auto tri2 : eit2->edgeTriangles) {  // the full interior connected search should have gotten everything
										if (eit->edgeTriangles.find(tri2) == eit->edgeTriangles.end()) {
											eit->edgeTriangles.insert(eit2->edgeTriangles.begin(), eit2->edgeTriangles.end());
											recollectInterior = true;
										}
									}
									eit->faces.splice(eit->faces.end(), eit2->faces);
									eit2 = eig.erase(eit2);
									notFound = false;
									break;
								}
							}
							if (notFound)
								++eit2;
						}
						// now have agglomerated connected faces in eit
						int thisTet = createTetAssignNodes(eit->faces, tc, firstAxis, firstAxisUp);
						// since there are likely multiple tets for this BCC locus, assign this tet to its internal vertices
						if (recollectInterior) {
							vertices.clear();
							tetConnectedSurface(tc, eit->edgeTriangles, vertices);
						}
						for (auto &v : vertices) {
							// ? comment out next assertion as overlapping undermines can violate it on a remake
							if (_vbt->_vertexTets[v] > -1 && _vbt->_vertexTets[v] != thisTet)  // should only happen in one surface tet
								assert(false);
							_vbt->_vertexTets[v] = thisTet;
						}
						eit = eig.erase(eit);
					}
				}
			}
		}
		// serial version
	//	}
		// parallel version
	});
	_stfVec.clear();
}

void vnBccTetCutterTbb::fillNonVnTetCenter()
{
	_vbt->_firstInteriorTet = _tets.size();
	// serial version
//	for (int k = 0; k < _interiorNodes.size(); ++k) {
//		tbb::blocked_range<std::size_t> r(k, k + 1);
	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, _interiorNodes.size()), [&](tbb::blocked_range<size_t> r) {
		for (int i = r.begin(); i != r.end(); ++i) {
			auto lpi = &_nodeLoci[i];
			auto lpj = _nodeLoci[i];
			lpj[2] += 2;
			auto nextZ = _interiorNodes.find(lpj);
			if (nextZ != _interiorNodes.end()) {
				int nn[4] = { -1, -1, -1, -1 };
				for (int j = 0; j < 4; ++j) {
					auto ll = _nodeLoci[i];
					ll[0] += j & 1 ? 1 : -1;
					ll[1] += j & 2 ? 1 : -1;
					++ll[2];
					auto nextNode = _interiorNodes.find(ll);
					if (nextNode != _interiorNodes.end())
						nn[j] = nextNode->second;
				}
				// See vnBccTetrahedra.h for how tetrahedral nodes are ordered.
				for (int j = 0; j < 3; j += 2) {
					if (nn[j] > -1 && nn[j + 1] > -1) {
						tetType tt;
						tt.tetNodes[0] = i;
						tt.tetNodes[1] = nextZ->second;
						tt.centroid.halfCoordAxis = 1;
						tt.centroid.xyz = _nodeLoci[i];
						tt.centroid.xyz[2] += 1;
						if (j) {
							tt.tetNodes[2] = nn[j + 1];
							tt.tetNodes[3] = nn[j];
						}
						else {
							--tt.centroid.xyz[tt.centroid.halfCoordAxis];
							tt.tetNodes[2] = nn[j];
							tt.tetNodes[3] = nn[j + 1];
						}
						if (_vbt->_tetHash.find(tt.centroid.ll) == _vbt->_tetHash.end())  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
							_tets.push_back(tt);
					}
				}
				for (int j = 0; j < 2; ++j) {
					if (nn[j] > -1 && nn[j + 2] > -1) {
						tetType tt;
						tt.tetNodes[0] = nn[j];
						tt.tetNodes[1] = nn[j + 2];
						tt.centroid.halfCoordAxis = 0;
						tt.centroid.xyz = _nodeLoci[i];
						tt.centroid.xyz[2] += 1;
						if (j) {
							tt.tetNodes[2] = i;
							tt.tetNodes[3] = i + 1;
						}
						else {
							--tt.centroid.xyz[tt.centroid.halfCoordAxis];
							tt.tetNodes[2] = i + 1;
							tt.tetNodes[3] = i;
						}
						if (_vbt->_tetHash.find(tt.centroid.ll) == _vbt->_tetHash.end())
							_tets.push_back(tt);
					}
				}
			}
			lpj = _nodeLoci[i];
			lpj[0] += 2;;
			auto in = _interiorNodes.find(lpj);
			if (in == _interiorNodes.end())
				continue;
			int x2 = in->second;
			tetType tt;
			tt.centroid.xyz = _nodeLoci[i];
			tt.centroid.xyz[0] += 1;
			for (int j = 1; j > -1; --j) {
				lpj = tt.centroid.xyz;
				lpj[2] += j ? 1 : -1;
				--lpj[1];
				in = _interiorNodes.find(lpj);
				if (in == _interiorNodes.end())
					continue;
				int y0 = in->second;
				lpj[1] += 2;
				in = _interiorNodes.find(lpj);
				if (in == _interiorNodes.end())
					continue;
				tt.centroid.halfCoordAxis = 2;
				if (j < 1)
					--tt.centroid.xyz[2];
				tt.tetNodes[0] = i;
				tt.tetNodes[1] = x2;
				tt.tetNodes[2] = j ? in->second : y0;
				tt.tetNodes[3] = j ? y0 : in->second;
				if (_vbt->_tetHash.find(tt.centroid.ll) == _vbt->_tetHash.end())
					_tets.push_back(tt);
			}
		}
	});
	_interiorNodes.clear();
//	}
}

