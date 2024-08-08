#ifndef __CLOCKWISE__
#define __CLOCKWISE__

#include <vector>
#include "Vec2d.h"
#include "Vec2f.h"

static int clockwise(std::vector<Vec2d>& polygon) {  // routine assumes non-self intersecting polygon.  No test is made.
	// returns 0 if polygon is collinear, 1 if clockwise, and 2 if counterclockwise.
	auto minCorner = polygon.begin();
	for (auto pit = polygon.begin(); pit != polygon.end(); ++pit) {
		if (pit->Y < minCorner->Y)
			minCorner = pit;
		else if (pit->Y == minCorner->Y && pit->X < minCorner->X)
			minCorner = pit;
		else;
	}
	auto next = minCorner;
	do{
		++next;
		if (next == polygon.end())
			next = polygon.begin();
		if (next->X != minCorner->X)
			break;
	} while (next != minCorner);
	auto prev = minCorner;
	do{
		if (prev == polygon.begin()) {
			prev = polygon.end();
			--prev;
		}
		else
			--prev;
		if (prev->X != minCorner->X && prev->X != next->X && prev->Y != next->Y)
			break;
	} while (prev != minCorner);
	if (prev == minCorner || next == minCorner)  // collinear
		return 0;
	if (minCorner->X >= prev->X && next->X >= minCorner->X)
		return(1);
	if (minCorner->X <= prev->X && next->X <= minCorner->X)
		return(2);
	// not easy so do the math
	Vec2d u = *next - *minCorner;
	Vec2d v = *prev - *minCorner;
	double det = (u.X - v.Y) * (u.Y - v.X);
	if (det > 0.0)
		return 2;
	else if (det < 0.0)
		return 1;
	else
		return 0;
}

#endif  // def __CLOCKWISE__