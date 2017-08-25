#ifndef __UTILITY_H
#define __UTILITY_H

#include "stdafx.h"

class Utility {
private:

public:
	Utility();
	virtual ~Utility();

	static Point_3 pointAddition(Point_3 p1, Point_3 p2);
	static Point_3 pointSubtraction(Point_3 p1, Point_3 p2);
	static Point_3 pointMulti(Point_3 p, double multi);
	static double pointSquaredDistance(Point_3 p1, Point_3 p2);
	static double pointDistance(Point_3 p1, Point_3 p2);
	static double pointLength(Point_3 p);
};

#endif