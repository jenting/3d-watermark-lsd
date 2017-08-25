#include "Utility.h"

Utility::Utility()
{
}

Utility::~Utility()
{
}

Point_3 Utility::pointAddition(Point_3 p1, Point_3 p2)
{
	Point_3 p_add(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());
	return p_add;
}

Point_3 Utility::pointSubtraction(Point_3 p1, Point_3 p2)
{
	Point_3 p_sub(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
	return p_sub;
}

Point_3 Utility::pointMulti(Point_3 p, double multi)
{
	Point_3 p_multi(p.x() * multi, p.y() * multi, p.z() * multi);
	return p_multi;
}

double Utility::pointSquaredDistance(Point_3 p1, Point_3 p2)
{
	return pow((p1.x() - p2.x()), 2) + pow((p1.y() - p2.y()), 2) + pow((p1.z() - p2.z()), 2);
}

double Utility::pointDistance(Point_3 p1, Point_3 p2)
{
	double dist = pointSquaredDistance(p1, p2);
	return sqrt(dist);
}

double Utility::pointLength(Point_3 p)
{
	return sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z());
}