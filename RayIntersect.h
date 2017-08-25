#ifndef __RAY_INTERSECT_H_
#define __RAY_INTERSECT_H_

#include "CGALTypes.h"
#include <list>
#include <set>
#include "tnt/tnt_array3d.h"

class Triangle_3_With_Normal : public Triangle_3 {

private:
	Enriched_kernel::Vector_3 m_normal;

public:
	Triangle_3_With_Normal() : Triangle_3() {}

	Triangle_3_With_Normal(
		const Enriched_kernel::Point_3& p,
		const Enriched_kernel::Point_3& q,
		const Enriched_kernel::Point_3& r) : Triangle_3(p,q,r) {}

		Triangle_3_With_Normal(
			const Enriched_kernel::Point_3& p,
			const Enriched_kernel::Point_3& q,
			const Enriched_kernel::Point_3& r,
			Enriched_kernel::Vector_3 normal) : Triangle_3(p,q,r), m_normal(normal) {}

		Enriched_kernel::Vector_3 normal() const { return m_normal; }

};

typedef std::list<int> triangleIndexList_t;
typedef std::vector<Triangle_3_With_Normal *> triangleVector_t;
typedef std::vector<std::vector<int> > cellVector_t;
typedef std::set<int> triangleIndexSet_t;

//forward declaration of mesh class
class Mesh;

/**
 *	A cell holding mesh triangles, searches for intersections of triangle with a given ray
 */
class RayIntersectCell {

private:

	/** hold a list of triangles in this cell */
	triangleIndexList_t m_triangleIndexList;
	/** hold a set of indexes already in this cell so as not to add twice*/
	triangleIndexSet_t m_triangleIndexSet;
	/** a reference to the original triangles structure for comparison purposes*/
	//const triangleVector_t& m_triangleVector;

	Point_3 m_points[8];
	Plane_3 m_planes[6];

public:
	
	RayIntersectCell() {};
	virtual ~RayIntersectCell() {};

	void initGeometry(const Point_3& p1, const double xsize, const double ysize, const double zsize);

	const Point_3* points() const { return m_points; }

	const Plane_3* planes() const { return m_planes; }

	int size() const { return m_triangleIndexList.size(); }

	void addTriangle(const int triangle)
	{
		if (m_triangleIndexSet.find(triangle) == m_triangleIndexSet.end())
		{
			m_triangleIndexList.push_back(triangle);
			m_triangleIndexSet.insert(triangle);
		}		
	}

	/** Searches which triangle intersects with the ray(closest) and returns the index */
	int intersects(const triangleVector_t& triangleVector, const Ray_3& ray, const bool debugMode = false);

};

class RayIntersect {

private:
	/** vector holding triangles and normals */
	triangleVector_t m_triangleVector;

	/** the grid holding the divided triangles*/
	TNT::Array3D<RayIntersectCell*> m_grid;

	/** have the structures been initialized */
	bool m_initialized;

	double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
	double m_xspread, m_yspread, m_zspread;
	double m_diagonal;

private:
	bool findNextCell(
		const Ray_3& ray, int& x, int& y, int& z, 
		const TNT::Array3D<RayIntersectCell*>& grid,
		FILE* debugFile = NULL) const;

	double ray2planeDistance(const Plane_3& plane, const Ray_3& ray, FILE* debugFile = NULL) const;

	Point_3 index2coordinates(const int x, const int y, const int z) const;
	Point_3 index2coordinates(const int x, const int y, const int z,
							const int& gridDim1, const int& gridDim2, const int& gridDim3, 
							const double& xmin, const double& xmax, 
							const double& ymin, const double& ymax, 
							const double& zmin, const double& zmax) const;

	void coordinates2index(const Point_3& p, int& x, int& y, int& z) const;
	void coordinates2index(const Point_3& p, int& x, int& y, int& z,
						const int& gridDim1, const int& gridDim2, const int& gridDim3, 
						const double& xmin, const double& ymin, const double& zmin,
						const double& xspread, const double& yspread, const double& zspread) const;

	std::vector<int> makeCellIndexVector(const int x, const int y, const int z) const {
		std::vector<int> cell;
		cell.push_back(x);
		cell.push_back(y);
		cell.push_back(z);

		return cell;
	}

	void renderSquare(const Point_3& p1, const Point_3& p2, const Point_3& p3, const Point_3& p4) const;

public:

	RayIntersect()  : m_initialized(false) {}
	virtual ~RayIntersect();

	bool isInitialized() { return m_initialized; }

	/** Read all triangles from mesh and put them in search structure */
	void Init(const Mesh& mesh, const int gridSize);

	void Init(Mesh& mesh, TNT::Array3D<RayIntersectCell*>& grid, triangleVector_t& triangleVector,
		double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax,
		double& xspread, double& yspread, double& zspread,
		const int gridSize = 40, const bool debugMode = false);

	/** release memory */
	void Terminate(TNT::Array3D<RayIntersectCell*>& grid, triangleVector_t& triangleVector);

	/** Queries closest triangle which the ray crosses */
	double Query(
		const Ray_3& ray, /*const bool& initialized, */
		const TNT::Array3D<RayIntersectCell*>& grid, const triangleVector_t& triangleVector,
		const double& xmin, const double& ymin, const double& zmin,
		const double& xspread, const double& yspread, const double& zspread,
		const bool& removeSameDirectionIntersection, cellVector_t* debugCellVector = NULL, int* debugTriangleIndex = NULL, 
		FILE* debugFile = NULL, const bool debugMode = false) const;


	/////////////////////////////////////////////////
	// Debug rendering procedures
	/////////////////////////////////////////////////
	void renderCell(const std::vector<int>& cellIndexes) const;
	void renderCell(const int x, const int y, const int z) const;
	void renderAllCells() const;

//	void renderRay(const Ray_3& ray) const;
//	void renderRayPath(const Ray_3& ray) const;

	void generateReportFile(const char* fileName);

};

#endif