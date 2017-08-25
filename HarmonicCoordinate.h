#ifndef __HARMONIC_COORDINATE_H
#define __HARMONIC_COORDINATE_H

#include "stdafx.h"

#include <list>
#include <set>
#include <vector>

#include "tnt/tnt_array3d.h"

enum CellTag{UNTYPED, BOUNDARY, INTERIOR, EXTERIOR};

class GridCell {

private:
	float m_value;
	CellTag m_tag;

public:

	GridCell() { m_tag = UNTYPED; }

	//GridCell(const triangleVector_t& triangleVector) : m_triangleVector(triangleVector){}

	~GridCell() {}
	
	/// value
	float& value() {
		return m_value;
	}
	const float& value() const { 
		return m_value; 
	}
	void value(float v) { 
		m_value = v;
	}

	/// tag
	CellTag& tag() {
		return m_tag;
	}
	const CellTag& tag() const {
		return m_tag;
	}
	void tag(CellTag t) {
		m_tag = t;
	}

};


class HarmonicCoordinate {

private:
	/// vector holding triangles and normals
	triangleVector_t m_triangleVector;

	/// the grid holding the divided triangles
	TNT::Array3D<GridCell*> m_grid;

	/// have the structures been initialized
	bool m_initialized;

	double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
	double m_xspread, m_yspread, m_zspread;
	double m_diagonal;

private:
	Point_3 index2coordinates(const int x, const int y, const int z) const;

	void coordinates2index(const Point_3& p, int& x, int& y, int& z) const;

public:
	HarmonicCoordinate() {}
	~HarmonicCoordinate() {}


	//initGeometry(const Point_3& p1, const double xsize, const double ysize, const double zsize);

	//bool isInitialized() { return m_initialized; }

	void init(const Mesh& mesh);

	void laplacianSmooth();

};


#endif