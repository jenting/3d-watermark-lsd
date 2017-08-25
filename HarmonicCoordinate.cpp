//#include "stdafx.h"

#include "HarmonicCoordinate.h"

/*vector<int> HarmonicCoordinate::makeCellIndexVector(const int x, const int y, const int z) const
{
	vector<int> cell;

	cell.push_back(x);
	cell.push_back(y);
	cell.push_back(z);

	return cell;
}*/

/*void HarmonicCoordinate::initGeometry(const Point_3& p1, const double xsize, const double ysize, const double zsize)
{
	m_points[0] = p1;
	m_points[1] = Point_3(p1.x(), p1.y(), p1.z() + zsize);
	m_points[2] = Point_3(p1.x(), p1.y() + ysize, p1.z());
	m_points[3] = Point_3(p1.x(), p1.y() + ysize, p1.z() + zsize);
	m_points[4] = Point_3(p1.x() + xsize, p1.y(), p1.z());
	m_points[5] = Point_3(p1.x() + xsize, p1.y(), p1.z() + zsize);
	m_points[6] = Point_3(p1.x() + xsize, p1.y() + ysize, p1.z());
	m_points[7] = Point_3(p1.x() + xsize, p1.y() + ysize, p1.z() + zsize);

	m_planes[0] = Plane_3(m_points[7],m_points[3],m_points[1]);
	m_planes[1] = Plane_3(m_points[1],m_points[3],m_points[2]);
	m_planes[2] = Plane_3(m_points[5],m_points[1],m_points[0]);
	m_planes[3] = Plane_3(m_points[7],m_points[5],m_points[4]);
	m_planes[4] = Plane_3(m_points[4],m_points[0],m_points[2]);
	m_planes[5] = Plane_3(m_points[6],m_points[2],m_points[3]);
}*/

Point_3 HarmonicCoordinate::index2coordinates(const int x, const int y, const int z) const
{
	Point_3 p( m_xmin + ((double) x / m_grid.dim1() * (m_xmax-m_xmin)), m_ymin + ((double) y / m_grid.dim2() * (m_ymax-m_ymin)), m_zmin + ((double) z / m_grid.dim3() * (m_zmax-m_zmin)));

	return p;
}

void HarmonicCoordinate::coordinates2index(const Point_3& p, int& x, int& y, int& z) const
{
	x = ((double) (p.x() - m_xmin) * m_xspread * m_grid.dim1());
	x = qMin(x,m_grid.dim1()-1);
	y = ((double) (p.y() - m_ymin) * m_yspread * m_grid.dim2());
	y = qMin(y,m_grid.dim2()-1);
	z = ((double) (p.z() - m_zmin) * m_zspread * m_grid.dim3());
	z = qMin(z,m_grid.dim3()-1);
}

void HarmonicCoordinate::init(const Mesh& mesh)
{
	//step 0 - make sure all data structures are empty
	if (m_initialized) {
		for (int x=0;x<m_grid.dim1();x++) {
			for (int y=0;y<m_grid.dim2();y++) {
				for (int z=0;z<m_grid.dim3();z++) {
					delete m_grid[x][y][z];
				}
			}
		}
	}

	// extended boundary grid cells

	const float delta = 1e-5;

	m_xmin = mesh.xmin() - delta;
	m_ymin = mesh.ymin() - delta;
	m_zmin = mesh.zmin() - delta;

	m_xmax = mesh.xmax() + delta;
	m_ymax = mesh.ymax() + delta;
	m_zmax = mesh.zmax() + delta;

	m_xspread = 1 / (m_xmax - m_xmin + 1e-5);
	m_yspread = 1 / (m_ymax - m_ymin + 1e-5);
	m_zspread = 1 / (m_zmax - m_zmin + 1e-5);
	
	m_diagonal = sqrt(pow(m_xmax-m_xmin, 2)+pow(m_ymax-m_ymin, 2) + pow(m_zmax-m_zmin, 2));

	float x_range = m_xmax - m_xmin;
	float y_range = m_ymax - m_ymin;
	float z_range = m_zmax - m_zmin;

	float max_range;
	
	if(x_range > y_range)
		max_range = x_range;
	else 
		max_range = y_range;
	if(max_range < z_range)
		max_range = z_range;

	printf("xrange: %g\n", x_range);
	printf("yrange: %g\n", y_range);
	printf("zrange: %g\n", z_range);

	printf("max range: %g\n", max_range);

	const int gridSize = 128;

	float cell_range = max_range / (float) gridSize;

	printf("cell range: %g\n", cell_range);

	m_grid = TNT::Array3D<GridCell*>(gridSize, gridSize, gridSize);
	for (int x=0;x<gridSize;x++) {
		for (int y=0;y<gridSize;y++) {
			for (int z=0;z<gridSize;z++) {
				m_grid[x][y][z] = new GridCell();
				//m_grid[x][y][z]->initGeometry(index2coordinates(x,y,z), cell_range, cell_range, cell_range);
			}
		}
	}

	// step 2a - tag all cells as UNTYPED
	// has been done in GridCell() constructor
	
	// step 2b - scan-convert boundary conditions into the grid, marking each scan converted cell with the BOUNDARY tag
	//m_grid[0][0][0]

	// step 2c - flood fill the exterior, marking each visited cell with the EXTERIOR tag. The flood fill recursion stops when BOUNDARY tags are reached
	

	// step 2d - mark remaining UNTYPED cells as INTERIOR with harmonic coordinate value equal to 0
	for (int x = 0; x < gridSize; x++) {
		for (int y = 0; y < gridSize; y++) {
			for (int z = 0; z < gridSize; z++) {
				if (m_grid[x][y][z]->tag() == UNTYPED) {
					m_grid[x][y][z]->tag(INTERIOR);
					m_grid[x][y][z]->value(0);
				}
			}
		}
	}

	// step 3 - Laplacian smooth
	laplacianSmooth();

	//step 4
	m_initialized = true;

}

void HarmonicCoordinate::laplacianSmooth()
{
	// create a template grid cell

	const int gridSize = 128;

	//TNT::Array3D<float*> temp_grid = TNT::Array3D<float*>(gridSize, gridSize, gridSize);
	TNT::Array3D<float> temp_grid = TNT::Array3D<float>(gridSize, gridSize, gridSize);
	for (int x=0;x<gridSize;x++)
		for (int y=0;y<gridSize;y++)
			for (int z=0;z<gridSize;z++)
				temp_grid[x][y][z] = 0.0;


	const float termin_thres = 1e-5;
	while(true) {
		// smoothing with 6-connected neighbors

		for (int x=0;x<gridSize;x++) {
			for (int y=0;y<gridSize;y++) {
				for (int z=0;z<gridSize;z++) {
					if (m_grid[x][y][z]->tag() == INTERIOR) {

						float avg = 0;

						avg += m_grid[x-1][y    ][z    ]->value();
						avg += m_grid[x    ][y-1][z    ]->value();
						avg += m_grid[x    ][y    ][z-1]->value();
						avg += m_grid[x+1][y    ][z    ]->value();
						avg += m_grid[x    ][y+1][z    ]->value();
						avg += m_grid[x    ][y     ][z+1]->value();
						
						avg /= 6;
						
						temp_grid[x][y][z] = avg;
					
					}
				}
			}
		}

		// storage temp_grid value into m_grid value, and calculate average cell change

		float avg_change = 0;
		int num_of_interior = 0;

		for (int x=0;x<gridSize;x++) {
			for (int y=0;y<gridSize;y++) {
				for (int z=0;z<gridSize;z++) {
					if (m_grid[x][y][z]->tag() == INTERIOR) {

						num_of_interior++;

						float new_value = temp_grid[x][y][z];
						
						avg_change += fabs(m_grid[x][y][z]->value() - new_value);

						m_grid[x][y][z]->value(new_value);

					}
				}
			}
		}

		avg_change /= num_of_interior;

		if (avg_change < termin_thres)
			break;
	}

}