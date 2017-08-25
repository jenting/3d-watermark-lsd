#include "stdafx.h"

#include <q3filedialog.h>
#include <q3url.h>
#include <qlabel.h>
#include <QSet>

#include "AppManager.h"
#include "ColorSchemeManager.h"
#include "CylLine.h"
#include "FeaturePainters.h"
#include "RayIntersect.h"
#include "WorldManager.h"

#ifndef MIN
#define MIN(a,b) (a<b?a:b)
#endif

void RayCellsPainter::paint(Mesh* mesh, MeshSelection* meshSelection, MeshSurfaceGraph* msg)
{
	if (!m_rayIntersect.isInitialized())
		m_rayIntersect.Init(*mesh, 10);
	
	m_rayIntersect.renderAllCells();
}


double compute_dihedral_angle(
	Mesh::Facet_const_handle& f1,
	Mesh::Facet_const_handle& f2,
	Mesh::Halfedge_const_handle& edgeBetween)
{
	float x = f1->normal() * f2->normal();
	// fix rounding errors which result in NANs
	if (x > 1.0)
		x = 1.0;
	if (x < -1.0)
		x = -1.0;
	double angle = acos(x);

	if (angle != angle)
		__asm int 3;

	//normalize angle	
	angle = angle * M_1_PI;
	if ((angle > 1.0) || (angle < 0.0))
		__asm int 3

	//angle = MYMIN(1, angle*1.3);
// 	Enriched_kernel::Segment_3 s(
// 		Mesh::Point_3(f1->normal()[0], f1->normal()[1], f1->normal()[2]), 
// 		Mesh::Point_3(f2->normal()[0], f2->normal()[1], f2->normal()[2]));
// 
// 	double dihedral = sqrt(s.squared_length());

	//////////////////////////////////////////////////////////////////////////
	// Concave angles
	//////////////////////////////////////////////////////////////////////////
	double concavityMultiplier = 0.1;
	Enriched_kernel::Plane_3 plane(edgeBetween->vertex()->point(), edgeBetween->next()->vertex()->point(), edgeBetween->next()->next()->vertex()->point());
	Enriched_kernel::Point_3 p = edgeBetween->opposite()->next()->vertex()->point();
	if (plane.oriented_side(p) == CGAL::ON_POSITIVE_SIDE)
	{
		concavityMultiplier = 1;		
	}

	//return (dihedral/2) * concavityMultiplier;
	return angle * concavityMultiplier;
}

/*
double CurvatureVertexPainter::getval(const Mesh::Vertex_const_handle& v)
{
	switch(m_type) {
		case 0:
			return v->curvature().kmin();
			break;
		case 1:
			return v->curvature().kmax();
			break;
		case 2:
			return v->curvature().kgauss();
			break;
		case 3:
			return v->curvature().kmean();
			break;	
		default:
			return 0.0;
	}	
}

bool CurvatureVertexPainter::fillBuffer( 
	const RenderingParams* renderingParams, 
	Mesh* mesh, MeshSurfaceGraph* msg,
	MeshSelection* meshSelection, 
	GLubyte* colorBuffer)
{
	double cmin = FLT_MAX;
	double cmax = -FLT_MAX;
	double inv_range;
	for (Mesh::Vertex_const_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end();; vit++) {
		double val = getval(vit);
		if (val<cmin) cmin = val;
		if (val>cmax) cmax = val;
	}
	inv_range = 1 / (cmax-cmin);
	vit = mesh->vertices_begin();

	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);

	int dataIndex = 0;
	for (;vit != vit_end; vit++) {
		double val = getval(vit);
		QColor c = cs.gradient((val-cmin) * inv_range);
		colorBuffer[dataIndex] = c.red();
		colorBuffer[dataIndex+1] = c.green();
		colorBuffer[dataIndex+2] = c.blue();
		colorBuffer[dataIndex+3] = renderingParams->m_alpha;
		dataIndex += 4;
	}

	return true;
}*/

unsigned int SdfVerticesPainter::requiredBufferSize(Mesh* mesh)
{
	return mesh->size_of_vertices() * 4;
}

bool SdfVerticesPainter::fillBuffer(
						const RenderingParams* renderingParams,
						Mesh* mesh,
						MeshSurfaceGraph* msg,
						MeshSelection* meshSelection,
						GLubyte* colorBuffer)
{
	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);

	int dataIndex = 0;
	for (Mesh::Vertex_const_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end(); vit++) {
		QColor c = cs.gradient(vit->volumeSDF());
		colorBuffer[dataIndex] = c.red();
		colorBuffer[dataIndex+1] = c.green();
		colorBuffer[dataIndex+2] = c.blue();
		colorBuffer[dataIndex+3] = renderingParams->m_alpha;
		dataIndex += 4;
	}

	return true;
}

bool SdfVerticesPainter::paintVertex(
						 const RenderingParams* renderingParams,
						 const Enriched_Mesh::Vertex_const_handle& pVertex,
						 GLubyte* color)
{
	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);
	QColor c = cs.gradient(pVertex->volumeSDF());
	color[0] = c.red();
	color[1] = c.green();
	color[2] = c.blue();
	color[3] = renderingParams->m_alpha;

	return true;
}


bool SdfVerticesPainter::paintVertexInFacet(
								const RenderingParams* renderingParams,
								const Enriched_Mesh::Facet_const_handle& pFacet,
								const Enriched_Mesh::Halfedge_const_handle& pHalfedge,
								GLubyte* color)
{
	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);
	QColor c = cs.gradient(pHalfedge->vertex()->volumeSDF());
	color[0] = c.red();
	color[1] = c.green();
	color[2] = c.blue();
	color[3] = 255;

	return true;
}

//////////////////////////////////////////////////////////////////////////

double SdfFacetsPainter::logit(const double& value, const int& range)
{
	return logf(value * range + 1)/logf(range+1);
}

unsigned int SdfFacetsPainter::requiredBufferSize(Mesh* mesh) {
	return mesh->size_of_facets() * 12;
}

bool SdfFacetsPainter::fillBuffer(
	const RenderingParams* renderingParams,
	Mesh* mesh,
	MeshSurfaceGraph* msg,
	MeshSelection* meshSelection,
	GLubyte* colorBuffer)
{
	/*double featureMin = FLT_MAX;
	double featureMax = 0.0;
	double inverseRange = 0.0;
	if (m_normalize) {
		for (Mesh::Facet_const_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) {
			double feature = fit->volumeSDF();
			if (feature	<featureMin) featureMin = feature;
			if (feature >featureMax) featureMax = feature;
		}
		inverseRange = 1 / (featureMax - featureMin);
	}*/

	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);

	int dataIndex = 0;
	for (Mesh::Facet_const_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) {
		double feature = fit->volumeNSDF();
		/*if (m_normalize) {
			feature = logit((feature - featureMin) * inverseRange, 5);
		}*/

		QColor c = cs.gradient(feature);
		for (int i=0;i<3;i++) {
			colorBuffer[dataIndex] = c.red();
			colorBuffer[dataIndex+1] = c.green();
			colorBuffer[dataIndex+2] = c.blue();
			colorBuffer[dataIndex+3] = renderingParams->m_alpha;
			dataIndex += 4;
		}		
	}

	return true;
}

bool SdfFacetsPainter::paintVertex(
	const RenderingParams* renderingParams,
	const Enriched_Mesh::Vertex_const_handle& pVertex,
	GLubyte* color)
{
	color[0] = 0;
	color[1] = 0;
	color[2] = 0;
	color[3] = 255;

	return true;
}


bool SdfFacetsPainter::paintVertexInFacet(
	const RenderingParams* renderingParams,
	const Enriched_Mesh::Facet_const_handle& pFacet,
	const Enriched_Mesh::Halfedge_const_handle& pHalfedge,
	GLubyte* color)
{
	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);
	QColor c = cs.gradient(pFacet->volumeSDF());
	color[0] = c.red();
	color[1] = c.green();
	color[2] = c.blue();
	color[3] = 255;

	return true;
}

bool PartitionPainter::fillBuffer(
	const RenderingParams* renderingParams,
	Mesh* mesh,
	MeshSurfaceGraph* msg,
	MeshSelection* meshSelection,
	GLubyte* colorBuffer)
{	
	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);

	int dataIndex = 0;
	for (Mesh::Facet_const_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) {
		int index = fit->cluster();
		int cluster = index % 6;

		//QColor c = cs.gradient(feature);
		//QColor c = cs.feature;
		switch(cluster)
		{
			case 0:
					for (int i=0;i<3;i++)
					{
						colorBuffer[dataIndex] = 255;
						colorBuffer[dataIndex+1] = 0;
						colorBuffer[dataIndex+2] = 0;
						colorBuffer[dataIndex+3] = renderingParams->m_alpha;
						dataIndex += 4;
					}
					break;
			case 1:
					for (int i=0;i<3;i++)
					{
						colorBuffer[dataIndex] = 0;
						colorBuffer[dataIndex+1] = 255;
						colorBuffer[dataIndex+2] = 0;
						colorBuffer[dataIndex+3] = renderingParams->m_alpha;
						dataIndex += 4;
					}
					break;
			case 2:
					for (int i=0;i<3;i++)
					{
						colorBuffer[dataIndex] = 0;
						colorBuffer[dataIndex+1] = 0;
						colorBuffer[dataIndex+2] = 255;
						colorBuffer[dataIndex+3] = renderingParams->m_alpha;
						dataIndex += 4;
					}
					break;
			case 3:
					for (int i=0;i<3;i++)
					{
						colorBuffer[dataIndex] = 255;
						colorBuffer[dataIndex+1] = 255;
						colorBuffer[dataIndex+2] = 0;
						colorBuffer[dataIndex+3] = renderingParams->m_alpha;
						dataIndex += 4;
					}
					break;
			case 4:
					for (int i=0;i<3;i++) 
					{
						colorBuffer[dataIndex] = 0;
						colorBuffer[dataIndex+1] = 255;
						colorBuffer[dataIndex+2] = 255;
						colorBuffer[dataIndex+3] = renderingParams->m_alpha;
						dataIndex += 4;
					}
					break;
			case 5:
					for (int i=0;i<3;i++)
					{
						colorBuffer[dataIndex] = 255;
						colorBuffer[dataIndex+1] = 0;
						colorBuffer[dataIndex+2] = 255;
						colorBuffer[dataIndex+3] = renderingParams->m_alpha;
						dataIndex += 4;
					}
					break;
		}
	}

	return true;
}

bool SdfDifferencesFacetsPainter::fillBuffer(
	const RenderingParams* renderingParams,
	Mesh* mesh,
	MeshSurfaceGraph* msg,
	MeshSelection* meshSelection,
	GLubyte* colorBuffer)
{	
	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);

	int dataIndex = 0;
	for (Mesh::Facet_const_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) {
		double feature = fit->volumeNSDF();

		//double sum = 0.0;
		//double sumweight = 0.0;
		double maxNeighborDiff = 0.0;
		Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();
		do {
			Mesh::Facet_const_handle neighbor = c->opposite()->facet();
			if (neighbor == NULL) continue;
			//double nvalue = neighbor->volumeSDF();
			double nvalue = neighbor->volumeNSDF();
			//sum += nvalue;
			//sumweight += 1.0;
			if (fabs(nvalue-feature) > maxNeighborDiff) maxNeighborDiff = fabs(nvalue-feature);
		} while (++c != fit->facet_begin());

		/*if (sumweight != 0.0)
			sum /= sumweight;*/

		//double diff = sqrt(fabs(feature - sum))/* / ((feature + sum)/2)*/;
		//double diff = 1 - (logf(sqrt(maxNeighborDiff) * 2 + 1) / logf(3));
		double diff = MIN(maxNeighborDiff / (feature + 1e-4), 1);

		QColor col = cs.gradient(diff);
		for (int i=0;i<3;i++) {
			colorBuffer[dataIndex] = col.red();
			colorBuffer[dataIndex+1] = col.green();
			colorBuffer[dataIndex+2] = col.blue();
			colorBuffer[dataIndex+3] = renderingParams->m_alpha;
			dataIndex += 4;
		}		
	}

	return true;
}

bool DeformationDegreeFacetsPainter::fillBuffer(
	const RenderingParams* renderingParams,
	Mesh* mesh,
	MeshSurfaceGraph* msg,
	MeshSelection* meshSelection,
	GLubyte* colorBuffer)
{
	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);

	int dataIndex = 0;
	for (Mesh::Facet_const_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) {
		double feature = fit->normalizeDeformationDegree();

		QColor c = cs.gradient(feature);
		for (int i=0;i<3;i++) {
			colorBuffer[dataIndex] = c.red();
			colorBuffer[dataIndex+1] = c.green();
			colorBuffer[dataIndex+2] = c.blue();
			colorBuffer[dataIndex+3] = renderingParams->m_alpha;
			dataIndex += 4;
		}
	}

	return true;
}

bool DeformationDegreeDifferencesFacetsPainter::fillBuffer(
	const RenderingParams* renderingParams,
	Mesh* mesh,
	MeshSurfaceGraph* msg,
	MeshSelection* meshSelection,
	GLubyte* colorBuffer)
{	
	ColorScheme cs = ColorSchemeManager::GetScheme(renderingParams->m_colorSchemeName);

	int dataIndex = 0;
	for (Mesh::Facet_const_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) {
		double feature = fit->normalizeDeformationDegree();

		double maxNeighborDiff = 0.0;
		Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();
		do {
			Mesh::Facet_const_handle neighbor = c->opposite()->facet();
			if (neighbor == NULL) continue;
			
			double nvalue = neighbor->normalizeDeformationDegree();
			
			if (fabs(nvalue-feature) > maxNeighborDiff) maxNeighborDiff = fabs(nvalue-feature);
		} while (++c != fit->facet_begin());

		double diff = MIN(maxNeighborDiff / (feature + 1e-4), 1);

		QColor col = cs.gradient(diff);
		for (int i=0;i<3;i++) {
			colorBuffer[dataIndex] = col.red();
			colorBuffer[dataIndex+1] = col.green();
			colorBuffer[dataIndex+2] = col.blue();
			colorBuffer[dataIndex+3] = renderingParams->m_alpha;
			dataIndex += 4;
		}
	}

	return true;
}

/*
void RayIntersect::renderSquare(const Point_3& p1, const Point_3& p2, const Point_3& p3, const Point_3& p4) const
{
	//render a cell nicely
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_QUADS);

	glVertex3f(p1.x(), p1.y(), p1.z());
	glVertex3f(p2.x(), p2.y(), p2.z());
	glVertex3f(p3.x(), p3.y(), p3.z());
	glVertex3f(p4.x(), p4.y(), p4.z());

	glEnd();
}

	Mesh::Vertex_handle v = sortedVertex[pointIndex];

			const Point_3 originalPoint = v->point();		
			const double originalNSDF = v->volumeNSDF();

			//////////////////////////////////////////////////////////////////////////
			// Get min and max distance of it`s 1-ring neighbor
			//////////////////////////////////////////////////////////////////////////
			double dist;
			double minDist = FLT_MAX;
			double maxDist = -FLT_MAX;
			double avgDist = 0;
			Mesh::Halfedge_around_vertex_circulator c = v->vertex_begin();
			do {
				Mesh::Vertex_handle vNeighbor = c->opposite()->vertex();
				Point_3 pNeighborPoint = vNeighbor->point();

				dist = distance(originalPoint, pNeighborPoint);

				avgDist += dist;

				if (dist < minDist) minDist = dist;
				if (dist > maxDist) maxDist = dist;
			} while (++c != v->vertex_begin());

			avgDist /= v->degree();
			
			double searchSpaceDist;
			switch (searchSpace)
			{
				case MIN_DIST: searchSpaceDist = minDist; break;
				case AVG_DIST: searchSpaceDist = avgDist; break;
				case MAX_DIST: searchSpaceDist = maxDist; break;
			}

			//printf("minDist = %g, avgDist = %g, maxDist = %g\n", minDist, avgDist, maxDist);
		
			//////////////////////////////////////////////////////////////////////////
			// Initialization
			//////////////////////////////////////////////////////////////////////////
			#pragma omp parallel num_threads(threadsNumber)
			{
				#pragma omp for
				for (int j = 0; j < particleNumbers; j++) /// Initialize Points ///
				{
					Point_3 randomPoint = randomPointInSphere(originalPoint, searchSpaceDist);
					
					deformedPointArray[j] = randomPoint;

					position[0][j] = individualOptimalPosition[0][j] = randomPoint.x();
					position[1][j] = individualOptimalPosition[1][j] = randomPoint.y();
					position[2][j] = individualOptimalPosition[2][j] = randomPoint.z();

					velocity[0][j] = randomNumber(-searchSpaceDist / 20.0, searchSpaceDist / 20.0);
					velocity[1][j] = randomNumber(-searchSpaceDist / 20.0, searchSpaceDist / 20.0);
					velocity[2][j] = randomNumber(-searchSpaceDist / 20.0, searchSpaceDist / 20.0);

					double newRadius = distance(originalPoint, randomPoint);

					if (newRadius > searchSpaceDist)
						cout << "\tError on initialize points: " << newRadius << "," << searchSpaceDist << endl;
				}
			} /// end OpenMP
*/

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////


SdfFacetDebugPainter::SdfFacetDebugPainter() : m_cones(3), m_coneSeperationDegree(5), m_raysInCone(8), m_gridSize(10) {}

void SdfFacetDebugPainter::paint(Mesh* mesh, MeshSelection* meshSelection, MeshSurfaceGraph* msg)
{
	if (meshSelection == NULL || meshSelection->selectedFacet() == -1)
		return;

	glDisable(GL_LIGHTING);

	Mesh::Facet_handle f = mesh->findFacet(meshSelection->selectedFacet());
	Point_3 p = mesh->getFacetCenter(f);
	Vector_3 n = f->normal();

	const double m_coneSeparationRadian = (double) m_coneSeperationDegree * DEG2RAD;
	std::list<std::pair<Ray_3, double> > rays;

	if (!m_rayIntersect.isInitialized())
	{
		//m_rayIntersect.Init(*mesh, 10);

		m_rayIntersect.Init(*mesh, m_grid, m_triangleVector,
						m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax,
						m_xspread, m_yspread, m_zspread, m_gridSize);
	}

	CalculateSDF csdf;
	csdf.computeVolumeSDF(p, n, m_gridSize, 
					m_cones, m_coneSeperationDegree, m_coneSeparationRadian, m_raysInCone, 
					true, true, false, 
					m_grid, m_triangleVector, m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax,
					m_xspread, m_yspread, m_zspread,
					false, &rays);

	//float dl = mesh->diagonalLength();

	//CylLine cyl(dl * 0.0015, 10);

	glBegin(GL_LINES);
	glLineWidth(5);

	//glColor3f(1.0, 0.1, 0.1);

	double diagonalLength = mesh->diagonalLength();

	for (std::list<std::pair<Ray_3, double> >::iterator it = rays.begin(); it != rays.end(); it++)
	{
		if(it->second == diagonalLength)
		{
			//glColor3f(0.85,0.1,0.05);
		}
		else
		{
			glColor3f(0.1,0.9,0.1);
		}

		Point_3 p1 = it->first.source();
		Point_3 p2 = p1 + it->first.to_vector()*it->second;

		glVertex3f(p1.x(), p1.y(), p1.z());
		glVertex3f(p2.x(), p2.y(), p2.z());

		//cyl(p1.x(), p1.y(), p1.z(),p2.x(), p2.y(), p2.z());
	}

	glEnd();
}

/*QString SdfFacetDebugPainter::command(const int key)
{
	switch(key) {
		case 0:
			m_cones--; break;
		case 1:
			m_cones++; break;
		case 2:
			m_angle--; break;
		case 3:
			m_angle++; break;
		case 4:
			m_raysInCone--; break;
		case 5:
			m_raysInCone++; break;
		case 9:
			m_cones = 3;
			m_angle = 20;
			m_raysInCone = 6;
	}

	m_cones = (m_cones<0? 0 : m_cones);
	m_angle = (m_angle<10? 10 : m_angle);
	m_raysInCone = (m_raysInCone<2? 2 : m_raysInCone);

	printf("Cones %d | Angle %d | Rays %d", m_cones, m_angle, m_raysInCone);

	QString temp;
	temp.sprintf("Cones %d | Angle %d | Rays %d", m_cones, m_angle, m_raysInCone);
	return temp;
}*/


//////////////////////////////////////////////////////////////////////////


SdfVertexDebugPainter::SdfVertexDebugPainter() : m_cones(4), m_coneSeperationDegree(5), m_raysInCone(16), m_gridSize(10) {}

void SdfVertexDebugPainter::paint(Mesh* mesh, MeshSelection* meshSelection, MeshSurfaceGraph* msg)
{
	//if (meshSelection == NULL || meshSelection->selectedVertex() == -1)
	//	return;
	
	if (meshSelection == NULL)
		return;

/*
	glDisable(GL_LIGHTING);

	Mesh::Vertex_handle v = mesh->findVertex(100);
	Point_3 p = v->point();
	Vector_3 n = v->normal();

	const double m_coneSeparationRadian = m_coneSeperationDegree * DEG2RAD;

	double debugMedian;
	double debugRange;
	std::list<std::pair<Ray_3, double> > rays;	

	mesh->computeVolumeSDF(
		p, n, 
		m_gridSize, m_cones, m_coneSeperationDegree, m_coneSeparationRadian, m_raysInCone, false, true, 
		false, &rays, &debugMedian, &debugRange);

	//float dl = mesh->diagonalLength();

	//CylLine cyl(dl * 0.0015, 10);

	glBegin(GL_LINES);
	glLineWidth(5);

	//glColor3f(1.0, 0.1, 0.1);

	double diagonalLength = mesh->diagonalLength();

	for (std::list<std::pair<Ray_3, double> >::iterator it = rays.begin(); it != rays.end(); it++)
	{
		//if(it->second == diagonalLength)
		//	glColor3f(0.85,0.1,0.05);
		//else
		//	glColor3f(0.1,0.9,0.1);

		if (fabs(it->second-debugMedian) < debugRange*1.5)
			glColor3f(0.1,0.9,0.1); // green
		else
			glColor3f(0.85,0.1,0.05); // red

		Point_3 p1 = it->first.source();
		Point_3 p2 = p1 + it->first.to_vector()*it->second;

		glVertex3f(p1.x(), p1.y(), p1.z());
		glVertex3f(p2.x(), p2.y(), p2.z());

		//cyl(p1.x(), p1.y(), p1.z(),p2.x(), p2.y(), p2.z());
	}

	glEnd();
*/
}


//////////////////////////////////////////////////////////////////////////


double DihedralDebugPainter::calcWeight(Mesh::Facet_const_handle& f1, Mesh::Facet_const_handle& f2, Mesh::Halfedge_const_handle& h)
{
	double dihedral = compute_dihedral_angle(f1,f2,h);
	return dihedral;
}

void DihedralDebugPainter::paint(Mesh* mesh, MeshSelection* MeshSelection, MeshSurfaceGraph* msg)
{
	if (!mesh) return;

	glDisable(GL_LINE_SMOOTH);
	glDisable(GL_BLEND);
	// enable polygon offset
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(3.0f,1.0f);
	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	ColorScheme cs;
//	cs.colors.push_back(Qt::red);
	cs.colors.push_back(Qt::green);
	cs.colors.push_back(Qt::blue);	

	glLineWidth(2);
	glBegin(GL_LINES);
	for (Mesh::Edge_const_iterator eit = mesh->edges_begin(); eit != mesh->edges_end(); eit++)
	{
		Mesh::Facet_const_handle f1 = eit->facet();
		Mesh::Facet_const_handle f2 = eit->opposite()->facet();
		if (f1 == NULL || f2 == NULL) continue;

		double weight = calcWeight(f1,f2, eit);

		if (weight > 0.1) 
		{
			QColor c = cs.gradient(weight);
			glColor3ub(c.red(), c.green(), c.blue());

			glVertex3f(eit->vertex()->point()[0], eit->vertex()->point()[1], eit->vertex()->point()[2]);
			glVertex3f(eit->opposite()->vertex()->point()[0], eit->opposite()->vertex()->point()[1], eit->opposite()->vertex()->point()[2]);
		}
	}

	glEnd();
}

//////////////////////////////////////////////////////////////////////////

/*
CagePainter::CagePainter() {}

void CagePainter::paint(Mesh* cage_mesh, MeshSelection* MeshSelection, MeshSurfaceGraph* msg)
{
	glDisable(GL_LIGHTING);

	glLineWidth(3);

	glBegin(GL_TRIANGLE_STRIP);

	Mesh::Vertex_const_iterator vit = cage_mesh->vertices_begin();
	Mesh::Vertex_const_iterator vit_end = cage_mesh->vertices_end();
	for (; vit != vit_end; vit++)
	{		
		glColor3f(0.1,0.9,0.1);

		Point_3 p = vit->point();
		glVertex3f(p.x(), p.y(), p.z());
	}

	glEnd();
}
*/