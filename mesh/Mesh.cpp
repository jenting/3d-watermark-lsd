#include "stdafx.h"

#include "Mesh.h"
//#include "curvature_estimator.h"

#include <list>

#include <QFile>
#include <QString>
#include <q3cstring.h>
#include <QTextStream>
#include <QBuffer>

#include <omp.h>

const float Mesh::cm_colorWhite[3] = {0.95,0.95,0.95};
const float Mesh::cm_colorGray[3] = {0.5,0.5,0.5};
const float Mesh::cm_colorBlack[3] = {0.0,0.0,0.0};

Mesh::Mesh() :
	//m_2dSearchStructure(NULL),
	m_totalSurface(FLT_MAX)/*,
	m_triangleSearch(NULL)*/
{

}

Mesh::~Mesh()
{

}

/*int Mesh::translatePart(int from, int to)
{
	int count = 0;
	for (Mesh::Facet_iterator it = facets_begin(); it!= facets_end(); ++it)
	{
		if (it->part() == from)
		{
			++count;
			it->part(to);
		}
	}
	return count;
}*/

double Mesh::computeTotalVolume()
{
	vector<Point_3> points(3);
	int index;
	double totalVolume = 0.0;

	Mesh::Facet_const_iterator fit = facets_begin();
	Mesh::Facet_const_iterator fit_end = facets_end();
	for (; fit != fit_end; fit++)
	{
		assert(fit->size() == 3);

		index = 0;

		Mesh::Halfedge_around_facet_const_circulator pHalfedge = fit->facet_begin();
		Mesh::Halfedge_around_facet_const_circulator end = pHalfedge;
		CGAL_For_all(pHalfedge, end)
		{
			points[index] = pHalfedge->vertex()->point();

			index++;
		}

		double det = (points[1].x() - points[0].x()) * (points[2].y() - points[0].y()) - (points[2].x() - points[0].x()) * (points[1].y() - points[0].y());
		double localVolume = (points[0].z() + points[1].z() + points[2].z()) * det;
		totalVolume += localVolume;
	}

	return totalVolume;
}

void Mesh::computeTriangleArea()
{
	double* edge_length = new double[3];
	int index;
	double sum;

	for (Mesh::Facet_iterator fit = facets_begin(); fit != facets_end(); fit++)
	{
		assert(fit->size() == 3);

		index = 0;
		sum = 0;

		Mesh::Halfedge_around_facet_circulator pHalfedge = fit->facet_begin();
		Mesh::Halfedge_around_facet_circulator end = pHalfedge;
		CGAL_For_all(pHalfedge,end)
		{
			sum += pHalfedge->length();
			edge_length[index] = pHalfedge->length();

			index++;
		}

		double s = sum / 2;
		double ar = sqrt(s * (s - edge_length[0]) * (s - edge_length[1]) * (s - edge_length[2]));

		fit->area(ar);
	}

	delete [] edge_length;
}

double Mesh::computeTriangleSurfaces()
{
	m_totalSurface = 0.0;
	for (Mesh::Facet_iterator it = facets_begin(); it != facets_end(); it++)
	{
		it->surface(sqrt(createTriangle(it->facet_begin()).squared_area()));
		m_totalSurface += it->surface();
	}

	return m_totalSurface;
}

double Mesh::computeSurfaceTotal()
{
	m_totalSurface = 0.0;
	for (Mesh::Facet_iterator it = facets_begin(); it!= facets_end(); it++)
	{
		m_totalSurface += it->surface();
	}

	return m_totalSurface;
}

double Mesh::getSurfaceTotal()
{
	return m_totalSurface;
}


void Mesh::compute2dSearchStructure()
{
	//if (m_2dSearchStructure) return;

	std::vector<Point_2_With_Vertex> points;
	points.reserve(size_of_vertices());

	NT_d p[2];
	for(Mesh::Vertex_iterator it = vertices_begin();
		it != vertices_end();
		it++)
	{
		if (it->vertex_begin() != NULL) {
			p[0] = it->point().x();
			p[1] = it->point().y();
			Point_2_With_Vertex pwv(2, p, p+2);
			pwv.handle = it;
			points.push_back(pwv);
		}
	}

	//m_2dSearchStructure = new KDTree(points.begin(),points.end());
}

/*void Mesh::compute2dSearchStructure()
{
	m_2dSearchStructure = new Enriched_point_set_2;

	std::ofstream f("c:/temp/compute2d.txt");

	f << "Computing 2d search structure on mesh with " << size_of_vertices() << std::endl << std::endl;

	//fill search structure
	int i=0;
	for (Vertex_iterator it = vertices_begin(); it != vertices_end(); it++) {
		Enriched_kernel::Point_2 p(it->point().x(), it->point().y());

		f << "[" << i << "] Adding " <<
			it->index() << " at point " <<
			p[0] << "," << p[1];

		Enriched_point_set_2::Vertex_handle p_vertex = m_2dSearchStructure->insert(p);
		assert(p_vertex != NULL);

		f << ", Is vertex handle null? " << (p_vertex == NULL) << std::endl;

		if (p_vertex != NULL && it->index() >= 0) p_vertex->info() = it;

		i++;
	}

	f << "=============================================" << std::endl;

	f.close();
}*/

bool Mesh::findClosestOnMesh(
	Enriched_kernel::Point_2 p,
	Enriched_Mesh::Vertex_handle& closestVertex,
	Enriched_Mesh::Facet_handle& containingFacet)
{
	/*assert(m_2dSearchStructure != NULL);

	NT_d p2[2];	p2[0] = p[0];	p2[1] = p[1];
	Point_2_With_Vertex queryPoint(2,p2,p2+2);

	std::vector<KDNeighborSearch::Point_with_distance> nearestNeighbors;

	KDNeighborSearch searcher(
		*m_2dSearchStructure, //kd-tree
		queryPoint, //point to look for
		KDEuclideanDistance(), //search type
		6 );//number of query points

	searcher.the_k_neighbors(std::back_inserter(nearestNeighbors));

	bool foundContainingFacet = false;
	for (int i=0; i<nearestNeighbors.size(); i++)
	{
		Vertex_handle v = nearestNeighbors[i].first->handle;

		//for debug reasons
		if (i==0) closestVertex = v;

		Halfedge_around_vertex_circulator c = v->vertex_begin();
		if (c != NULL ) do {
			Enriched_kernel::Triangle_2 t = createTriangle2(c);

			//if (t.orientation() == CGAL::RIGHT_TURN) t = t.opposite();
			if (t.has_on_bounded_side(p) || t.has_on_boundary(p)) {
				foundContainingFacet = true;
				containingFacet = c->facet();
				break;
			}
		} while (++c != v->vertex_begin());

		if (foundContainingFacet) break;
	}

	if (foundContainingFacet) {
		Point_3 p3(p[0],p[1],0.0);
		float vdistance = FLT_MAX;
		Halfedge_around_facet_circulator c = containingFacet->facet_begin();
		do {
			Enriched_kernel::Segment_3 s(p3, c->vertex()->point());
			if (s.squared_length() < vdistance) {
				vdistance = s.squared_length();
				closestVertex = c->vertex();
			}
		} while (++c != containingFacet->facet_begin());

		return true;

	} else {
		return false;
	}	*/

	return false;
}

/*bool Mesh::findClosestOnMesh(
																Enriched_kernel::Point_3 p,
																Enriched_Mesh::Vertex_handle& closestVertex,
																Enriched_Mesh::Facet_handle& containingFacet)
{
	assert(m_2dSearchStructure != NULL);

	Enriched_kernel::Point_2 p2(p.x(),p.y());

	std::list<Enriched_point_set_2::Vertex_handle> closest;
	m_2dSearchStructure->nearest_neighbors(p2, 3, std::back_inserter(closest));
	bool foundContainingFacet = false;
	do {
		Vertex_handle v = closest.front()->info();
		closest.pop_front();

		Halfedge_around_vertex_circulator c = v->vertex_begin();
		do {
			Enriched_kernel::Triangle_3 t = createTriangle(c);
			if (t.has_on(p)) {
				foundContainingFacet = true;
				containingFacet = c->facet();
				break;
			}
		} while (++c != v->vertex_begin());
	} while (closest.size() > 0 && !foundContainingFacet);

	if (foundContainingFacet) {
		float vdistance = FLT_MAX;
		Halfedge_around_facet_circulator c = containingFacet->facet_begin();
		do {
			Enriched_kernel::Segment_3 s(p, c->vertex()->point());
			if (s.squared_length() < vdistance) {
				vdistance = s.squared_length();
				closestVertex = c->vertex();
			}
		} while (++c != containingFacet->facet_begin());
		return true;

	} else {
		return false;
	}
}*/

Enriched_kernel::Triangle_3 Mesh::createTriangle(const Mesh::Halfedge_handle& h)
{
	if (h == NULL) {
		return Enriched_kernel::Triangle_3();
	}

	Enriched_kernel::Triangle_3 t(
		h->vertex()->point(),
		h->next()->vertex()->point(),
		h->next()->next()->vertex()->point());

	return t;
}

Enriched_kernel::Triangle_2 Mesh::createTriangle2(const Mesh::Halfedge_handle& h)
{
	Enriched_kernel::Point_2 p1(h->vertex()->point().x(),h->vertex()->point().y());
	Enriched_kernel::Point_2 p2(h->next()->vertex()->point().x(),
		h->next()->vertex()->point().y());
	Enriched_kernel::Point_2 p3(h->next()->next()->vertex()->point().x(),
		h->next()->next()->vertex()->point().y());

	Enriched_kernel::Triangle_2 t(p1,p2,p3);

	return t;
}

void Mesh::computeVertexFastindex()
{
	if (m_fastIndexMap.size()) m_fastIndexMap.clear();
	m_fastIndexMap.reserve(size_of_vertices());

	for (Vertex_iterator it = vertices_begin(); it != vertices_end(); it++)
	{
		m_fastIndexMap.push_back(it);
	}
}

void Mesh::computeFacetFastindex()
{
	//if (m_fastFacetMap) return;

	if (m_fastFacetMap.size()) m_fastFacetMap.clear();

	int facetTag = 0;
	m_fastFacetMap.reserve(size_of_facets());
	for (Facet_iterator it = facets_begin(); it != facets_end(); it++) {
		it->index(facetTag);
		m_fastFacetMap.push_back(it);
		facetTag++;
	}
}

Enriched_Mesh::Vertex_handle Mesh::findVertex(const int index)
{
	assert(m_fastIndexMap);

	if (!m_fastIndexMap.size()) computeVertexFastindex();

	if (index<0 || index>=m_fastIndexMap.size()) return Enriched_Mesh::Vertex_handle();

	return m_fastIndexMap[index];
}

Enriched_Mesh::Facet_handle Mesh::findFacet(const int index)
{
	assert(m_fastFacetMap);
	if (!m_fastFacetMap.size()) computeFacetFastindex();

	if (index<0 || index>=m_fastFacetMap.size()) return Enriched_Mesh::Facet_handle();

	return m_fastFacetMap[index];
}

void Mesh::computeBoundingBox()
{
	Enriched_kernel::FT xmin,xmax,ymin,ymax,zmin,zmax;
	Enriched_Mesh::computeBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);

	/*printf("xmax: %f\n", xmax);
	printf("ymax: %f\n", ymax);
	printf("zmax: %f\n", zmax);
	printf("xmin: %f\n", xmin);
	printf("ymin: %f\n", ymin);
	printf("zmin: %f\n", zmin);*/
}

double Mesh::diagonalLength()
{
	Enriched_kernel::FT xmin,xmax,ymin,ymax,zmin,zmax;
	Enriched_Mesh::computeBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);

	Point_3 m(xmin, ymin, zmin);
	Point_3 M(xmax, ymax, zmax);

	Vector_3 diag(m,M);
	return sqrt(diag*diag);
}

Point_3 Mesh::centerOfMass()
{
	Point_3 centerOfMass(0,0,0);
	for (Mesh::Vertex_iterator v = this->vertices_begin();
		 v != this->vertices_end();
		 v++)
	{
		centerOfMass = centerOfMass + (v->point() - CGAL::ORIGIN);
	}

	unsigned int size = this->size_of_vertices();

	m_computedCenterOfMass = Point_3(centerOfMass.x()/size, centerOfMass.y()/size, centerOfMass.z()/size);

	return m_computedCenterOfMass;
}

Point_3 Mesh::polyCenterOfMass()
{
	double c[3] = { 0.0, 0.0, 0.0 };
	double centerWeight = 0.0;

	for ( Mesh::Facet_iterator f = this->facets_begin(); f != this->facets_end(); ++f)
	{
		//Mesh::Facet_handle f = mesh->findFacet(m_facets[i]);
		Mesh::Point_3 fcenter;

		compute_facet_center(f, fcenter);
		double facetSurface = f->surface();
		c[0] += fcenter.x() * facetSurface;
		c[1] += fcenter.y() * facetSurface;
		c[2] += fcenter.z() * facetSurface;
		centerWeight += facetSurface;

	}

	c[0] /= centerWeight;
	c[1] /= centerWeight;
	c[2] /= centerWeight;
	return Mesh::Point_3(c[0], c[1], c[2]);

}

void Mesh::translate(const Vector_3& translateTo)
{
	Enriched_kernel::Aff_transformation_3 transformation(CGAL::TRANSLATION, translateTo);
	this->transform(transformation);
}

void Mesh::scale(const double& scaleBy)
{
	Enriched_kernel::Aff_transformation_3 transformation(CGAL::SCALING, scaleBy);
	this->transform(transformation);
}

void Mesh::transform(const Enriched_kernel::Aff_transformation_3 &transformation)
{
	for (Mesh::Vertex_iterator ivertex = this->vertices_begin(); ivertex != this->vertices_end(); ivertex++)
	{
		ivertex->point() = ivertex->point().transform(transformation);
	}
}

/*void Mesh::estimateCurvature()
{
	std::cerr << "Estimate curvature tensor...";
	double start = clock();
	typedef CCurvature_estimator<typename Mesh,Enriched_kernel> Estimator;
	Estimator estimator(this);
	estimator.run(Estimator::ONE_RING);
	double duration = (double)((clock()-start)/CLOCKS_PER_SEC);
	std::cerr << "Estimate curvature tensor...done " << duration << " seconds" << std::endl;
}*/

Mesh::Vertex_handle Mesh::splitEdge(Mesh::Halfedge_handle h)
{
	Point_3 p1 = h->vertex()->point();
	Point_3 p2 = h->opposite()->vertex()->point();

	Point_3 center((p1.x() + p2.x())/2,
		(p1.y() + p2.y())/2,
		(p1.z() + p2.z())/2);

	//delete edge h
	Mesh::Halfedge_handle g = this->join_facet(h);

	//create a center vertex
	g = this->create_center_vertex(g);
	g->vertex()->point() = center;

	return g->vertex();
}

void Mesh::splitEdges(std::vector<Mesh::Halfedge_handle> &halfEdges)
{
	for (std::vector<Mesh::Halfedge_handle>::iterator it = halfEdges.begin();
		 it != halfEdges.end();
		 it++)
	{
		splitEdge(*it);
	}

#ifdef _DEBUG
	//check that all facets are of size 3
	for ( Mesh::Facet_iterator f = this->facets_begin(); f != this->facets_end(); ++f){
		//all facets must be of size 3
		CGAL_assertion( f->size() ==  3);
	}
#endif

}

Mesh::Halfedge_handle Mesh::flipEdge(Mesh::Halfedge_handle h)
{
	h = this->join_facet(h);
	return this->split_facet(h->prev(), h->next());
}

void Mesh::flipEdges(std::vector<Mesh::Halfedge_handle> &halfEdges)
{
	std::cout << "flipping " << halfEdges.size() << " edges"<<std::endl;
	for (std::vector<Mesh::Halfedge_handle>::iterator it = halfEdges.begin();
		 it != halfEdges.end();
		 it++)
	{
		flipEdge(*it);
	}

#ifdef _DEBUG
	//check that all facets are of size 3
	for ( Mesh::Facet_iterator f = this->facets_begin(); f != this->facets_end(); ++f){
		//all facets must be of size 3
		CGAL_assertion( f->size() ==  3);
	}
#endif

}

void Mesh::flipLongEdges()
{
	//naive (slow) implementation
	int flipCount = 0;

	Mesh::Halfedge_handle edgeToFlip;

	do
	{
		edgeToFlip = NULL;
		double longest = 0;

		for (Mesh::Edge_iterator h = this->edges_begin();
			h != this->edges_end();
			h++)
		{
			if (h->vertex()->degree() < 3
				||
				h->opposite()->vertex()->degree() < 3)
				continue;

			double edgeLenSquared = edgeLengthSquared(h);

			if (edgeLenSquared > longest)
			{
				Vector_3 altVec = h->next()->vertex()->point() - h->opposite()->next()->vertex()->point();
				double alternativeEdgeLenSquared = altVec*altVec;

				if (edgeLenSquared > alternativeEdgeLenSquared)
				{
					longest = edgeLenSquared;
					edgeToFlip = h;
				}
			}
		}

		if (edgeToFlip != NULL)
		{
			flipEdge(edgeToFlip);
			flipCount++;
		}

	}
	while(edgeToFlip != NULL);;

	std::cerr << "flipped "<<flipCount<<" edges"<<std::endl;
}

void Mesh::flipEdgesAroundSmallTriangles()
{
	//naive (slow) implementation
	int flipCount = 0;

	Mesh::Halfedge_handle edgeToFlip;

	do
	{
		edgeToFlip = NULL;
		double best = 0;

		for (Mesh::Edge_iterator h = this->edges_begin();
			h != this->edges_end();
			h++)
		{
			if (h->vertex()->degree() < 3
				||
				h->opposite()->vertex()->degree() < 3)
				continue;

			Point_3 a = h->vertex()->point();
			Point_3 b = h->next()->vertex()->point();
			Point_3 c = h->opposite()->vertex()->point();
			Point_3 d = h->opposite()->next()->vertex()->point();

			Triangle_3 existing1(a,b,c);
			Triangle_3 existing2(a,c,d);
			double area1 = existing1.squared_area();
			double area2 = existing2.squared_area();

			double ratio1 = qMin(area1, area2) / qMax(area1,area2);

			if (ratio1 > 0.1)
				continue;

			Triangle_3 alternative1(b,d,a);
			Triangle_3 alternative2(b,d,c);

			area1 = alternative1.squared_area();
			area2 = alternative2.squared_area();

			double ratio2 = qMin(area1, area2) / qMax(area1,area2);

			if (ratio2 > ratio1 && (ratio2 * ratio1) > best)
			{
				best = ratio2 * ratio1;
				edgeToFlip = h;
			}
		}

		if (edgeToFlip != NULL)
		{
			flipEdge(edgeToFlip);
			flipCount++;
		}

	}
	while(edgeToFlip != NULL);;

	std::cerr << "flipped "<<flipCount<<" edges"<<std::endl;
}
double maxAngleInTriangle(const Point_3 &a, const Point_3 &b, const Point_3 &c)
{
	double ab = sqrt((b-a).squared_length());
	double ac = sqrt((c-a).squared_length());
	double bc = sqrt((b-c).squared_length());

	double a1 = (b-a)*(c-a) / (ab*ac);
	double a2 = (c-b)*(a-b) / (bc*ab);
	double a3 = (a-c)*(b-c) / (ac*bc);

	double angle = qMin(a1, qMin(a2,a3));
	return angle;
}

void Mesh::flipEdgesAroundBigAngledTriangles()
{
	//naive (slow) implementation
	int flipCount = 0;

	Mesh::Halfedge_handle edgeToFlip;

	do
	{
		edgeToFlip = NULL;
		double best = 0;

		for (Mesh::Edge_iterator h = this->edges_begin();
			h != this->edges_end();
			h++)
		{
			if (h->vertex()->degree() < 3
				||
				h->opposite()->vertex()->degree() < 3)
				continue;

			Point_3 a = h->vertex()->point();
			Point_3 b = h->next()->vertex()->point();
			Point_3 c = h->opposite()->vertex()->point();
			Point_3 d = h->opposite()->next()->vertex()->point();

			double mAngle1 = maxAngleInTriangle(a,b,c);
			double mAngle2 = maxAngleInTriangle(a,c,d);

			//min, because we work with the cosine of the angle
			double maxAngle = qMin(mAngle1, mAngle2);

			mAngle1 = maxAngleInTriangle(a,b,d);
			mAngle2 = maxAngleInTriangle(c,d,b);

			if ((b-c)*(d-c) < 0
				||
				(b-a)*(d-a) < 0)
				continue;

			//min, because we work with the cosine of the angle
			double maxAlternativeAngle = qMin(mAngle1, mAngle2);

			if (maxAlternativeAngle > maxAngle)
			{
				edgeToFlip = h;
				break;
			}
		}

		if (edgeToFlip != NULL)
		{
			flipEdge(edgeToFlip);
			flipCount++;
		}
	}
	while(edgeToFlip != NULL);;

	std::cerr << "flipped "<<flipCount<<" edges"<<std::endl;
}

double Mesh::edgeLengthSquared(Mesh::Halfedge_handle halfedge)
{
	Vector_3 vec = halfedge->vertex()->point() - halfedge->prev()->vertex()->point();
	return vec*vec;
}

void Mesh::triangulateFacet(Mesh::Facet_handle facet)
{
	unsigned int facetSize = facet->size();
	if (facetSize < 4)
		return;

	Mesh::Halfedge_handle splitStart = NULL, splitEnd = NULL;

	Mesh::Halfedge_handle hStart = facet->halfedge();

	Mesh::Halfedge_handle h = hStart;

	float best = -1;

	do
	{
		Mesh::Halfedge_handle g = h->next();

		do
		{
			if (
				g != h->next()
				&&
				h != g->next()
				)
			{
				Enriched_kernel::Vector_3 dVec = g->vertex()->point() - h->vertex()->point();
				float d = dVec*dVec;

				//avoid illegal triangulations
				bool illegalTriangulation = false;
				if (g == h->next()->next()
					&&
					h->next()->opposite()->facet() == g->opposite()->facet())
					illegalTriangulation = true;
				else
				if (h == g->next()->next()
					&&
					g->next()->opposite()->facet() == h->opposite()->facet())
					illegalTriangulation = true;
				else
				if (h->vertex() == g->vertex())
					illegalTriangulation = true;

				if (!illegalTriangulation && (best == -1 || best > d))
				{
					splitStart = h;
					splitEnd = g;
					best = d;
				}
			}
			g = g->next();
		}
		while (g != h);

		h = h->next();
	}
	while (h != hStart);

	if (best == -1)
	{
		//there is no way to triangulate the polygon, so add a center vertex instead

		//find average of neighbouring vertices
		float x=0,y=0,z=0;
		Mesh::Halfedge_handle g = hStart;
		do
		{
			Enriched_kernel::Point_3 p = g->vertex()->point();
			x += p.x();
			y += p.y();
			z += p.z();

			g = g->next();
		}while (g != hStart);

		Enriched_kernel::Point_3 center(x/facetSize,y/facetSize,z/facetSize);

		//re-create a center vertex at middle of the polygon
		g = this->create_center_vertex(hStart);
		g->vertex()->point() = center;
	}
	else
	{
		Mesh::Halfedge_handle splitter = this->split_facet(splitStart, splitEnd);
		Mesh::Halfedge_handle splitterOp = splitter->opposite();

		triangulateFacet(splitter->facet());
		triangulateFacet(splitterOp->facet());
	}
}

void Mesh::removeVertex(Mesh::Vertex_handle vertex)
{
	Mesh::Halfedge_handle h = vertex->halfedge();
    Mesh::Halfedge_handle g = this->erase_center_vertex(h);
	triangulateFacet(g->facet());
}

void Mesh::removeVertices(std::vector<Mesh::Vertex_handle> &vertices)
{
	std::cout << "removing " << vertices.size() << " vertices" << std::endl;
	for (unsigned int i=0; i<vertices.size(); i++)
	{
		removeVertex(vertices[i]);
	}

#ifdef _DEBUG
	//sanity checks
	//check that all facets are of size 3
	for ( Mesh::Facet_iterator f = this->facets_begin(); f != this->facets_end(); ++f){
		//all facets must be of size 3
		CGAL_assertion( f->size() ==  3);
	}

	for ( Mesh::Vertex_iterator v = this->vertices_begin(); v != this->vertices_end(); ++v){
		//all vertices must be of degree 3 or more
		CGAL_assertion( v->degree() >=  3);
	}

#endif
}

void Mesh::splitEdgesLongerThan(const double maxEdgeLength)
{
	double M = maxEdgeLength * maxEdgeLength;

	Mesh::Halfedge_handle best;
	do{
		//find the longest edge
		double bestM = M;
		best = NULL;

		for ( Mesh::Edge_iterator h = this->edges_begin(); h != this->edges_end(); ++h)
		{
			//meet facet-join preconditions
			if (h->vertex()->degree() < 3
				||
				h->opposite()->vertex()->degree() < 3)
				continue;

			double hVecLen = edgeLengthSquared(h);

			if (hVecLen > bestM)
			{
				bestM = hVecLen;
				best = h;
			}
		}

		//split the longest edge (if found) - by adding a vertex t its center
		if (best != NULL)
		{
			Enriched_kernel::Point_3 p1 = best->vertex()->point();
			Enriched_kernel::Point_3 p2 = best->opposite()->vertex()->point();

			Enriched_kernel::Point_3 center((p1.x() + p2.x())/2,
				(p1.y() + p2.y())/2,
				(p1.z() + p2.z())/2);

			//delete edge h
			Mesh::Halfedge_handle g = this->join_facet(best);

			//create a center vertex
			g = this->create_center_vertex(g);
			g->vertex()->point() = center;

		}

	} while (best != NULL);		//continue until no long edge is found


/*	double len2 = len*len;

	std::vector<Mesh::Halfedge_handle> edgesToSplit;
	do
	{
		edgesToSplit.clear();

		for (Mesh::Halfedge_iterator h = this->halfedges_begin();
			 h != this->halfedges_end();
			 h++)
		{
			//measure edge length
			Vector_3 vec(h->vertex()->point(), h->prev()->vertex()->point());
			int edgeLen = vec*vec;

			//add the edge into the list
			if (edgeLen > len2)
				edgesToSplit.push_back(h);
		}

		splitEdges(edgesToSplit);
	}
	while (edgesToSplit.size() > 0);*/
}


void Mesh::createTetrahedron(const double radius, const Point_3 center)
{
	double r = radius / 1.73205;

	Point_3 a( r, r, r );
	Point_3 b( -r, r, -r);
	Point_3 c( -r, -r, r);
	Point_3 d( r, -r, -r);
	this->make_tetrahedron(a,b,c,d);

	this->translate(center - CGAL::ORIGIN);
}

void Mesh::createCube(const double radius, const Point_3 center)
{
	createSphere(radius, 1, center);
}

void Mesh::createSphere(const double radius, const unsigned int subDivisions, const Point_3 center)
{
	createTetrahedron(radius);

	std::vector<Mesh::Facet_handle> facetsToSubDiv;
	std::vector<Mesh::Halfedge_handle> edgesToFlip;
	for (unsigned int level = 0; level < subDivisions; level++)
	{
		facetsToSubDiv.clear();
		edgesToFlip.clear();

		//must put all faces into a vector
		for (Mesh::Facet_iterator f = this->facets_begin();
			 f != this->facets_end();
			 f++)
			facetsToSubDiv.push_back(f);
		for (Mesh::Edge_iterator h = this->edges_begin();
			 h != this->edges_end();
			 h++)
			edgesToFlip.push_back(h);

		for (unsigned int i=0; i<facetsToSubDiv.size(); i++)
		{
			Mesh::Facet_handle f = facetsToSubDiv[i];

			Point_3 center;
			this->compute_facet_center(f,center);
			Vector_3 rad(CGAL::ORIGIN, center);
			double radLen = sqrt(rad*rad);
			rad = radius * rad / radLen;

			Mesh::Halfedge_handle h = f->halfedge();
			Mesh::Halfedge_handle g = this->create_center_vertex(h);

			g->vertex()->point() = CGAL::ORIGIN + rad;
		}

		flipEdges(edgesToFlip);
	}

	this->translate(center - CGAL::ORIGIN);
}

void Mesh::computeVolumeRangeOnFacets(double& minVal, double& maxVal)
{
	minVal = FLT_MAX;
	maxVal = -FLT_MAX;

	Mesh::Facet_const_iterator it = facets_begin();
	Mesh::Facet_const_iterator it_end = facets_end();
	for (;it != it_end; it++)
	{
		if (it->volumeSDF()<minVal) minVal = it->volumeSDF();
		if (it->volumeSDF()>maxVal) maxVal = it->volumeSDF();
	}
}

void Mesh::averageVolume()
{
	std::vector<double> avgVolume;
	avgVolume.reserve(size_of_vertices());

	for (Mesh::Vertex_iterator it = vertices_begin(); it != vertices_end(); it++)
	{
		Mesh::Halfedge_around_vertex_circulator c = it->vertex_begin();

		double avg = 0;
		double weights = 0.0;

		do {
			double weight = c->length();
			weights += weight;
			avg += weight * c->opposite()->vertex()->volumeSDF();
		} while (++c != it->vertex_begin());

		avg  /= weights;

		avg = (avg + it->volumeSDF()) / 2;

		avgVolume.push_back(avg);
	}

	int index = 0;
	for (Mesh::Vertex_iterator it = vertices_begin(); it != vertices_end(); it++)
	{
		it->volumeSDF(avgVolume[index]);
		index++;
	}
}

void Mesh::fillNormalizedFacetVolume()
{
	const static float logModifier = 4.0;
	const static double mulBy = 1 / logf(logModifier+1);

	double minVolume = FLT_MAX;
	double maxVolume = 0.0;
	Mesh::Facet_iterator fit = facets_begin();
	Mesh::Facet_iterator fit_end = facets_end();
	for (; fit != fit_end; fit++)
	{
		if (fit->volumeSDF() < minVolume)
			minVolume = fit->volumeSDF();
		if (fit->volumeSDF() > maxVolume)
			maxVolume = fit->volumeSDF();
	}

	//double range = (double) 1 / (maxVolume - minVolume);
	double range = (double) 1 / computeTotalVolume();

	fit = facets_begin();
	for (; fit != fit_end; fit++)
	{
		float nvolume = (fit->volumeSDF() - minVolume) * range;
		//nvolume = logf(nvolume * logModifier + 1) * mulBy;
		fit->volumeNSDF() = nvolume;
	}
}

void Mesh::fillNormalizedVertexVolume()
{
	const float logModifier = 4.0;
	const double mulBy = 1 / logf(logModifier+1);

	double minVolume = FLT_MAX;
	double maxVolume = 0.0;
	Mesh::Vertex_iterator vit = vertices_begin();
	Mesh::Vertex_iterator vit_end = vertices_end();
	for (; vit != vit_end; vit++)
	{
		if (vit->volumeSDF() < minVolume) minVolume = vit->volumeSDF();
		if (vit->volumeSDF() > maxVolume) maxVolume = vit->volumeSDF();
	}

	double range = (double) 1 / maxVolume;
	//double range = (double) 1 / (maxVolume - minVolume);
	//double range = (double) 1 / computeTotalVolume();

	vit = vertices_begin();
	for (; vit != vit_end; vit++)
	{
		double nvolume = vit->volumeSDF() * range;
		//float nvolume = (vit->volumeSDF() - minVolume) * range;
		//nvolume = logf(nvolume * logModifier + 1) * mulBy;
		vit->volumeNSDF() = nvolume;
	}
}

void Mesh::smoothVolume(const bool anisotropic, const float windowSize, const int iterations)
{
	const static float reduced_weight = 0.3;

	Mesh::Facet_iterator fit = facets_begin();
	Mesh::Facet_iterator fit_end = facets_end();
	for (; fit != fit_end; fit++)
	{
		double initialValue = fit->volumeSDF();
		double initialNValue = fit->volumeNSDF();

		for (int i=0;i<iterations;i++)
		{
			double smoothedValue = initialValue;
			double smoothedNValue = initialNValue;
			double weights = (initialValue==0.0 ? 0.0 : 1.0);

			Mesh::Halfedge_around_facet_circulator c = fit->facet_begin();
			do {
				if (c->opposite()->facet() != NULL)
				{
					double sneighborVolume = c->opposite()->facet()->volumeSDF();
					double sneighborNVolume = c->opposite()->facet()->volumeNSDF();
					if ( !anisotropic || (sneighborNVolume - initialNValue) < windowSize)
					{
						smoothedValue += sneighborVolume * reduced_weight;
						smoothedNValue += sneighborNVolume * reduced_weight;
						weights += reduced_weight;
					}
				}
			} while (++c != fit->facet_begin());

			//fit->volumeSDF(safeDiv(smoothedValue, weights));
			fit->volumeNSDF(safeDiv(smoothedNValue, weights));
			initialValue = fit->volumeSDF();
			initialNValue = fit->volumeNSDF();
		}
	}
}


Enriched_kernel::Segment_3 shortenRay(
	const Enriched_kernel::Segment_3& ray,
	const float distance)
{
	Enriched_kernel::Vector_3 v = ray.to_vector();
	v = v * 1/sqrt(v*v);
	Enriched_kernel::Segment_3 s(ray[0], ray[0] + v * distance);

	return s;
}

void Mesh::fillHoles()
{
	//find one hole:
	Mesh::Halfedge_handle holeEdge;
	do
	{
		holeEdge = NULL;
		Mesh::Halfedge_iterator h = this->halfedges_begin();
		while (h != this->halfedges_end() && holeEdge == NULL)
		{
			if (h->facet() == NULL && h->is_border())
				holeEdge = h;
			h++;
		}

		//fill the hole
		if (holeEdge != NULL)
		{
			holeEdge = this->fill_hole(holeEdge);
			std::cerr << "filling hole of size "<<holeEdge->facet()->size()<<std::endl;
			this->triangulateFacet(holeEdge->facet());
		}
	}
	while (holeEdge != NULL); //repeat until no holes found

}

// superimpose vertices
void Mesh::superimpose_vertices()
{
	::glBegin(GL_POINTS);
	for(Point_iterator pPoint = points_begin();
		pPoint !=	points_end();
		pPoint++)
		::glVertex3d(pPoint->x(),pPoint->y(),pPoint->z());
	::glEnd(); //	// end point assembly
}

// superimpose vertices
void Mesh::superimpose_spheres(double scale)
{
	/*GLUquadricObj* pQuadric = gluNewQuadric();
	::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

	for(Vertex_iterator pVertex = vertices_begin();
		pVertex !=	vertices_end();
		pVertex++)
	{
		if (m_renderingParams->m_renderModeVertices ==
			RenderingParamsManager::RENDER_VERTICES_IMPORTANCE_ONLY &&
			!pVertex->important())
			continue;

		if (!m_renderingParams->m_displayFiller && pVertex->index()==-1)
			continue;

		::glPushMatrix();
		double radius = average_edge_length_around(pVertex);
		::glTranslated(pVertex->point().x(),
			pVertex->point().y(),
			pVertex->point().z());
		setSphereColor(pVertex);
		::gluSphere(pQuadric,scale*radius*2,10,10);
		::glPopMatrix();
	}
	gluDeleteQuadric(pQuadric);*/
}

// superimpose edges
void Mesh::superimpose_edges(
	bool	skip_ordinary_edges,
	bool	skip_control_edges,
	bool voronoi_edge)
{
	/*if (m_renderingParams->m_renderModeEdges == RenderingParamsManager::RENDER_EDGES_NONE)
		::glBegin(GL_LINES);
	for(Edge_iterator h = edges_begin();
		h != edges_end();
		h++)
	{
		// ignore	this edges
		if(skip_ordinary_edges &&	!h->control_edge())
			continue;

		// ignore	control	edges
		if(skip_control_edges	&& h->control_edge())
			continue;

		if (m_renderingParams->m_renderModeEdges == RenderingParamsManager::RENDER_EDGES_ONLY_BOUNDARY &&
			!h->is_border_edge())
		{
			continue;
		}

		if (m_renderingParams->m_renderModeEdges == RenderingParamsManager::RENDER_EDGES_CLUSTER_SEPERATION &&
			((h->facet() != NULL && h->opposite()->facet() != NULL &&
			h->facet()->cluster() == h->opposite()->facet()->cluster()) ||
			(h->facet() == NULL || h->opposite()->facet() == NULL)
			))
		{
			continue;
		}

		//ignore filler edges if setting is so...
		if (!m_renderingParams->m_displayFiller &&
			(h->vertex()->index() == -1 ||
			h->opposite()->vertex()->index() == -1))
		{
			continue;
		}

		// assembly	and	draw line	segment
		if (m_renderingParams->m_renderModeEdges != RenderingParamsManager::RENDER_EDGES_NONE) {
			setEdgeColor(h);
			glBegin(GL_LINES);
		}
		const Point& p1 = h->prev()->vertex()->point();
		const Point& p2 = h->vertex()->point();
		::glVertex3d(p1[0],p1[1],p1[2]);
		::glVertex3d(p2[0],p2[1],p2[2]);

		if (m_renderingParams->m_renderModeEdges != RenderingParamsManager::RENDER_EDGES_NONE)
			glEnd();
	}

	if (m_renderingParams->m_renderModeEdges == RenderingParamsManager::RENDER_EDGES_NONE)
		::glEnd();*/
}

//draw facet
void Mesh::gl_draw_facet(Facet_handle pFacet,
				bool smooth_shading,
				bool use_normals)
{
	/*
	// one normal	per	face
	if(use_normals &&	!smooth_shading)
	{
		const Facet::Normal_3& normal = pFacet->normal();
		::glNormal3d(normal[0],normal[1],normal[2]);
	}

	bool highlightFacet = true;
	bool fillerFacet = false;

	int vcount  = 0;
	Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
	if (pHalfedge != NULL)
	do {
		vcount++;
		if (m_renderingParams->m_renderModeFacets == RenderingParamsManager::RENDER_FACETS_HIGHLIGHT_PATCH &&
			  pHalfedge->vertex()->distance() == FLT_MAX)
		{
			highlightFacet = false;
		}

		if (pHalfedge->vertex()->index() == -1)
			fillerFacet = true;
	} while(++pHalfedge	!= pFacet->facet_begin());

	if (!m_renderingParams->m_displayFiller && fillerFacet)
		return;

	bool useGray = false;
	if (m_renderingParams->m_renderModeFacets == RenderingParamsManager::RENDER_FACETS_DIFFERENCE_DISTANCE ||
		m_renderingParams->m_renderModeFacets == RenderingParamsManager::RENDER_FACETS_TRIANGLE_ERROR)
		useGray = true;

	if (m_renderingParams->m_renderModeFacets == RenderingParamsManager::RENDER_FACETS_HIGHLIGHT_PATCH) {
		if (highlightFacet) {
			glColor3f(0.7,0.7,0.5);
		} else {
			glColor3ubv(m_renderingParams->m_meshColor);
		}
	}

	pHalfedge = pFacet->facet_begin();
	if (pHalfedge != NULL) do
	{
		// one normal	per	vertex
		if(use_normals &&	smooth_shading)
		{
			const Facet::Normal_3& normal	=	pHalfedge->vertex()->normal();
			::glNormal3d(normal[0],normal[1],normal[2]);
		}

		// polygon assembly	is performed per vertex
		const Point& point	=	pHalfedge->vertex()->point();

		setFacetVertexColor(pFacet, pHalfedge, fillerFacet, vcount, useGray, painter);

		::glVertex3d(point[0],point[1],point[2]);
	}
	while(++pHalfedge	!= pFacet->facet_begin());
	*/
}

// draw	using	OpenGL commands	(display lists)
void Mesh::gl_draw(
	bool	smooth_shading,
	bool	use_normals)
{
	/*
	Facet_iterator pFacet	=	facets_begin();
	//find correct painter
	int color[3];
	FacetPainter* painter = NULL;
	for (int i=0; i<m_facetPainters.size(); i++)
	{
		if (m_facetPainters[i]->queryColor(pFacet, pFacet->halfedge(), m_renderingParams->m_renderModeFacets, color) != -999)
		{
			painter = m_facetPainters[i];
			break;
		}
	}

	// draw	polygons

	for(;pFacet	!= facets_end();pFacet++)
	{
		// begin polygon assembly
		::glBegin(GL_POLYGON);
		gl_draw_facet(pFacet,smooth_shading,use_normals, painter);
		::glEnd(); // end polygon assembly
	}
	glFlush();
	*/
}

void Mesh::setIndexVertices()
{
	int	index	=	0;
	for(Vertex_iterator	pVertex	=	vertices_begin();
		pVertex	!= vertices_end();
		pVertex++)
	{
		if (pVertex->index() >= 0)
			pVertex->tag(index++);
	}
}

void Mesh::writePatch(
	const char *pFilename,
	std::set<Halfedge_handle, HalfedgeCompare>& seamTree,
	int seamLoops,
	int saveMode /* = SAVEMODE_COMPLETE */,
	bool saveFillerVertices /* = false */,
	bool saveDebugInfo /* = false */)
{
	std::ofstream	stream(pFilename);

	stream << saveMode << std::endl;

	stream << m_radius << std::endl;

	switch(saveMode) {
	case SAVEMODE_COMPLETE:
		stream << size_of_vertices() << ' ' << size_of_facets() << std::endl;
		break;
	case SAVEMODE_MINIMAL:
		writeMinimalPatch(stream, seamTree, seamLoops, saveDebugInfo);
		stream.close();
		return;
		break;
	case SAVEMODE_NOCONNECTIVITY:
		stream << size_of_vertices() << ' ' << size_of_facets() << std::endl;
		break;
	}

	// output	vertices
	for(Vertex_iterator pVertex	=	vertices_begin();
		pVertex !=	vertices_end();
		pVertex++)
	{
		if (saveFillerVertices || pVertex->index() >= 0) {
			const Point_3 p = pVertex->point();
			stream <<	'v'	<< ' ' <<	p.x()	<< ' ' <<
			p.y()	<< ' ' <<
			pVertex->index();

			if (saveDebugInfo) {
				stream << ' ' << pVertex->distance();
			}

			stream << std::endl;
		}
	}

	// precompute	vertex indices
	this->setIndexVertices();

	// output	facets
	for(Facet_iterator pFacet	=	facets_begin(); pFacet !=	facets_end(); pFacet++)
	{
		bool fillerFacet = false;
		QString s("f");
		Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
		do {
			s += " " + QString::number(pHalfedge->vertex()->tag());
			if (!saveFillerVertices && pHalfedge->vertex()->index() == -1) {
				fillerFacet = true;
				break;
			}
		} while(++pHalfedge	!= pFacet->facet_begin());

		if (!fillerFacet) {
			stream <<	std::string(s.toAscii()) << std::endl;
		}
	}

	stream.close();
}

void Mesh::writeMinimalPatch(
	std::ofstream& stream,
	halfedgeSet_t& seamTree,
	int seamLoops,
	bool saveDebugInfo)
{
	std::map<int, float> index_distanceMap;

	if (seamLoops == 1 && seamTree.size() == 0) {
		for (Halfedge_iterator pHalfedge = border_halfedges_begin();
			pHalfedge != halfedges_end();
			pHalfedge++)
		{
			seamTree.insert(pHalfedge);
		}
	}

	//save all unique indexes
	for(Vertex_iterator pVertex	=	vertices_begin();
		pVertex !=	vertices_end();
		pVertex++)
	{
		if (pVertex->index() != -1)
			index_distanceMap.insert(std::pair<int, float>(pVertex->index(), pVertex->distance()));
	}

	std::set<int> indexesOnBorderSet;
	for (halfedgeSet_t::iterator it = seamTree.begin(); it != seamTree.end(); it++) {
		indexesOnBorderSet.insert((*it)->vertex()->index());
		indexesOnBorderSet.insert((*it)->opposite()->vertex()->index());
	}

	stream << index_distanceMap.size();
	if (seamLoops != 1)
		stream << ' ' << seamTree.size() << std::endl;
	else
		stream << " 0" << std::endl;

	for (std::map<int, float>::iterator it = index_distanceMap.begin();
		it!= index_distanceMap.end();
		it++)
	{
		stream /*<< "v "*/ << it->first;
		if (
			indexesOnBorderSet.find(it->first) != indexesOnBorderSet.end())
			stream << ' ' << it->second;
		stream << std::endl;
	}

	if (seamLoops != 1)
		for (halfedgeSet_t::iterator it = seamTree.begin(); it != seamTree.end(); it++) {
			stream /*<< "h "*/ << (*it)->vertex()->index() << ' ' << (*it)->opposite()->vertex()->index() << std::endl;
		}

	return;
}

bool Mesh::saveVolume(const char* fileName, const bool saveForClustering)
{
	Vertex_iterator vIt = vertices_begin();

	QByteArray array;
	//QBuffer buf( array );
	//buf.open( QIODevice::WriteOnly );
	QTextStream s( &array );

	QString str;
	for (Vertex_iterator it = vertices_begin(); it != vertices_end(); it++) {
		if (saveForClustering) {
			Point_3 p = it->point();
			//normalize point
			Point_3 newp(
				(p.x() - xmin()) / (xmax() - xmin()),
				(p.y() - ymin()) / (ymax() - ymin()),
				(p.z() - zmin()) / (zmax() - zmin()));
			str.sprintf("%f %f %f %f\n", newp[0], newp[1], newp[2], it->volumeSDF());
		} else {
			str.sprintf("%f\n", it->volumeSDF());
		}
		s << str;
	}
//	buf.close();
	s.flush();

	QFile f(QString::fromAscii(fileName));
	if (!f.open(QIODevice::WriteOnly))
		return false;

	f.writeBlock(array);
	f.close();

	return true;
}

bool Mesh::loadVolume(const char* fileName)
{
	if (!QFile::exists(fileName))
		return false;

	Vertex_iterator vIt = vertices_begin();
	QFile f(QString::fromAscii(fileName));
	if (f.open(QIODevice::ReadOnly))
	{
		QByteArray qba = f.readAll();
		QTextIStream s(&qba);
		while (!s.atEnd())
		{
			float volume = s.readLine().toFloat();
			if (vIt != vertices_end())
			{
				if (volume==volume)
					vIt->volumeSDF(volume);
				else
					vIt->volumeSDF(0.0);
				vIt++;
			}
			else
			{
				break;
			}
		}
		f.close();
		return true;
	}
	else
	{
		return false;
	}
}

void Mesh::reverseNormals()
{
	Enriched_Mesh::Vertex_iterator vit = vertices_begin();
	Enriched_Mesh::Vertex_iterator vit_end = vertices_end();

	for (;vit != vit_end; vit++) {
		vit->normal() = vit->normal()*(-1);
	}

	Enriched_Mesh::Facet_iterator fit = facets_begin();
	Enriched_Mesh::Facet_iterator fit_end = facets_end();

	for (;fit != fit_end; fit++) {
		fit->normal() = fit->normal()*(-1);
	}
}

Point_3 Mesh::getFacetCenter(const Enriched_Mesh::Facet_const_handle& f)
{
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	Enriched_Mesh::Halfedge_around_facet_const_circulator c = f->facet_begin();
	do {
		Point_3 p = c->vertex()->point();
		x += p.x();
		y += p.y();
		z += p.z();
	} while (++c != f->facet_begin());

	Point_3 center((double) x/3, (double) y/3, (double) z/3);

	return center;
}