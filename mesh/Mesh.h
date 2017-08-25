#pragma once

#include "RayIntersect.h"

//#include <CGAL/Kd_tree.h>
//#include <CGAL/Orthogonal_standard_search.h>

#include <hash_map>
#include <queue>

//Enriched mesh definition
typedef Enriched_polyhedron<Enriched_kernel,Enriched_items> Enriched_Mesh;

/** \brief Class representing meshes
 * 
 *	The Mesh class extends the already rich Polyhedron class "enriched_polyhedron" by Pierre Alliez.
 *	The class provides additional mesh operations.
 */

struct HalfedgeCompare;

//define types for 2d searchstructure
/*typedef CGAL::Triangulation_vertex_base_with_info_2<Enriched_Mesh::Vertex_handle, Enriched_kernel> Vb;
typedef CGAL::Triangulation_face_base_2<Enriched_kernel> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Point_set_2<Enriched_kernel, Tds> Enriched_point_set_2;*/

//Some definition for the KD-Tree (points including curvature and normals)
typedef CGAL::Cartesian_d<double> R;
typedef R::Point_d Point_d;
typedef Point_d::R::FT NT_d;

/**
* Derives from point_D class, used for KD-Tree searches                                                                     
*/
class Point_2_With_Vertex : public Point_d
{
public:
	Enriched_Mesh::Vertex_handle handle;	
	const int dimension() {return 2;}
	Point_2_With_Vertex()	{	}
	Point_2_With_Vertex(int d, NT_d* first, NT_d* last)
		: Point_d(d,first,last)	{	}
};

/**
 *	Class which holds information as to what parts of the mesh are selected
 */
class MeshSelection 
{
public:
	MeshSelection() {}
private:
	/** true if whole mesh is selected */
	bool m_meshSelected;
	/** vector of selected facets */
	TIndexList m_selectedFacets;
	/** vector of selected vertices */
	TIndexList m_selectedVertices;

	TIndexList m_selectedParts; // index of part

protected:

	/** clear selection and set index as selected */
	void select(const int& index, TIndexList& selectedVector) {
		selectedVector.clear();
		selectedVector.push_back(index);
	}

	void selectAdd(const int& index, TIndexList& selectedVector) { selectedVector.push_back(index); }

	/** add index as selected */
	void addToSelection(const int& index, TIndexList& selectedVector, const bool& checkForDuplicates) {
		if (checkForDuplicates) {
			for (unsigned int i=0;i<selectedVector.size();i++)
				if (selectedVector[i] == index) return;
		}
		selectedVector.push_back(index);
	}

	/** add several indices to selection */
	void addToSelection(const TIndexList& indices, TIndexList& selectedVector) { selectedVector += selectedVector; }

	int selected(const TIndexList& selectedVector) const
	{
		if (selectedVector.size()>0) 
			return selectedVector.back();
		else
			return -1;
	}

public:
	/** clear selection */
	void clear() 
	{
		m_selectedFacets.clear();
		m_selectedVertices.clear();
		m_selectedParts.clear();
	}


	void selectPart(int pindex) {
		select(pindex, m_selectedParts);
	}
	void unselectPart(int pindex)
	{
		int i = m_selectedParts.indexOf(pindex);
		if (i != -1)
			m_selectedParts.remove(i);
	}
	void selectPartAdd(int pindex) {
		selectAdd(pindex, m_selectedParts);
	}

	/** select one vertex */
	void selectVertex(const int& index) {
		select(index, m_selectedVertices);
	}

	void selectVertexAdd(int index) {
		selectAdd(index, m_selectedVertices);
	}

	/** select one facet */
	void selectFacet(const int& index) {
		select(index, m_selectedFacets);
	}

	void selectFacetAdd(int index) {
		selectAdd(index, m_selectedFacets);
	}

	/** select number of vertices */
	void selectVertices(const TIndexList& indices) {
		m_selectedVertices.clear();
		addToSelection(indices, m_selectedVertices);
	}

	/** select number of facets */
	void selectFacets(const TIndexList& indices) {
		m_selectedFacets.clear();
		addToSelection(indices, m_selectedFacets);
	}

	/** all selected facets */
	const TIndexList& selectedFacets() const {
		return m_selectedFacets;
	}

	/** all selected vertices */
	const TIndexList& selectedVertices() const {
		return m_selectedVertices;
	}

	const TIndexList& selectedParts() const {
		return m_selectedParts;
	}

	int selectedPart() const { 
		return selected(m_selectedParts);
	}

	/** first selected facet */
	int selectedFacet() const {
		return selected(m_selectedFacets);
	}

	/** first selected vertex */
	int selectedVertex() const {
		return selected(m_selectedVertices);
	}

	/** is this mesh selected */
	bool isSelected() const {
		return m_meshSelected;
	}

	/** set if this mesh is selected or not */
	void setSelected(const bool& selected) {
		m_meshSelected = selected;
	}
};


class Mesh : public Enriched_Mesh
{
private:
	static const float cm_colorWhite[3];
	static const float cm_colorGray[3];
	static const float cm_colorBlack[3];

	double m_totalSurface;

	Point_3 m_computedCenterOfMass;

	//for fast 2d search structure
	//Enriched_point_set_2* m_2dSearchStructure;
	//KDTree* m_2dSearchStructure;

	//for fast index mapping
	std::vector<Enriched_Mesh::Vertex_handle> m_fastIndexMap;

	//for fast facet mapping
	std::vector<Enriched_Mesh::Facet_handle> m_fastFacetMap;

	double m_minVolume;
	double m_maxVolume;
	double m_avgVolume;
	double m_minNVolume;
	double m_maxNVolume;
	double m_avgNVolume;
	double m_minStdDevVolume;
	double m_maxStdDevVolume;
	double m_avgStdDevVolume;
	double m_minStdDevNVolume;
	double m_maxStdDevNVolume;
	double m_avgStdDevNVolume;

	/**
	* Pointer to triangle search data structure                                                                     
	*/
	//TriangleSearch * m_triangleSearch;
	//RayIntersect m_rayIntersect;

	QString meshName;

private:
	void setIndexVertices();

	void writeMinimalPatch(
		std::ofstream& stream,
		std::set<Mesh::Halfedge_handle, HalfedgeCompare>& seamTree,
		int seamLoops,
		bool saveDebugInfo = false);	

public:
	Mesh();
	~Mesh();

	void setMeshName (QString name) { meshName = name; }
	QString getMeshName() { return meshName; }

	HDS& get_hds() { return hds; } // screw you guys.


	/*int translatePart(int from, int to);*/

	/**
	 * Computes bounding box on the mesh                                                                     
	 */
	void computeBoundingBox();

	/**
	 * Prepare 2d search structure which can return facet containing specific point...                                                                     
	 */
	void compute2dSearchStructure();

	static Enriched_kernel::Triangle_3 createTriangle(const Halfedge_handle& h);
	static Enriched_kernel::Triangle_2 createTriangle2(const Halfedge_handle& h);

	void computeVertexFastindex();

	void computeFacetFastindex();

	Enriched_Mesh::Vertex_handle findVertex(const int index);

	Enriched_Mesh::Facet_handle findFacet(const int index);

	bool findClosestOnMesh(
		Enriched_kernel::Point_2 p, 
		Enriched_Mesh::Vertex_handle& closestVertex,
		Enriched_Mesh::Facet_handle& containingFacet);

	/** Returns length of diagonal on the mesh */
	double diagonalLength();

	/**Calculates mesh's center mass */
	Point_3 centerOfMass();
	Point_3 polyCenterOfMass();
	Point_3 computedCenterOfMass()  { return m_computedCenterOfMass; }	


	void makeFacesNVolume(bool smoothing, bool smoothingAnisotropic, int smoothingIterations)
	{
		fillNormalizedFacetVolume();
		if (smoothing)
			smoothVolume(smoothingAnisotropic, 0.1, smoothingIterations);	
	}

	/** Estimates min, max, gauss and mean curvature on the mesh vertices */
	void estimateCurvature();

	/**
	 * Translates the mesh by a vector
	 * @param translateTo vector by which to translate the mesh
	 */
	void translate(const Vector_3 &translateTo);

	/**
	 * Scales the mesh by a scalar
	 * @param scaleBy multiplier by which to scale the mesh vertices
	 */
	void scale(const double &scaleBy);

	/**
	 * Transforms the mesh by an affine transformation
	 * @param transformation the transformation to apply to the mesh
	 */
	void transform(const Enriched_kernel::Aff_transformation_3 &transformation);

	//shape creation

	/**
	 * Creates a simple tetrahedron shaped mesh                                                                     
	 */
	void createTetrahedron(const double radius, const Point_3 center = CGAL::ORIGIN);
	/**
	 * Creates a simple sphere shaped mesh                                                                     
	 */
	void createSphere(const double radius, const unsigned int subDivisions, const Point_3 center = CGAL::ORIGIN);
	/**
	 * Creates a simple cube shaped mesh                                                                     
	 */
	void createCube(const double radius, const Point_3 center = CGAL::ORIGIN);

	//topological manipulations
	void splitEdges(std::vector<Mesh::Halfedge_handle> &edges);
	Mesh::Vertex_handle splitEdge(Mesh::Halfedge_handle edge);
	void flipEdges(std::vector<Mesh::Halfedge_handle> &edges);
	Mesh::Halfedge_handle flipEdge(Mesh::Halfedge_handle edge);
	void removeVertices(std::vector<Mesh::Vertex_handle> &vertices);
	void removeVertex(Mesh::Vertex_handle vertex);
	void triangulateFacet(Mesh::Facet_handle facet);

	//advanced topological manipulations
	void splitEdgesLongerThan(const double len);
	void flipLongEdges();
	void flipEdgesAroundSmallTriangles();
	void flipEdgesAroundBigAngledTriangles();

	/// various calculations
	double edgeLengthSquared(Mesh::Halfedge_handle halfedge);

	static Point_3 getFacetCenter(const Enriched_Mesh::Facet_const_handle& f);

	/** Go over all triangles and compute the volume */
	double computeTotalVolume();
	/** Go over all triangles and compute the area for each one, put it in area */
	void computeTriangleArea();

	/** Go over all triangles and compute the surface for each one, put it in ratio */
	double computeTriangleSurfaces();

	/** Go over all triangles and sum up the surface */
	double computeSurfaceTotal();

	/**Retrieve total surface if already calculated */
	double getSurfaceTotal();

	//const RayIntersect& rayIntersectError() { return m_rayIntersect; }

	void averageVolume();

	void fillNormalizedFacetVolume();
	void fillNormalizedVertexVolume();

	void smoothVolume(const bool anisotropic, const float windowSize, const int iterations);

	bool saveVolume(const char* fileName, const bool saveForClustering = false);

	bool loadVolume(const char* fileName);

	double getMinVolume() const { return m_minVolume; }
	double getMaxVolume() const { return m_maxVolume; }
	double getAvgVolume() const { return m_avgVolume; }
	double getMinNVolume() const { return m_minNVolume; }
	double getMaxNVolume() const { return m_maxNVolume; }
	double getAvgNVolume() const { return m_avgNVolume; }
	double getMinStdDevVolume() const { return m_minStdDevVolume; }
	double getMaxStdDevVolume() const { return m_maxStdDevVolume; }
	double getAvgStdDevVolume() const { return m_avgStdDevVolume; }
	double getMinStdDevNVolume() const { return m_minStdDevNVolume; }
	double getMaxStdDevNVolume() const { return m_maxStdDevNVolume; }
	double getAvgStdDevNVolume() const { return m_avgStdDevNVolume; }

	void setMinVolume(const double& v) { m_minVolume = v; }
	void setMaxVolume(const double& v) { m_maxVolume = v; }
	void setAvgVolume(const double& v) { m_avgVolume = v; }
	void setMinNVolume(const double& v) { m_minNVolume = v; }
	void setMaxNVolume(const double& v) { m_maxNVolume = v; }
	void setAvgNVolume(const double& v) { m_avgNVolume = v; }
	void setMinStdDevVolume(const double& v) { m_minStdDevVolume = v; }
	void setMaxStdDevVolume(const double& v) { m_maxStdDevVolume = v; }
	void setAvgStdDevVolume(const double& v) { m_avgStdDevVolume = v; }
	void setMinStdDevNVolume(const double& v) { m_minStdDevNVolume = v; }
	void setMaxStdDevNVolume(const double& v) { m_maxStdDevNVolume = v; }
	void setAvgStdDevNVolume(const double& v) { m_avgStdDevNVolume = v; }

	void computeVolumeRangeOnFacets(double& minVal, double& maxVal);

	/**
	 * Fill holes in mesh to make it water tight                                                                     
	 */
	void fillHoles();

	
	// draw	using	OpenGL commands	(display lists)
	void gl_draw(bool	smooth_shading, bool	use_normals);

	void gl_draw_facet(Facet_handle pFacet, bool smooth_shading, bool use_normals);

	// superimpose vertices
	void superimpose_vertices();
	// superimpose edges
	void superimpose_edges(bool skip_ordinary_edges	=	true, bool skip_control_edges = false, bool voronoi_edge = false);
	
	void superimpose_spheres(double scale);

	static const int SAVEMODE_COMPLETE = 0;
	static const int SAVEMODE_NOCONNECTIVITY = 1;
	static const int SAVEMODE_MINIMAL = 2;

	void writePatch(
		const char *pFilename,
		std::set<Halfedge_handle, HalfedgeCompare>& seamTree,
		int seamLoops,
		int saveMode = SAVEMODE_COMPLETE,
		bool saveFillerVertices = false,
		bool saveDebugInfo = false);

	void reverseNormals();

};

struct HalfedgeCompare {
	bool operator()(const Mesh::Halfedge_handle& _Left, const Mesh::Halfedge_handle& _Right) const
	{	// apply operator< to operands		
		Mesh::Vertex_handle lv = _Left->vertex();
		Mesh::Vertex_handle rv = _Right->vertex();
		if (*lv != *rv)
			return (*lv < *rv);
		else
			return (*(_Left->opposite()->vertex()) < *(_Right->opposite()->vertex()));
	}
};

typedef std::set<Mesh::Halfedge_handle, HalfedgeCompare> halfedgeSet_t;

struct VertexCompare {
	bool operator()(const Mesh::Vertex_handle& _Left, const Mesh::Vertex_handle& _Right) const
	{	// apply operator< to operands
		return (*_Left < *_Right);
	}
};

typedef std::set<Mesh::Vertex_handle, VertexCompare> vertexSet_t;

struct VertexTag2Compare {
	bool operator()(const Mesh::Vertex_handle& _Left, const Mesh::Vertex_handle& _Right) const
	{	// apply operator< to operands
		return (_Left->tag2() < _Right->tag2());
	}
};

typedef std::map<Mesh::Vertex_handle, Mesh::Point_3, VertexTag2Compare> vertexPointMap_t;

typedef std::set<Mesh::Vertex_handle, VertexTag2Compare> vertexTag2Set_t;

struct FacetLess {
	bool operator()(const Mesh::Facet_handle& _Left, const Mesh::Facet_handle& _Right) const
	{	// apply operator< to operands
		return (_Left->index() < _Right->index());
	}
};

/*struct FacetDistanceLess {
	bool operator()(const Mesh::Facet_handle& _Left, const Mesh::Facet_handle& _Right) const
	{	// apply operator< to operands
		return (_Left->ratio() < _Right->ratio());
	}
};*/

typedef std::set<Mesh::Facet_handle, FacetLess> facetSet_t;
typedef std::vector<Mesh::Facet_handle> facetVector_t;

struct vertexDistanceLess {
	bool operator() (const Mesh::Vertex_handle& v1, const Mesh::Vertex_handle& v2) {
		return v1->distance() > v2->distance();
	}
};

typedef std::priority_queue<Mesh::Vertex_handle, std::vector<Mesh::Vertex_handle>, vertexDistanceLess> vertexDistanceHeap_t;

inline float safeDiv(float a, float b)
{
	if (b != 0.0)
		return a/b;
	return 0.0;
}