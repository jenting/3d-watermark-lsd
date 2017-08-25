///////////////////////////////////////////////////////////////////////////
//																																			 //
//	Class: Enriched_polyhedron                                           //
//																																			 //
///////////////////////////////////////////////////////////////////////////

#ifndef	_POLYGON_MESH_
#define	_POLYGON_MESH_

// CGAL	stuff
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <list>


// a refined facet with a normal and a tag
template <class	Refs,	class	T, class P,	class	Norm>
class	Enriched_facet : public	CGAL::HalfedgeDS_face_base<Refs, T>
{
private:
	int m_tag; // tag
	int m_index;
	Norm m_normal; // normal

	int m_cluster;
	//short m_part;
	
	double m_area; // area

	double m_volume; // volume

	double m_volume_sdf; // SDF volume
	double m_volume_nsdf; //normalized SDF volume

	double m_volume_vsi; // VSI volume
	double m_volume_nvsi; // normalized NVSI volume

	double m_minimal_dist; // minimal distance 
	double m_maximal_dist; // maximal distance
	double m_volume_minimal_dist; // minimal distance as height
	double m_volume_maximal_dist; // maximal distance as height
	double m_volume_norm_minimal_dist; // normalized minimal distance
	double m_volume_norm_maximal_dist; // normalized maximal distance

	double m_dd; // deformation degree
	double m_ndd; // normalize deformation degree

	double m_rp_x;
	double m_rp_y;
	double m_rp_z;

	double m_surface;

	bool m_ray_intersect_err;

	//double m_ratio;
	//double m_feature;
	
public:
	// life	cycle
	// no	constructors to	repeat,	since	only
	// default constructor mandatory

	//Enriched_facet() : m_ratio(1), m_cluster(-1), m_surface(0), m_index(-1), m_volume_sdf(0.0), m_volume_nsdf(0.0), m_tag(0), m_part(-1), m_dd(0.0), m_ndd(0.0) {}
	Enriched_facet() : m_cluster(-1), m_index(-1), m_volume_sdf(0.0), m_volume_nsdf(0.0), m_tag(0), m_dd(0.0), m_ndd(0.0), m_rp_x(0.0), m_rp_y(0.0), m_rp_z(0.0), m_ray_intersect_err(false) {}

	// tag
	const int&	tag()	const {	return m_tag; }
	void tag(const int& t) { m_tag	=	t; }

	// index
	int& index() { return	m_index; }
	const int&	index() const { return	m_index; }
	void index(int v) {m_index = v;}

	// normal
	typedef	Norm Normal_3;
	const Normal_3&	normal() const { return	m_normal; }
	Normal_3&	normal() { return	m_normal; }

	// cluster
	int& cluster() { return	m_cluster;	}
	const int&	cluster() const { return	m_cluster;	}
	void cluster(int v) {m_cluster = v;}

	// part
	/*short& part() { return	m_part;	}
	const short& part() const { return	m_part;	 }
	void part(short v) {m_part = v;}*/

	// area
	double& area() { return m_area; }
	const double& area() const { return m_area; }
	void area(const  double a) { m_area = a; }

	// Volume
	double& volume() { return m_volume; }
	const double& volume() const { return m_volume; }
	void volume(const double v) { m_volume = v; }

	// SDF
	double& volumeSDF() { return m_volume_sdf; }
	const double& volumeSDF() const { return	m_volume_sdf; }
	void volumeSDF(double v) {m_volume_sdf = v;}

	double& volumeNSDF() { return	m_volume_nsdf;	}
	const double& volumeNSDF() const { return m_volume_nsdf;	}
	void volumeNSDF(const double v) {m_volume_nsdf = v;}

	// VSI
	double& volumeVSI() { return m_volume_vsi; }
	const double& volumeVSI() const { return m_volume_vsi; }
	void volumeVSI(const double v)	{ m_volume_vsi = v; }
	
	double& volumeNVSI() { return m_volume_nvsi; }
	const double& volumeNVSI() const { return m_volume_nvsi; }
	void volumeNVSI(const double v) { m_volume_nvsi = v; }

	// minimal distance
	double& minDist() { return m_minimal_dist; }
	const double& minDist() const { return m_minimal_dist; }
	void minDist(const double v) { m_minimal_dist = v; }

	// maximal distance
	double& maxDist() { return m_maximal_dist; }
	const double& maxDist() const { return m_maximal_dist; }
	void maxDist(const double v) { m_maximal_dist = v; }

	// minimal distance triangle area
	double& volumeMinDist() { return m_volume_minimal_dist; }
	const double& volumeMinDist() const { return m_volume_minimal_dist; }
	void volumeMinDist(const double v) { m_volume_minimal_dist = v; }
	
	double& volumeNormMinDist() { return m_volume_norm_minimal_dist; }
	const double& volumeNormMinDist() const { return m_volume_norm_minimal_dist; }
	void volumeNormMinDist(const double v) { m_volume_norm_minimal_dist = v; }	

	// maximal distance triangle area
	double& volumeMaxDist() { return m_volume_maximal_dist; }
	const double& volumeMaxDist() const { return m_volume_maximal_dist; }
	void volumeMaxDist(const double v) { m_volume_maximal_dist = v; }

	double& volumeNormMaxDist() { return m_volume_norm_maximal_dist; }
	const double& volumeNormMaxDist() const { return m_volume_norm_maximal_dist; }
	void volumeNormMaxDist(const double v) { m_volume_norm_maximal_dist = v; }

	// deformation degree
	double& deformationDegree() { return m_dd; }
	const double& deformationDegree() const { return m_dd; }
	void deformationDegree(double dd) { m_dd = dd; }

	// normalize deformation degree
	double& normalizeDeformationDegree() { return m_ndd; }
	const double& normalizeDeformationDegree() const { return m_ndd; }
	void normalizeDeformationDegree(double ndd) { m_ndd = ndd; }

	// VSI reference point x-coordinate
	double referencePointX() { return m_rp_x; };
	const double referencePointX() const { return m_rp_x; }
	void referencePointX(double p_x) { m_rp_x = p_x; }

	// VSI reference point y-coordinate
	double referencePointY() { return m_rp_y; };
	const double referencePointY() const { return m_rp_y; }
	void referencePointY(double p_y) { m_rp_y = p_y; }

	// VSI reference point z-coordinate
	double referencePointZ() { return m_rp_z; };
	const double referencePointZ() const { return m_rp_z; }
	void referencePointZ(double p_z) { m_rp_z = p_z; }

	// surface
	double& surface() { return m_surface; }
	const double&	surface() const { return	m_surface; }
	void surface(double v) {m_surface = v;}

	// ray
	bool& rayIntersectError() { return m_ray_intersect_err; }
	const bool & rayIntersectError() const { return m_ray_intersect_err; }
	void rayIntersectError(bool r) { m_ray_intersect_err = r; }

	// ratio
	//double& ratio() { return	m_ratio; }
	//const double& ratio() const { return	m_ratio; }
	//void ratio(double v) {m_ratio = v;}

	// feature
	//double& feature() { return	m_feature; }
	//const double&	feature() const { return	m_feature; }
	//void feature(double v) {m_feature = v;}

};


// a refined vertex with a normal and a tag
template <class	Refs,	class	T, class P,	class N, class Kernel>
class	Enriched_vertex	:	public CGAL::HalfedgeDS_vertex_base<Refs,	T, P>
{
public:
	typedef	typename N Normal;
//	typedef typename CCurvature<Kernel> Curvature;

private:
	int m_tag; 
	int m_tag2;
	int m_index;
	Normal m_normal; /// normal

	int m_cluster;

	double m_volume; /// volume

	double m_volume_sdf; /// SDF volume
	double m_volume_nsdf; /// normalized SDF volume

	double m_volume_vsi; /// VSI volume
	double m_volume_nvsi; /// normalized VSI volume

	double m_minimal_dist; /// minimal distance 
	double m_maximal_dist; /// maximal distance
	double m_volume_minimal_dist; /// minimal distance as height
	double m_volume_maximal_dist; /// maximal distance as height
	double m_volume_norm_minimal_dist; /// normalized minimal distance
	double m_volume_norm_maximal_dist; /// normalized maximal distance

	double m_square_errors; /// square errors

	double m_distance; /// distance
	double m_curvature; /// curvature
	
	double m_diff_length; /// orthogonal projection difference length

	double m_Kmin; /// min curvature
	double m_Kmax; /// max curvature
	double m_Kg; /// mean curvature
	double m_Kh;
	Normal m_VKmin;
	Normal m_VKmax;

	double m_rp_x;
	double m_rp_y;
	double m_rp_z;

	//double m_centricity;
	//bool m_important;

public:

	double curvatureMean;
	double curvatureVariance;
	double curvatureCovariance;

	// life	cycle
	//Enriched_vertex() : m_important(false), m_tag(0), m_index(0), m_distance(0), m_cluster(-1), m_volume_sdf(0) {}
	Enriched_vertex() : m_tag(0), m_index(0), m_distance(0), m_cluster(-1), m_volume_sdf(0) {}

	/// repeat	mandatory	constructors ///
	Enriched_vertex(const	P& pt)  :	CGAL::HalfedgeDS_vertex_base<Refs, T,	P>(pt) {}

	bool operator<(const Enriched_vertex& v) { return index() < v.index();	}
	bool operator==(const Enriched_vertex& v) { return index() == v.index(); }
	bool operator!=(const Enriched_vertex& v) { return index() != v.index(); }

	/// normal ///
	Normal&	normal() { return	m_normal; }
	const Normal&	normal() const { return	m_normal; }

	/// tag ///
	int& tag() {	return m_tag;	}
	const int& tag() const {	return m_tag;	}
	void tag(const int& t)	{	m_tag	=	t; }

	/// index ///
	int& index() {	return m_index;	}
	const int& index() const {	return m_index;	}
	void index(const int& t)	{	m_index	=	t; }

	// tag2 ///
	int& tag2() {	return m_tag2;	}
	const int& tag2() const {	return m_tag2;	}
	void tag2(const int& t)	{	m_tag2	=	t; }

	// cluster
	int& cluster() {	return m_cluster; }
	const int& cluster() const { return m_cluster; }
	void cluster(const int& t) { m_cluster	=	t; }

	/// Volume ///
	double& volume() { return m_volume; }
	const double& volume() const { return m_volume; }
	void volume(const double& v) { m_volume = v; }

	/// SDF volume ///
	double& volumeSDF() {	return m_volume_sdf; }
	const double& volumeSDF() const {	return m_volume_sdf;	}
	void volumeSDF(const double& t)	{	m_volume_sdf	=	t; }

	double& volumeNSDF() { return m_volume_nsdf; }
	const double& volumeNSDF() const {	return m_volume_nsdf; }
	void volumeNSDF(const double& t)	{ m_volume_nsdf =	t; }

	/// VSI ///
	double& volumeVSI() { return m_volume_vsi; }
	const double& volumeVSI() const {	return m_volume_vsi; }
	void volumeVSI(const double& v)	{ m_volume_vsi = v; }

	double& volumeNVSI() { return m_volume_nvsi; }
	const double& volumeNVSI() const { return m_volume_nvsi; }
	void volumeNVSI(const double& v) { m_volume_nvsi = v; }
	
	/// minimal distance ///
	double& minDist() { return m_minimal_dist; }
	const double& minDist() const { return m_minimal_dist; }
	void minDist(const double& v) { m_minimal_dist = v; }

	/// maximal distance ///
	double& maxDist() { return m_maximal_dist; }
	const double& maxDist() const { return m_maximal_dist; }
	void maxDist(const double& v) { m_maximal_dist = v; }

	/// minimal distance triangle area ///
	double& volumeMinDist() { return m_volume_minimal_dist; }
	const double& volumeMinDist() const { return m_volume_minimal_dist; }
	void volumeMinDist(const double& v) { m_volume_minimal_dist = v; }

	double& volumeNormMinDist() { return m_volume_norm_minimal_dist; }
	const double& volumeNormMinDist() const { return m_volume_norm_minimal_dist; }
	void volumeNormMinDist(const double& v) { m_volume_norm_minimal_dist = v; }	

	/// maximal distance triangle area ///
	double& volumeMaxDist() { return m_volume_maximal_dist; }
	const double& volumeMaxDist() const { return m_volume_maximal_dist; }
	void volumeMaxDist(const double& v) { m_volume_maximal_dist = v; }
	
	double& volumeNormMaxDist() { return m_volume_norm_maximal_dist; }
	const double& volumeNormMaxDist() const { return m_volume_norm_maximal_dist; }
	void volumeNormMaxDist(const double& v) { m_volume_norm_maximal_dist = v; }

	/// square errors ///
	double& squareErrors() { return m_square_errors; }
	const double& squareErrors() const { return m_square_errors; }
	void squareErrors(const double& e) { m_square_errors = e; }

	/// vertex curvature ///
	double& curvature() { return m_curvature; }
	const double& curvature() const { return m_curvature; }
	void curvature(const double& t) { m_curvature = t; }

	/// distance ///
	double& distance() { return m_distance; }
	const double& distance() const {	return m_distance;	}
	void distance(const double& t)	{	m_distance	=	t; }

	/// degree ///
	int degree() {
		int result = 0;
		Halfedge_around_vertex_circulator c = vertex_begin();
		if (c == NULL) return 0;
		do {
			result++;
		} while (++c != vertex->end());

		return result;
	}

	/// difference length ///
	double & diff_length() { return m_diff_length; }
	const double& diff_length() const { return m_diff_length; }
	void diff_length(const double& t) { m_diff_length = t; }
	
	/// mean curvature ///
	double& Kh() {	return m_Kh;	}
	const double& Kh() const {	return m_Kh;	}
	void Kh(const double& Kh)	{	m_Kh	=	Kh; }

	double& Kg() {	return m_Kg;	}
	const double& Kg() const {	return m_Kg;	}
	void Kg(const double& kg)	{	m_Kg	=	kg; }
	double KgNrm(){return fabs(m_Kg);}

	/// min curvature ///
	double& Kmin() {	return m_Kmin;	}
	const double& Kmin() const {	return m_Kmin;	}
	void Kmin(const double& Kmin)	{	m_Kmin	=	Kmin; }
	double KminNrm(){return fabs(m_Kmin);}

	/// max curvature ///
	double& Kmax() {	return m_Kmax;	}
	const double& Kmax() const {	return m_Kmax;	}
	void Kmax(const double& Kmax)	{	m_Kmax =	Kmax; }
	
	Normal& VKmin() {	return m_VKmin;	}
	const Normal& VKmin() const {	return m_VKmin;	}
	void VKmin(const Normal& VKmin)	{	m_VKmin	=	VKmin; }

	Normal& VKmax() {	return m_VKmax;	}
	const Normal& VKmax() const {	return m_VKmax;	}
	void VKmax(const Normal& VKmax)	{	m_VKmax	=	VKmax; }

	/// VSI reference point x-coordinate ///
	double referencePointX() { return m_rp_x; };
	const double referencePointX() const { return m_rp_x; }
	void referencePointX(double p_x) { m_rp_x = p_x; }

	/// VSI reference point y-coordinate ///
	double referencePointY() { return m_rp_y; };
	const double referencePointY() const { return m_rp_y; }
	void referencePointY(double p_y) { m_rp_y = p_y; }

	/// VSI reference point z-coordinate ///
	double referencePointZ() { return m_rp_z; };
	const double referencePointZ() const { return m_rp_z; }
	void referencePointZ(double p_z) { m_rp_z = p_z; }

	//const double& logVolume(const double range) const {	return logf(m_volume_sdf*range+1) / logf(range+1); }

	// centricity
	//double& centricity() {	return m_centricity;	}
	//const double& centricity() const {	return m_centricity;	}
	//void centricity(const double& t)	{	m_centricity	=	t; }

	// important
	//bool& important() {	return m_important;	}
	//const bool& important() const { return m_important;	}
	//void important(const bool& t)	{	m_important	=	t; }

};


// a refined halfedge with a general tag and 
// a binary tag to indicate wether it belongs 
// to the control mesh or not
template <class	Refs,	class	Tprev, class Tvertex,	class	Tface, class Norm>
class	Enriched_halfedge	:	public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
private:
	int	m_tag; // tag
	double m_weight;
	bool m_control_edge; // option for edge superimposing

public:
	// life	cycle
	Enriched_halfedge() { m_control_edge = true; }

	// tag
	const int& tag() const { return m_tag; }
	int& tag() { return m_tag;	}
	void tag(const int& t) {	m_tag	=	t; }

	// weight
	const double& weight() const { return m_weight;	}
	double& weight() { return m_weight;	}
	void weight(const double& t)	{	m_weight	=	t; }

	// control edge	
	bool& control_edge()	{ return m_control_edge; }
	const bool& control_edge()	const { return m_control_edge; }
	void control_edge(const bool& flag) { m_control_edge	=	flag;	}

	double length() const {
		Vertex_const_handle v1 = this->vertex();
		Vertex_const_handle v2 = this->opposite()->vertex();
		Segment_3 seg(v1->point(), v2->point());
		return sqrt(seg.squared_length());
	}

	/*Norm::Vector_3 vector() {
		Vertex_handle v1 = this->vertex();
		Vertex_handle v2 = this->opposite()->vertex();
		Vector_3 vec(v1->point(), v2->point());
		return vec;
	}*/
};

// A redefined items class for the Polyhedron_3	
// with	a	refined	vertex class that	contains a 
// member	for	the	normal vector	and	a	refined
// facet with	a	normal vector	instead	of the 
// plane equation	(this	is an	alternative	
// solution	instead	of using 
// Polyhedron_traits_with_normals_3).

struct Enriched_items	:	public CGAL::Polyhedron_items_3
{
		// wrap	vertex
		template <class	Refs,	class	Traits>
		struct Vertex_wrapper
		{
				typedef	typename Traits::Point_3 Point;
				typedef	typename Traits::Vector_3	Normal;
				typedef	Enriched_vertex<Refs,
													CGAL::Tag_true,
													Point,
													Normal,
				                  Traits>	Vertex;
		};

		// wrap	face
		template <class	Refs,	class	Traits>
		struct Face_wrapper
		{
				typedef	typename Traits::Point_3	Point;
				typedef	typename Traits::Vector_3	Normal;
				typedef	Enriched_facet<Refs,
												 CGAL::Tag_true,
												 Point,
												 Normal> Face;
		};

		// wrap	halfedge
		template <class	Refs,	class	Traits>
		struct Halfedge_wrapper
		{
				typedef	typename Traits::Vector_3	Normal;
				typedef	Enriched_halfedge<Refs,
														CGAL::Tag_true,
														CGAL::Tag_true,
														CGAL::Tag_true,
														Normal>	Halfedge;
		};
};

//*********************************************************
template <class	kernel,	class	items>
class	Enriched_polyhedron	:	public CGAL::Polyhedron_3<kernel,items>
{
public :
	typedef	typename kernel::FT	FT;
	typedef	typename kernel::Point_3 Point;
	typedef	typename kernel::Vector_3	Vector;
	/** typedef typename kernel::Iso_cuboid_3 Iso_cuboid; */

private	:
	// bounding box
	FT m_min[3];
	FT m_max[3];

	// type
	bool m_pure_quad;
	bool m_pure_triangle;

	//courbure extrema
	double m_MinNrmMinCurvature;
	double m_MaxNrmMinCurvature;
	double m_MinNrmMaxCurvature;
	double m_MaxNrmMaxCurvature;
	double m_MinNrmRoughCurvature;
	double m_MaxNrmRoughCurvature;
	double m_MinNrmMeanCurvature;
	double m_MaxNrmMeanCurvature;
	double m_MinNrmGaussCurvature;
	double m_MaxNrmGaussCurvature;
	double m_MinNrmSalCurvature;
	double m_MaxNrmSalCurvature;


public:
	//curvature min and max
	FT m_gaussCurvatureMin, m_gaussCurvatureMax, m_gaussCurvatureMean, m_gaussCurvatureStdDeviation;
	FT m_meanCurvatureMin, m_meanCurvatureMax, m_meanCurvatureMean, m_meanCurvatureStdDeviation;

	FT m_radius;

	/*
	int m_renderModeFacets;
	int m_renderModeEdges;
	int m_renderModeVertices;
	*/
	
	//bool m_displayFiller;

	Vertex_handle m_centerVertex;

public :

	// life	cycle
	Enriched_polyhedron()	:
		m_pure_quad(false),
		m_pure_triangle(false),
		m_radius(0),				
		m_centerVertex(NULL)		
	{}

	virtual	~Enriched_polyhedron() {	}

	// type
	bool is_pure_triangle() { return m_pure_triangle; }
	bool is_pure_quad() { return m_pure_quad; }
	
	double& MinNrmMeanCurvature() {	return m_MinNrmMeanCurvature;	}
	const double& MinNrmMeanCurvature() const {	return m_MinNrmMeanCurvature;	}
	void MinNrmMeanCurvature(const double& t)	{	m_MinNrmMeanCurvature	=	t; }

	double& MaxNrmMeanCurvature() {	return m_MaxNrmMeanCurvature;	}
	const double& MaxNrmMeanCurvature() const {	return m_MaxNrmMeanCurvature;	}
	void MaxNrmMeanCurvature(const double& t)	{	m_MaxNrmMeanCurvature	=	t; }

	double& MinNrmMinCurvature() {	return m_MinNrmMinCurvature;	}
	const double& MinNrmMinCurvature() const {	return m_MinNrmMinCurvature;	}
	void MinNrmMinCurvature(const double& t)	{	m_MinNrmMinCurvature	=	t; }

	double& MaxNrmMinCurvature() {	return m_MaxNrmMinCurvature;	}
	const double& MaxNrmMinCurvature() const {	return m_MaxNrmMinCurvature;	}
	void MaxNrmMinCurvature(const double& t)	{	m_MaxNrmMinCurvature	=	t; }

	double& MinNrmMaxCurvature() {	return m_MinNrmMaxCurvature;	}
	const double& MinNrmMaxCurvature() const {	return m_MinNrmMaxCurvature;	}
	void MinNrmMaxCurvature(const double& t)	{	m_MinNrmMaxCurvature	=	t; }

	double& MaxNrmMaxCurvature() {	return m_MaxNrmMaxCurvature;	}
	const double& MaxNrmMaxCurvature() const {	return m_MaxNrmMaxCurvature;	}
	void MaxNrmMaxCurvature(const double& t)	{	m_MaxNrmMaxCurvature	=	t; }

	double& MaxNrmRoughCurvature() {	return m_MaxNrmRoughCurvature;	}
	const double& MaxNrmRoughCurvature() const {	return m_MaxNrmRoughCurvature;	}
	void MaxNrmRoughCurvature(const double& t)	{	m_MaxNrmRoughCurvature	=	t; }

	double& MinNrmRoughCurvature() {	return m_MinNrmRoughCurvature;	}
	const double& MinNrmRoughCurvature() const {	return m_MinNrmRoughCurvature;	}
	void MinNrmRoughCurvature(const double& t)	{	m_MinNrmRoughCurvature	=	t; }

	double& MaxNrmSalCurvature() {	return m_MaxNrmSalCurvature;	}
	const double& MaxNrmSalCurvature() const {	return m_MaxNrmSalCurvature;	}
	void MaxNrmSalCurvature(const double& t)	{	m_MaxNrmSalCurvature	=	t; }

	double& MinNrmSalCurvature() {	return m_MinNrmSalCurvature;	}
	const double& MinNrmSalCurvature() const {	return m_MinNrmSalCurvature;	}
	void MinNrmSalCurvature(const double& t)	{	m_MinNrmSalCurvature	=	t; }

	double& MinNrmGaussCurvature() {	return m_MinNrmGaussCurvature;	}
	const double& MinNrmGaussCurvature() const {	return m_MinNrmGaussCurvature;	}
	void MinNrmGaussCurvature(const double& t)	{	m_MinNrmGaussCurvature	=	t; }

	double& MaxNrmGaussCurvature() {	return m_MaxNrmGaussCurvature;	}
	const double& MaxNrmGaussCurvature() const {	return m_MaxNrmGaussCurvature;	}
	void MaxNrmGaussCurvature(const double& t)	{	m_MaxNrmGaussCurvature	=	t; }


	// normals (per	facet, then	per	vertex)
	void compute_normals_per_facet() { std::for_each(facets_begin(),facets_end(),Facet_normal()); }
	void compute_normals_per_vertex() { std::for_each(vertices_begin(),vertices_end(),Vertex_normal()); }
	void computeNormals()
	{
		compute_normals_per_facet();
		compute_normals_per_vertex();
	}

	// compute bounding	box
	void computeBoundingBox(FT &xmin, FT &xmax, FT &ymin, FT &ymax, FT &zmin, FT &zmax)
	{
		if(size_of_vertices()	== 0)
			return;

		Vertex_iterator	pVertex	=	vertices_begin();
		xmin = xmax = pVertex->point().x();
		ymin = ymax = pVertex->point().y();
		zmin = zmax = pVertex->point().z();
		for(;pVertex !=	vertices_end();pVertex++)
		{
			const Point& p = pVertex->point();
			xmin	=	std::min(xmin,p.x());
			ymin	=	std::min(ymin,p.y());
			zmin	=	std::min(zmin,p.z());
			xmax	=	std::max(xmax,p.x());
			ymax	=	std::max(ymax,p.y());
			zmax	=	std::max(zmax,p.z());
		}
		m_min[0] = xmin; m_max[0] = xmax;
		m_min[1] = ymin; m_max[1] = ymax;
		m_min[2] = zmin; m_max[2] = zmax;
	}

	// bounding box
	FT xmin() const { return m_min[0]; }
	FT xmax() const { return m_max[0]; }
	FT ymin() const { return m_min[1]; }
	FT ymax() const { return m_max[1]; }
	FT zmin() const { return m_min[2]; }
	FT zmax() const { return m_max[2]; }

	// copy bounding box
	void copy_bounding_box(Enriched_polyhedron<kernel,items> *pMesh)
	{
		m_min[0] = pMesh->xmin(); m_max[0] = pMesh->xmax();
		m_min[1] = pMesh->ymin(); m_max[1] = pMesh->ymax();
		m_min[2] = pMesh->zmin(); m_max[2] = pMesh->zmax();
	}

	// degree	of a face
	static unsigned int degree(Facet_handle pFace) { return CGAL::circulator_size(pFace->facet_begin()); }
	// valence of	a	vertex
	static unsigned int valence(Vertex_handle pVertex) { return CGAL::circulator_size(pVertex->vertex_begin()); }

	// check wether	a	vertex is	on a boundary	or not
	static bool	is_border(Vertex_handle	pVertex)
	{
		Halfedge_around_vertex_circulator	pHalfEdge	=	pVertex->vertex_begin();
		if(pHalfEdge ==	NULL)	// isolated	vertex
			return true;
		Halfedge_around_vertex_circulator	d	=	pHalfEdge;
		CGAL_For_all(pHalfEdge,d)
			if(pHalfEdge->is_border())
				return true;
		return false;
	}

	// get any border	halfedge attached	to a vertex
	Halfedge_handle	get_border_halfedge(Vertex_handle	pVertex)
	{
		Halfedge_around_vertex_circulator	pHalfEdge	=	pVertex->vertex_begin();
		Halfedge_around_vertex_circulator	d	=	pHalfEdge;
		CGAL_For_all(pHalfEdge,d)
			if(pHalfEdge->is_border())
				return pHalfEdge;
		return NULL;
	}

	// tag all halfedges
	void tag_halfedges(const int tag)
	{
		for(Halfedge_iterator pHalfedge = halfedges_begin();
				pHalfedge	!= halfedges_end();
				pHalfedge++)
			pHalfedge->tag(tag);
	}

	// tag all facets
	void tag_facets(const	int	tag)
	{
		for(Facet_iterator pFace	=	facets_begin();
				pFace	!= facets_end();
				pFace++)
			pFace->tag(tag);
	}

	// set index for all vertices
	void setIndexVertices()
	{
		int	index	=	0;
		for(Vertex_iterator	pVertex	=	vertices_begin();
				pVertex	!= vertices_end();
				pVertex++)
			pVertex->tag(index++);
	}

	// is pure degree ?
	bool is_pure_degree(unsigned int d)
	{
		for(Facet_iterator pFace	=	facets_begin();
				pFace	!= facets_end();
				pFace++)
			if(degree(pFace) != d)
				return false;
		return true;
	}

	// compute type
	void compute_type()
	{
		m_pure_quad = is_pure_degree(4);
		m_pure_triangle = is_pure_degree(3);
	}

	// compute facet center
	void compute_facet_center(Facet_handle pFace,
														Point& center)
	{
		Halfedge_around_facet_circulator pHalfEdge = pFace->facet_begin();
		Halfedge_around_facet_circulator end = pHalfEdge;
		Vector vec(0.0,0.0,0.0);
		int	degree = 0;
		CGAL_For_all(pHalfEdge,end)
		{
			vec	=	vec	+	(pHalfEdge->vertex()->point()-CGAL::ORIGIN);
			degree++;
		}
		center = CGAL::ORIGIN	+	(vec/(kernel::FT)degree);
	}

	// computer average curvature around a vertex
	FT average_curvature_around(Vertex_handle pVertex)
	{
		Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
		Halfedge_around_vertex_circulator end = pHalfEdge;

		FT sum = 0.0;
		int num = 0;
		CGAL_For_all(pHalfEdge, end)
		{
			Vector vec1 = pHalfEdge->facet()->normal();
			Vector vec2 = pHalfEdge->opposite()->facet()->normal();

			double length1 = sqrt(vec1.squared_length());
			double length2 = sqrt(vec2.squared_length());
			
			double angle = acos((vec1*vec2) / (length1*length2));
			sum += angle;
			num++;
		}

		return sum / (FT) num;
	}

	// compute average edge length around a vertex
	FT average_edge_length_around(Vertex_handle pVertex)
	{
		FT sum = 0.0;
		Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
		Halfedge_around_vertex_circulator end = pHalfEdge;
		Vector vec(0.0, 0.0, 0.0);
		int	degree = 0;
		CGAL_For_all(pHalfEdge, end)
		{
			Vector vec = pHalfEdge->vertex()->point()-
				           pHalfEdge->opposite()->vertex()->point();
			sum += std::sqrt(vec*vec);
			degree++;
		}
		return sum / (FT) degree;
	}

	// compute average edge length around a vertex
	FT min_edge_length_around(Vertex_handle pVertex)
	{
		FT min_edge_length = 1e38;
		Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
		Halfedge_around_vertex_circulator end = pHalfEdge;
		Vector vec(0.0,0.0,0.0);
		int	degree = 0;
		CGAL_For_all(pHalfEdge,end)
		{
			Vector vec = pHalfEdge->vertex()->point()-
				           pHalfEdge->opposite()->vertex()->point();
			FT len = std::sqrt(vec*vec);
			if(len < min_edge_length)
				min_edge_length = len;
		}
		return min_edge_length;
	}

	// draw principal direction
	void gl_draw_directions(bool kmin, bool kmax,	bool asymptotic, bool exclude_asymptotic)
	{
		if(kmin)
			for(Vertex_iterator pVertex = vertices_begin(); pVertex !=	vertices_end(); pVertex++)
			{
				::glPushMatrix();
					::glTranslated(pVertex->point().x(),
												pVertex->point().y(),
												pVertex->point().z());
					pVertex->curvature().gl_draw_kmin(exclude_asymptotic);
				::glPopMatrix();
			}
		if(kmax)
			for(Vertex_iterator pVertex = vertices_begin();
					pVertex !=	vertices_end();
					pVertex++)
			{
				::glPushMatrix();
					::glTranslated(pVertex->point().x(),
												pVertex->point().y(),
												pVertex->point().z());
				pVertex->curvature().gl_draw_kmax(exclude_asymptotic);
				::glPopMatrix();
			}
		if(asymptotic)
			for(Vertex_iterator pVertex = vertices_begin();
					pVertex !=	vertices_end();
					pVertex++)
			{
				::glPushMatrix();
					::glTranslated(pVertex->point().x(),
												pVertex->point().y(),
												pVertex->point().z());
					pVertex->curvature().gl_draw_asymptotic();
				::glPopMatrix();
			}
	}	
		
	bool edgeOnBoundary(Halfedge_handle& h) {
		if (h->vertex()->distance() != FLT_MAX &&
				h->opposite()->vertex()->distance() != FLT_MAX &&
				(h->next()->vertex()->distance() == FLT_MAX ||
				h->opposite()->next()->vertex()->distance() == FLT_MAX))
		{
			return true;
		} else {
			return false;
		}
	}	
	
	// set random index for all vertices
	void set_random_index(unsigned int seed)
	{
		int numVertex = (int)this->size_of_vertices();
		int* random_sequence_vertex = new int[numVertex];
		int i;
		for (i=0;i<numVertex;i++)
		  random_sequence_vertex[i] = i;

	  int numFacet = (int)this->size_of_facets();
	  int* random_sequence_facet = new int[numFacet];
	  for (i=0;i<numFacet;i++)
		  random_sequence_facet[i] = i;

	  srand(seed+(unsigned int)time(NULL));
	  int temp_space, index1, index2;

	  for (i=0;i<2*numVertex;i++)
	  {
		  index1 = floor(1.0*rand()/RAND_MAX*numVertex);
		  if (index1<0)
			  index1 = 0;
		  if (index1>=numVertex)
			  index1 = numVertex-1;
		  index2 = floor(1.0*rand()/RAND_MAX*numVertex);
		  if (index2<0)
			  index2 = 0;
		  if (index2>=numVertex)
			  index2 = numVertex-1;
		  temp_space = random_sequence_vertex[index1];
		  random_sequence_vertex[index1] = random_sequence_vertex[index2];
		  random_sequence_vertex[index2] = temp_space;
	  }

	  for (i=0;i<2*numFacet;i++)
	  {
		  index1 = floor(1.0*rand()/RAND_MAX*numFacet);
		  if (index1<0)
			  index1 = 0;
		  if (index1>=numFacet)
			  index1 = numFacet-1;
		  index2 = floor(1.0*rand()/RAND_MAX*numFacet);
		  if (index2<0)
			  index2 = 0;
		  if (index2>=numFacet)
			  index2 = numFacet-1;
		  temp_space = random_sequence_facet[index1];
		  random_sequence_facet[index1] = random_sequence_facet[index2];
		  random_sequence_facet[index2] = temp_space;
	  }

	  i = 0;
      for(Vertex_iterator pVertex = vertices_begin(); pVertex != vertices_end(); pVertex++)
	  {
		  pVertex->tag(random_sequence_vertex[i]);
		  i++;
	  }

	  i = 0;
      for(Facet_iterator pFacet = facets_begin(); pFacet != facets_end(); pFacet++)
	  {
		  pFacet->tag(random_sequence_facet[i]);
		  i++;
	  }

	  delete random_sequence_vertex;
	  random_sequence_vertex = 0;
	  delete random_sequence_facet;
	  random_sequence_facet = 0;
	}

	// write in off file format (OFF).
	bool write_off(const QString& exportFileName)
	{
		cout << "export file to: " << qPrintable(exportFileName) << endl;
		
		QFile file(exportFileName);
		if (!file.open(QIODevice::WriteOnly))
			return false;

		QTextStream ts(&file);

		ts << "OFF\n";

		ts.setRealNumberPrecision(10);
		
		ts << size_of_vertices() << ' ' << size_of_facets() << ' ' << size_of_halfedges() / 2 << '\n';

		for(Point_iterator pPoint	=	points_begin(); pPoint !=	points_end(); pPoint++)
			ts << pPoint->x()	<< ' ' << pPoint->y()	<< ' ' << pPoint->z() << '\n';

		this->setIndexVertices();

		// output	facets
		for(Facet_iterator pFacet	=	facets_begin(); pFacet !=	facets_end();	 pFacet++)
		{
			ts << pFacet->facet_degree();
			Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
			do {
				ts << ' ' << pHalfedge->vertex()->tag();
			} while(++pHalfedge != pFacet->facet_begin());
			ts << '\n';
		}	 

		file.close();
		return true;
	}

	// write in	obj	file format	(OBJ).
	bool write_obj(const QString& exportFileName, int incr = 1) //	1-based	by default
	{
		cout << "export file to: " << qPrintable(exportFileName) << endl;

		QFile file(exportFileName);
		if (!file.open(QIODevice::WriteOnly))
			return false;

		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		for(Point_iterator pPoint	=	points_begin(); pPoint !=	points_end(); pPoint++)
			ts << 'v' << ' ' <<	pPoint->x()	<< ' ' << pPoint->y()	<< ' ' << pPoint->z() << '\n';

		this->setIndexVertices();

		// output	facets
		for(Facet_iterator pFacet	=	facets_begin(); pFacet !=	facets_end();	 pFacet++)
		{
			ts << 'f';
			Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
			do {
				ts << ' ' << pHalfedge->vertex()->tag()+incr;
			} while(++pHalfedge != pFacet->facet_begin());
			ts << '\n';
		}	 

		file.close();
		return true;
	}

	// write in obj file format with random reordering (OBJ).
	bool writeObjRandom(const QString& exportFileName, unsigned int currentSeed,  int incr  = 1) // 1-based by default
	{
		QFile file(exportFileName);
		if (!file.open(QIODevice::WriteOnly))
			return false;

		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		// set random indices for vertices and facets
		this->set_random_index(currentSeed);

		int numVertex = (int)this->size_of_vertices();
		double ** vertexCoordinates;
		vertexCoordinates = new double * [numVertex];
		int * vertexOrder = new int[numVertex];
		Vertex_iterator pVert;
		int i = 0;
		for (pVert = vertices_begin(); pVert != vertices_end(); pVert++)
		{
			vertexCoordinates[i] = new double[3];
			vertexCoordinates[i][0] = pVert->point().x();
			vertexCoordinates[i][1] = pVert->point().y();
			vertexCoordinates[i][2] = pVert->point().z();
			vertexOrder[pVert->tag()] = i;
			i++;
		}

		// output vertices
		for (i=0;i<numVertex;i++)
		{
			int indexTempVertex = vertexOrder[i];
			 ts << 'v' << ' ' << vertexCoordinates[indexTempVertex][0] << ' ' << vertexCoordinates[indexTempVertex][1] << ' ' << vertexCoordinates[indexTempVertex][2] << '\n';
		}

		int numFacet = (int)this->size_of_facets();
		Facet_iterator pFacet;
		int * facetOrder = new int[numFacet];
		vector<int> * indicesFacet = new vector<int>[numFacet];
		i = 0;
		for(pFacet = facets_begin(); pFacet != facets_end(); pFacet++)
		{
			Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
			do
			{
				indicesFacet[i].push_back(pHalfedge->vertex()->tag());
			}
			while(++pHalfedge != pFacet->facet_begin());
			facetOrder[pFacet->tag()] = i;
			i++;
		}

		// output facets
		for (int i=0;i<numFacet;i++)
		{
			int indexTempFacet = facetOrder[i]; 
			ts << 'f';
			for (unsigned int j=0; j<indicesFacet[indexTempFacet].size(); j++)
				ts << " " << indicesFacet[indexTempFacet][j]+incr;
			ts << '\n';
		}

		file.close();

		delete [] vertexCoordinates;
		vertexCoordinates = 0;
		delete [] vertexOrder;
		vertexOrder = 0;
		delete [] facetOrder;
		facetOrder = 0;
		delete [] indicesFacet;
		indicesFacet = 0;

		return true;
	}

	void triangulate()
    {
        Facet_iterator f = this->facets_begin();
        Facet_iterator f2 = this->facets_begin();
        do //for (; f != this->facets_end(); f++)
        {
            f = f2;
            if (f == this->facets_end())
            {
                break;
            }
            f2++;

            if (!(f->is_triangle()))
            {
                int num = (int)(f->facet_degree() - 3);
                Halfedge_handle h = f->halfedge();

                h = this->make_hole(h);

                Halfedge_handle g = h->next();
                g = g->next();

                g = this->add_facet_to_border (h, g);

                num--;
                while (num != 0)
                {
                    g = g->opposite();
                    g = g->next();
                    g = this->add_facet_to_border (h, g);
                    num--;
                }

                this->fill_hole(h);
            }

        } while (true);

		this->computeNormals();
		this->compute_type();
    }

	// draw bounding box
	void gl_draw_bounding_box()
	{
		::glBegin(GL_LINES);
			// along x axis
			::glVertex3d(m_min[0],m_min[1],m_min[2]);
			::glVertex3d(m_max[0],m_min[1],m_min[2]);
			::glVertex3d(m_min[0],m_min[1],m_max[2]);
			::glVertex3d(m_max[0],m_min[1],m_max[2]);
			::glVertex3d(m_min[0],m_max[1],m_min[2]);
			::glVertex3d(m_max[0],m_max[1],m_min[2]);
			::glVertex3d(m_min[0],m_max[1],m_max[2]);
			::glVertex3d(m_max[0],m_max[1],m_max[2]);
			// along y axis
			::glVertex3d(m_min[0],m_min[1],m_min[2]);
			::glVertex3d(m_min[0],m_max[1],m_min[2]);
			::glVertex3d(m_min[0],m_min[1],m_max[2]);
			::glVertex3d(m_min[0],m_max[1],m_max[2]);
			::glVertex3d(m_max[0],m_min[1],m_min[2]);
			::glVertex3d(m_max[0],m_max[1],m_min[2]);
			::glVertex3d(m_max[0],m_min[1],m_max[2]);
			::glVertex3d(m_max[0],m_max[1],m_max[2]);
			// along z axis
			::glVertex3d(m_min[0],m_min[1],m_min[2]);
			::glVertex3d(m_min[0],m_min[1],m_max[2]);
			::glVertex3d(m_min[0],m_max[1],m_min[2]);
			::glVertex3d(m_min[0],m_max[1],m_max[2]);
			::glVertex3d(m_max[0],m_min[1],m_min[2]);
			::glVertex3d(m_max[0],m_min[1],m_max[2]);
			::glVertex3d(m_max[0],m_max[1],m_min[2]);
			::glVertex3d(m_max[0],m_max[1],m_max[2]);
		::glEnd();
	}

	// count #boundaries
	unsigned int nb_boundaries()
	{
		unsigned int nb = 0;
		tag_halfedges(0);
		Halfedge_handle	seed_halfedge	=	NULL;
		while((seed_halfedge = get_border_halfedge_tag(0)) !=	NULL)
		{
			nb++;
			seed_halfedge->tag(1);
			Vertex_handle	seed_vertex	=	seed_halfedge->prev()->vertex();
			Halfedge_handle	current_halfedge = seed_halfedge;
			Halfedge_handle	next_halfedge;
			do
			{
				next_halfedge	=	current_halfedge->next();
				next_halfedge->tag(1);
				current_halfedge = next_halfedge;
			}
			while(next_halfedge->prev()->vertex()	!= seed_vertex);
		}
		return nb;
	}

	// get any border	halfedge with	tag
	Halfedge_handle get_border_halfedge_tag(int	tag)
	{
		for(Halfedge_iterator pHalfedge	=	halfedges_begin();
				pHalfedge	!= halfedges_end();
				pHalfedge++)
			if(pHalfedge->is_border()	&&
				 pHalfedge->tag()	== tag)
				return pHalfedge;
		return NULL;
	}

	// get any facet with	tag
	Facet_handle get_facet_tag(const int tag)
	{
		for(Facet_iterator pFace	=	facets_begin();
				pFace	!= facets_end();
				pFace++)
			if(pFace->tag()	== tag)
				return pFace;
		return NULL;
	}

	// tag component 
	void tag_component(Facet_handle	pSeedFacet,
										 const int tag_free,
										 const int tag_done)
{
		pSeedFacet->tag(tag_done);
		std::list<Facet_handle> facets;
		facets.push_front(pSeedFacet);
		while(!facets.empty())
		{
			Facet_handle pFacet = facets.front();
			facets.pop_front();
			pFacet->tag(tag_done);
			Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
			Halfedge_around_facet_circulator end = pHalfedge;
			CGAL_For_all(pHalfedge,end)
			{
				Facet_handle pNFacet = pHalfedge->opposite()->facet();
				if(pNFacet !=	NULL && pNFacet->tag() == tag_free)
				{
					facets.push_front(pNFacet);
					pNFacet->tag(tag_done);
				}
			}
		}
}

	// count #components
	unsigned int nb_components()
	{
		unsigned int nb = 0;
		tag_facets(0);
		Facet_handle seed_facet	=	NULL;
		while((seed_facet	=	get_facet_tag(0))	!= NULL)
		{
			nb++;
			tag_component(seed_facet,0,1);
		}
		return nb;
	}

	// compute the genus
	// V - E + F + B = 2 (C	-	G)
	// C ->	#connected components
	// G : genus
	// B : #boundaries
	int	genus()
	{
		int	c	=	nb_components();
		int	b	=	nb_boundaries();
		int	v	=	size_of_vertices();
		int	e	=	size_of_halfedges()/2;
		int	f	=	size_of_facets();
		return genus(c,v,f,e,b);
	}
	int	genus(int	c,
						int	v,
						int	f,
						int	e,
						int	b)
	{
		return (2*c+e-b-f-v)/2;
	}

	double getMaxDim()
	{
		FT xmin,xmax,ymin,ymax,zmin,zmax;
		computeBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
		
		FT max = xmax - xmin;
		if(ymax-ymin>max)
			max = ymax-ymin;
		if(zmax-zmin>max)
			max = zmax-zmin;

		return max;
	}
	
	void KmaxKmean(double coef)
	{
		double max=0;
		double min=10000000;
		for(Vertex_iterator	pVertex	=	vertices_begin();pVertex!= vertices_end();pVertex++)
		{
			double kmax=pVertex->Kmax()*coef;
			double kmin=pVertex->Kmin()*coef;

			pVertex->Kmax((kmax+kmin)/2.);

			if(pVertex->Kmax()>max)
				max=pVertex->Kmax();
			if(pVertex->Kmax()<min)
				min=pVertex->Kmax();

		}

		MinNrmMaxCurvature(min);
		MaxNrmMaxCurvature(max);

		//GaussianNormalisationKmax();
	}

};


// compute facet normal	
struct Facet_normal	// (functor)
{
	template <class	Facet>
	void operator()(Facet& f)
	{
		typename Facet::Normal_3 sum = CGAL::NULL_VECTOR;
		typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
		if (h != NULL && h->next() != NULL && h->next()->next() != NULL) 
			do {
				assert(h == NULL) ;
				typename Facet::Normal_3 normal = CGAL::cross_product(
					h->next()->vertex()->point() - h->vertex()->point(),
					h->next()->next()->vertex()->point() - h->next()->vertex()->point());
				double sqnorm = normal * normal;
				if(sqnorm	!= 0)
					normal = normal	/	(double)std::sqrt(sqnorm);
				sum	=	sum	+	normal;
			} while(++h	!= f.facet_begin());

		double sqnorm = sum * sum;
		if(sqnorm	!= 0.0)
			f.normal() = sum / std::sqrt(sqnorm);
		else
		{
			f.normal() = CGAL::NULL_VECTOR;
			//TRACE("degenerate face\n");
			//std::cerr << "degenerate face" << std::endl;
		}
	}
};


// compute vertex	normal 
struct Vertex_normal //	(functor)
{
		template <class	Vertex>
		void operator()(Vertex&	v)
		{
				typename Vertex::Normal	normal = CGAL::NULL_VECTOR;
				Vertex::Halfedge_around_vertex_const_circulator	pHalfedge	=	v.vertex_begin();
				Vertex::Halfedge_around_vertex_const_circulator	begin	=	pHalfedge;
				CGAL_For_all(pHalfedge,begin)	
					if(!pHalfedge->is_border())
						normal = normal	+	pHalfedge->facet()->normal();
				double	sqnorm = normal * normal;
				if(sqnorm != 0.0f)
					v.normal() = normal	/	(double)std::sqrt(sqnorm);
				else
					v.normal() = CGAL::NULL_VECTOR;
		}
};

#endif //	_POLYGON_MESH_
