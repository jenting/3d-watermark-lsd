#pragma once

//Define CGAL classes

typedef CGAL::Simple_cartesian<double> Enriched_kernel;
typedef CGAL::Polyhedron_3<Enriched_kernel, Enriched_items> Surface;
typedef CGAL::Bbox_3 Bbox_3;

typedef Enriched_kernel::Vector_3 Vector_3;
typedef Enriched_kernel::Point_3 Point_3;
typedef Enriched_kernel::Triangle_3 Triangle_3;
typedef Enriched_kernel::Segment_3 Segment_3;
typedef Enriched_kernel::Plane_3 Plane_3;
typedef Enriched_kernel::Ray_3 Ray_3;
typedef Enriched_kernel::Line_3 Line_3;
typedef Enriched_kernel::Aff_transformation_3 Aff_transformation_3;
typedef Enriched_kernel::Sphere_3 Sphere_3;