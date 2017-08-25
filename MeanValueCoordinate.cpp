/** Raif Rustamov, Drew University
  *  evaluate 3D mean value interpolant
  *  inputs must be transposed!!!! remember fortran C differenc row first vs columns first
  *  follows the pseudocode in Ju et al.*/

//#include "mex.h"

#include "stdafx.h"

#include "MeanValueCoordinate.h"

MeanValueCoordinate::MeanValueCoordinate()
{
    //ctor
}

MeanValueCoordinate::~MeanValueCoordinate()
{
    //dtor
}

Point_3 MeanValueCoordinate::Point_3_plus(Point_3 p1, Point_3 p2)
{
	Point_3 p12(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());

	return p12;
}

Point_3 MeanValueCoordinate::Point_3_scalar_multi(float scalar, Point_3 p1)
{
	Point_3 p2(scalar * p1.x(), scalar * p1.y(), scalar * p1.z());

	return p2;
}

Point_3 MeanValueCoordinate::Point_3_scalar_div(float scalar, Point_3 p1)
{
	Point_3 p2(p1.x() / scalar, p1.y() / scalar, p1.z() / scalar);

	return p2;
}

float MeanValueCoordinate::distance(float pt1[3], float pt2[3])
{
	float res=0;
	for (int i=0; i<3; i++)
		res += (pt2[i]-pt1[i])*(pt2[i]-pt1[i]);
	
	return sqrt(res);
}

float MeanValueCoordinate::length(float pt[3])
{
	float res=0;
	for (int i=0; i<3; i++)
		res += pt[i]*pt[i];

	return sqrt(res);
}

float MeanValueCoordinate::sign_det(float v0[3], float v1[3], float v2[3])
{
	float det = v0[0]*v1[1]*v2[2] + v0[1]*v1[2]*v2[0] + v0[2]*v1[0]*v2[1] - v0[2]*v1[1]*v2[0] - v0[0]*v1[2]*v2[1]-v0[1]*v1[0]*v2[2];

	if(det >= 0.0) return 1.0;
	return -1.0;
}

bool MeanValueCoordinate::saveMVC(const QString& exportFileName, float* mvc, int no_inter_pts, int no_cage_pts)
{
	cout << "export file to: " << qPrintable(exportFileName) << endl;

	QFile file(exportFileName);
	if (file.open(QIODevice::WriteOnly)) {
		QTextStream ts(&file);
		
		for(int i=0; i<no_inter_pts; i++)
		{
			for (int j=0; j<no_cage_pts; j++)
				ts << mvc[j + no_cage_pts*i] << ' ';
			ts << '\n';
		} 

		file.close();
		return true;
	}

	return false;
}

void MeanValueCoordinate::preprocess_vertex(Mesh* mesh, float (*origins)[3])
{
	int ind = 0;
	Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
	Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
	for (; vit != vit_end; vit++)
	{
		origins[ind][0] = vit->point().x();
		origins[ind][1] = vit->point().y();
		origins[ind][2] = vit->point().z();
		ind++;
	}
}

vector<Point_3> MeanValueCoordinate::preprocess_vertex(Mesh* mesh)
{
	int no_pts = mesh->size_of_vertices();
	vector<Point_3> v_array(no_pts);
	
	int ind = 0;
	Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
	Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
	for (; vit != vit_end; vit++)
	{
		v_array[ind] = vit->point();
		ind++;
	}

	return v_array;
}

void MeanValueCoordinate::preprocess_facet(Mesh *mesh, int (*origins)[3])
{
	int ind = 0;
	Mesh::Facet_iterator fit = mesh->facets_begin();
	Mesh::Facet_iterator fit_end = mesh->facets_end();
	for (; fit != fit_end; fit++) {

		CGAL_assertion(fit->size() == 3);

		int j = 0;
		Mesh::Halfedge_around_facet_circulator fvait = fit->halfedge()->facet_begin();
		Mesh::Halfedge_around_facet_circulator begin = fvait;
		CGAL_For_all(fvait, begin)
		{
			origins[ind][j] = fvait->vertex()->index();

			j++;
		}

		ind++;
	}
}

vector<vector<int> > MeanValueCoordinate::preprocess_facet(Mesh *mesh)
{
	int no_trgs = mesh->size_of_facets();
	vector<vector<int> > facet_vertex_index(no_trgs, vector<int>(3,0));

	int ind = 0;
	Mesh::Facet_iterator fit = mesh->facets_begin();
	Mesh::Facet_iterator fit_end = mesh->facets_end();
	for (; fit != fit_end; fit++) {
		int j = 0;
		Mesh::Halfedge_around_facet_circulator fvait = fit->halfedge()->facet_begin();
		Mesh::Halfedge_around_facet_circulator begin = fvait;
		CGAL_For_all(fvait, begin)
		{
			facet_vertex_index[ind][j] = fvait->vertex()->index();
			j++;
		}

		ind++;
	}

	return facet_vertex_index;
}

float* MeanValueCoordinate::mvc(
	float x[3], int point_index,
	/*float mesh_funcs[],  int nofuncs, */
	float mesh_coord[][3], int nopts, 
	int mesh_triang[][3], int notrg) 
{
	/// allocate the variables used in the function
	//float *output = new float[nofuncs];
	//float *totalF = new float[nofuncs];
	
	// totalW = 0
	float* mvc = new float[nopts];
	for (int i = 0; i < nopts; i++)
		mvc[i] = 0;

	// distance from vertices to the x
	float *d = new float[nopts];

	// unit vector from mesh x to a cage mesh vertex
	float (*u)[3] =  new float[nopts][3];

	const float tresh1 = 0.0001;
	const float tresh2 = 0.001;
	const float tresh3 = 0.0001;

	int i, t, f, k;
	//float totalW;

	float dist;
	float l[3], theta[3], /*w[3],*/ s[3], c[3];
	int vert[3];
	float h, sss;

	//printf("vert=%g %g %g\n", x[0], x[1], x[2]);

	for(i = 0; i < nopts; i++)
	{
		dist = distance(x, mesh_coord[i]);

		/*if (dist<tresh1)
		{
			//very close to a vertex, simply return the value at the vertex
			for (f=0; f<nofuncs; f++)
			{
				output[f] = mesh_funcs[f+ nofuncs*i];
			}
			return output;
		}*/

		if (dist < tresh1)
		{
			printf("threshold1!!! mesh point index: %d, cage mesh point index: %d\n", point_index, i);

			 for (f = 0; f < nopts; f++)
			 {
				 if (f==i)
					mvc[f] = 1;
				 else
					mvc[f] = 0;
			 }

			 return mvc;
		}

		u[i][0] = (mesh_coord[i][0] - x[0]) / dist;
		u[i][1] = (mesh_coord[i][1] - x[1]) / dist;
		u[i][2] = (mesh_coord[i][2] - x[2]) / dist;

		d[i] = dist;
	}


	/*for (f=0; f<nofuncs; f++)
	{
		totalF[f] = 0;
	}
	totalW = 0;*/

	//printf("notrg = %d\n", notrg);

	for(t = 0; t < notrg; t++)
	{
		//printf("t=%d\n", t);

		for (k = 0; k < 3; k++)
		{
			//vert[k] = (int) mesh_triang[t][k] - 1; //triangle vertex indices
			vert[k] = mesh_triang[t][k]; //triangle vertex indices
			if (vert[k] < 0)
				printf("Error!!! Array index out of range\n");
		}

		//printf("vert=%d %d %d\n", vert[0], vert[1], vert[2]);

		l[0] = distance(u[vert[1]], u[vert[2]]);
		l[1] = distance(u[vert[2]], u[vert[0]]);
		l[2] = distance(u[vert[0]], u[vert[1]]);

		//printf("l are %g %g %g \n", l[0], l[1], l[2]);

		for (k = 0; k < 3; k++)
			theta[k] = 2*asin(l[k]/2);

		h = (theta[0] + theta[1] + theta[2]) / 2;

		/*if((M_PI - h) < tresh2)
		{
			printf("threshold2!!! point index:%d , triangle index: %d\n", point_index, t);

			//x lies within the triangle, use 2D MVcoords

			//w[0] = sin(theta[0])*d[vert[1]]*d[vert[2]];
			//w[1] = sin(theta[1])*d[vert[0]]*d[vert[2]];
			//w[2] = sin(theta[2])*d[vert[1]]*d[vert[0]];

			for (i = 0; i < nopts; i++)
				mvc[i] = 0;

			mvc[vert[0]] = sin(theta[0])*d[vert[1]]*d[vert[2]];
			mvc[vert[1]] = sin(theta[1])*d[vert[0]]*d[vert[2]];
			mvc[vert[2]] = sin(theta[2])*d[vert[1]]*d[vert[0]];

			//for (f = 0; f<nofuncs; f++)
				//output[f] = (w[0]*mesh_funcs[f+ nofuncs*vert[0]] + w[1]*mesh_funcs[f+ nofuncs*vert[1]]+w[2]*mesh_funcs[f+ nofuncs*vert[2]])/(w[0]+w[1]+w[2]);
			//}
			//return output;

			return mvc;
		}*/

		c[0] =2*sin(h)*sin(h-theta[0]) / (sin(theta[1])*sin(theta[2])) - 1;
		c[1] =2*sin(h)*sin(h-theta[1]) / (sin(theta[2])*sin(theta[0])) - 1;
		c[2] =2*sin(h)*sin(h-theta[2]) / (sin(theta[0])*sin(theta[1])) - 1;

		sss = sign_det(u[vert[0]], u[vert[1]], u[vert[2]]);

		for(k = 0; k < 3; k++)
			s[k] = sss * sqrt(1 - c[k] * c[k]);

		//printf("s are %g %g %g \n", s[0], s[1], s[2]);

		if( (fabs(s[0]) > tresh3) && (fabs(s[1]) > tresh3) && (fabs(s[2]) > tresh3) )
		{
			//if any is less thatn tresh then no contribution
			//belongs to the same plane but outside
			/*w[0] = (theta[0] - c[1]*theta[2] - c[2]*theta[1])/(d[vert[0]]*sin(theta[1])*s[2]);
			w[1] = (theta[1] - c[0]*theta[2] - c[2]*theta[0])/(d[vert[1]]*sin(theta[2])*s[0]);
			w[2] = (theta[2] - c[1]*theta[0] - c[0]*theta[1])/(d[vert[2]]*sin(theta[0])*s[1]);*/

			mvc[vert[0]] += (theta[0] - c[1]*theta[2] - c[2]*theta[1])/(d[vert[0]]*sin(theta[1])*s[2]);
			mvc[vert[1]] += (theta[1] - c[0]*theta[2] - c[2]*theta[0])/(d[vert[1]]*sin(theta[2])*s[0]);
			mvc[vert[2]] += (theta[2] - c[1]*theta[0] - c[0]*theta[1])/(d[vert[2]]*sin(theta[0])*s[1]);

			//totalW += w[0]+w[1]+w[2];

			//printf("totalW is %g  \n", totalW);

			/*for(f=0; f<nofuncs; f++)
			{
				totalF[f] += w[0]*mesh_funcs[f+ nofuncs*vert[0]] + w[1]*mesh_funcs[f+ nofuncs*vert[1]]+w[2]*mesh_funcs[f+ nofuncs*vert[2]];
			}*/
		}
	}

	//printf("totalW is %g  \n", totalW);

	/*for(f=0; f<nofuncs; f++)
		output[f] = totalF[f]/totalW;*/

	//delete [] output;
	//delete [] totalF;
	delete [] d;
	delete [] u;

	//return output;

	return mvc;
}

float* MeanValueCoordinate::mvc_all(
	float (*interior_points)[3], int no_inter_pts, 
	float (*mesh_coord)[3], int nopts, 
	int (*mesh_triang)[3], int notrg)
{
	// allocate memory
	float *result = new float[no_inter_pts * nopts];
	float *mvc_out;
	int point_index;


	// now for each interior point will evaluate the MVC interpolant, and put into result
	int i, j;
	for(i = 0; i < no_inter_pts; i++)
	{
		point_index = i;

		mvc_out = mvc(interior_points[i], point_index, /*mesh_funcs, nofuncs,*/ mesh_coord, nopts, mesh_triang, notrg);

		for (j=0; j<nopts; j++)
			result[j + nopts*i] = mvc_out[j];

		delete [] mvc_out;
	}

	return result;
}


vector<Point_3> MeanValueCoordinate::postprocess(int no_inter_pts, int no_cage_pts, float *mvc_coord, vector<Point_3> deform_cage_point)
{
	vector<Point_3> deform_coord(no_inter_pts);

	for(int i=0; i<no_inter_pts; i++)
	{
		float new_x = 0, new_y = 0, new_z = 0;
		float sum_weight = 0;

		for (int j=0; j<no_cage_pts; j++)
		{
			float weight = mvc_coord[j + no_cage_pts*i];

			Point_3 p1 = Point_3_scalar_multi(weight, deform_cage_point[j]);

			new_x += p1.x();
			new_y += p1.y();
			new_z += p1.z();

			sum_weight += weight;
		}

		Point_3 p(new_x, new_y, new_z);

		p = Point_3_scalar_div(sum_weight, p);

		deform_coord[i] = p;
	}

	return deform_coord;
}

bool MeanValueCoordinate::go(AppObject* appObject, Mesh* origin_mesh, Mesh* cage_mesh, Mesh* deform_cage_mesh)
{
	//PreProcessingDeformMesh();

	/// convert origin mesh vertices to interior_points array
	int no_inter_pts =  origin_mesh->size_of_vertices();
	float (*interior_points)[3] =  new float[no_inter_pts][3];
	preprocess_vertex(origin_mesh, interior_points);

	/*vector<Point_3> interior_points = preprocess_vertex(origin_mesh);*/

	printf("preprocess_vertex(origin_mesh, interior_points);\n");

	/// convert cage mesh vertices to cage_points array
	int no_cage_pts = cage_mesh->size_of_vertices();
	float (*cage_points)[3] = new float[no_cage_pts][3];
	preprocess_vertex(cage_mesh, cage_points);

	/*vector<Point_3> cage_points = preprocess_vertex(cage_mesh);*/

	printf("preprocess_vertex(cage_mesh, cage_points);\n");

	/// convert cage mesh facets to cage_facets array
	int no_cage_trg = cage_mesh->size_of_facets();
	int (*cage_facets)[3] = new int[no_cage_trg][3];
	preprocess_facet(cage_mesh, cage_facets);

	/*vector<vector<int> > no_cage_trg = preprocess_facet(cage_mesh);*/
	
	printf("preprocess_facet(cage_mesh, cage_facets);\n");

	/// create value f array
	int no_funcs = no_cage_pts;
	//int no_funcs = cage_points.size();
	//float *mesh_funcs = (float*) malloc(no_funcs);
	/*float *mesh_funcs = new float[no_funcs];
	for (int i = 0; i < no_funcs; i++)
		mesh_funcs[i] = 1;

	vector<float> mesh_funcs(no_funcs);
	for (int i = 0; i < no_funcs; i++)
		mesh_funcs[i] = 1;

	printf("create value array f done!\n");*/

	/// calculate local coordinate accordinate to cage mesh

	printf("calculate mvc\n");

	float *mvc_coord = mvc_all(interior_points, no_inter_pts, /*mesh_funcs, no_funcs,*/ cage_points, no_cage_pts, cage_facets, no_cage_trg);

	printf("mvc_all done!\n");

	/// save the mvc
	QString mvcTxtFileName = appObject->fileDir + "\\" + appObject->fileName + ".mvc.txt"; // load cage mesh
	saveMVC(mvcTxtFileName, mvc_coord, no_inter_pts, no_cage_pts);

	printf("save mvc result\n");

	delete [] interior_points;
	delete [] cage_points;
	delete [] cage_facets;
	
	//delete [] mesh_funcs;
	
	//printf("release memory done!\n");

	/// convert deform cage mesh vertices to deform_cage_points array
	int no_deform_cage_pts = deform_cage_mesh->size_of_vertices();
	//float (*deform_cage_points)[3] = new float[no_deform_cage_pts][3];
	//preprocess_vertex(deform_cage_mesh, deform_cage_points);

	vector<Point_3> deform_cage_points = preprocess_vertex(deform_cage_mesh);

	printf("preprocess_vertex(deform_cage_mesh, deform_cage_points);\n");

	/// calculate the deformed vertices position

	//float *deform_result = postprocess(no_inter_pts, no_deform_cage_pts, mvc_coord, deform_cage_points);
	vector<Point_3> deform_result = postprocess(no_inter_pts, no_deform_cage_pts, mvc_coord, deform_cage_points);

	printf("postprocess(no_inter_pts, no_deform_cage_pts, mvc_coord, deform_cage_points);\n");

	//delete [] deform_cage_points;
	delete [] mvc_coord;
	
	printf("release memory done!\n");

	int ind = 0;
	for (Mesh::Vertex_iterator ivertex = origin_mesh->vertices_begin(); ivertex != origin_mesh->vertices_end(); ivertex++)
	{
		Mesh::Vertex_handle v = ivertex;
		v->point() = deform_result[ind];
		ind++;
	}

	//delete [] deform_result;

	return true;
}