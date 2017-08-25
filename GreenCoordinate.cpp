#include "stdafx.h"

#include "GreenCoordinate.h"

GreenCoordinate::GreenCoordinate()
{
}

GreenCoordinate::~GreenCoordinate()
{
}

Point_3 GreenCoordinate::Point_3_plus(Point_3 p1, Point_3 p2)
{
	Point_3 p12(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());

	return p12;
}

Point_3 GreenCoordinate::Point_3_minus(Point_3 p1, Point_3 p2)
{
	Point_3 p12(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());

	return p12;
}

Point_3 GreenCoordinate::Point_3_scalar_multi(float scalar, Point_3 p1)
{
	Point_3 p2(scalar * p1.x(), scalar * p1.y(), scalar * p1.z());

	return p2;
}

Point_3 GreenCoordinate::Point_3_scalar_div(float scalar, Vector_3 p1)
{
	Point_3 p2(p1.x() / scalar, p1.y() / scalar, p1.z() / scalar);

	return p2;
}

bool GreenCoordinate::saveGC(
	const QString& vertexGCTxtFileName, 
	const QString& facetGCTxtFileName, 
	float* gc_coord, 
	int no_inter_pts, int no_cage_pts, int no_cage_trg)
{
	cout << "export file to: " << qPrintable(vertexGCTxtFileName) << endl;

	QFile file_v(vertexGCTxtFileName);
	if (file_v.open(QIODevice::WriteOnly)) {
		Q3TextStream ts_v(&file_v);
		
		for(int i=0; i<no_inter_pts; i++)
		{
			for (int j=0; j<no_cage_pts; j++)
				ts_v << gc_coord[j + no_cage_pts*i] << ' ';
			ts_v << '\n';
		} 

		file_v.close();
	}
	else {
		return false;
	}

	cout << "export file to: " << qPrintable(facetGCTxtFileName) << endl;

	QFile file_f(facetGCTxtFileName);
	if (file_f.open(QIODevice::WriteOnly)) {
		Q3TextStream ts_f(&file_f);
		
		int size_result_points = no_inter_pts * no_cage_pts;
		for(int i=0; i<no_inter_pts; i++)
		{
			for (int j=0; j<no_cage_trg; j++)
				ts_f << gc_coord[size_result_points + j + no_cage_trg*i] << ' ';
			ts_f << '\n';
		} 

		file_f.close();

		return true;
	}
	
	return false;
}

bool GreenCoordinate::saveVertexGC(const QString& exportFileName, float* gc_v, int no_inter_pts, int no_cage_pts)
{
	cout << "export file to: " << qPrintable(exportFileName) << endl;

	QFile file(exportFileName);
	if (file.open(QIODevice::WriteOnly)) {
		QTextStream ts(&file);
		
		for(int i=0; i<no_inter_pts; i++)
		{
			for (int j=0; j<no_cage_pts; j++)
				ts << gc_v[j + no_cage_pts*i] << ' ';
			ts << '\n';
		} 

		file.close();
		return true;
	}

	return false;
}
bool GreenCoordinate::saveFacetGC(const QString& exportFileName, float* gc_f, int no_inter_pts, int no_cage_trg)
{
	cout << "export file to: " << qPrintable(exportFileName) << endl;

	QFile file(exportFileName);
	if (file.open(QIODevice::WriteOnly)) {
		QTextStream ts(&file);
		
		for(int i=0; i<no_inter_pts; i++)
		{
			for (int j=0; j<no_cage_trg; j++)
				ts << gc_f[j + no_cage_trg*i] << ' ';
			ts << '\n';
		} 

		file.close();
		return true;
	}

	return false;
}

vector<Vector_3> GreenCoordinate::preprocess_facet_norm(Mesh* mesh)
{
	int no_trgs = mesh->size_of_facets();
	vector<Vector_3> f_n_array(no_trgs);

	int ind = 0;
	Mesh::Facet_iterator fit = mesh->facets_begin();
	Mesh::Facet_iterator fit_end = mesh->facets_end();
	for (; fit != fit_end; fit++)
	{	
		Vector_3 normal = fit->normal();
		
		f_n_array[ind] = fit->normal();

		ind++;
	}

	return f_n_array;
}

vector<Point_3> GreenCoordinate::preprocess_vertex(Mesh* mesh)
{
	int no_pts = mesh->size_of_vertices();
	vector<Point_3> p_array(no_pts);

	int ind = 0;
	Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
	Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
	for (; vit != vit_end; vit++)
	{
		p_array[ind] = vit->point();

		ind++;
	}

	return p_array;
}

vector<vector<int> > GreenCoordinate::preprocess_facet(Mesh *mesh)
{
	int no_trgs = mesh->size_of_facets();
	vector<vector<int> > f_array(no_trgs, vector<int>(3,0));

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
			//f_array[ind][j] = fvait->vertex()->index() - 1;

			f_array[ind][j] = fvait->vertex()->index();

			if (f_array[ind][j] < 0)
				printf("Error! vertex index less than 0\n");

			j++;
		}

		ind++;
	}

	return f_array;
}

vector<Point_3> GreenCoordinate::postprocess(
		float* gc_coord, int no_inter_pts, 
		vector<Point_3> deform_cage_vertices, 
		vector<Vector_3> deform_cage_trg_norm, 
		vector<float> scaling_factor)
{
	int no_cage_pts = deform_cage_vertices.size();
	int no_cage_trg = deform_cage_trg_norm.size();
	int no_pts = no_inter_pts * no_cage_pts;

	vector<Point_3> deform_coord(no_inter_pts);
	for(int i = 0; i < no_inter_pts; i++)
	{
		float deform_v_x = 0, deform_v_y = 0, deform_v_z = 0;
		float deform_f_x = 0, deform_f_y = 0, deform_f_z = 0;

		for (int j = 0; j < no_cage_pts; j++)
		{
			float gc_v = gc_coord[j + no_cage_pts * i];

			deform_v_x += gc_v * deform_cage_vertices[j].x();
			deform_v_y += gc_v * deform_cage_vertices[j].y();
			deform_v_z += gc_v * deform_cage_vertices[j].z();
		}

		Point_3 deform_v(deform_v_x, deform_v_y, deform_v_z);

		for (int j = 0; j < no_cage_trg; j++)
		{
			float gc_f = gc_coord[no_pts + j + no_cage_trg * i];
			float factor = scaling_factor[j]; 
			
			deform_f_x += gc_f  * factor * deform_cage_trg_norm[j].x();
			deform_f_y += gc_f  * factor * deform_cage_trg_norm[j].y();
			deform_f_z += gc_f  * factor * deform_cage_trg_norm[j].z();
		}

		Point_3 deform_f(deform_f_x, deform_f_y, deform_f_z);

		Point_3 deform_p(deform_v.x() + deform_f.x(), deform_v.y() + deform_f.y(), deform_v.z() + deform_f.z());

		deform_coord[i] = deform_p;
	}

	return deform_coord;
}

vector<float> GreenCoordinate::get_scaling_factor(vector<Vector_3> cage_trg_norm, vector<Vector_3> deform_cage_trg_norm)
{
	int no_trgs = cage_trg_norm.size();
	vector<float> scaling_factor(no_trgs);

	for (int i = 0; i < no_trgs; i++)
	{
		/*float trg_norm_length = sqrt(cage_trg_norm[i].squared_length());
		float deform_trg_norm_length = sqrt(deform_cage_trg_norm[i].squared_length());

		scaling_factor[i] = deform_trg_norm_length / trg_norm_length;*/

		scaling_factor[i] = 1;
	}

	return scaling_factor;
}


float GreenCoordinate::TriInt(Point_3 p, Point_3 v1, Point_3 v2, Point_3 eta)
{
		/// alpha = acos( (v2-v1) * (p-v1) / (||v2-v1|| * ||p-v1||) )
	Vector_3 v21(v2.x() - v1.x(), v2.y() - v1.y(), v2.z() - v1.z());
	Vector_3 vp1(p.x() - v1.x(), p.y() - v1.y(), p.z() - v1.z());
	float v21_length = sqrt(v21.squared_length());
	float vp1_length = sqrt(vp1.squared_length());
	float alpha = acos( v21*vp1 / v21_length*vp1_length );


	/// beta = acos( (v1-p) * (v2-p) / (||v1-p|| * ||v2-p||) )
	Vector_3 v1p(v1.x() - p.x(), v1.y() - p.y(), v1.z() - p.z());
	Vector_3 v2p(v2.x() - p.x(), v2.y() - p.y(), v2.z() - p.z());
	float v1p_length = sqrt(v1p.squared_length());
	float v2p_length = sqrt(v2p.squared_length());
	float beta = acos( v1p*v2p / v1p_length*v2p_length );

	/// lambda = || p - v1 || ^ 2 * sin(alpha) ^2
	float lambda = vp1.squared_length() * sin(alpha) * sin(alpha);


	/// c = || p - eta || ^ 2
	Vector_3 vpeta(p.x() - eta.x(), p.y() - eta.y(), p.z() - eta.z());
	float c = vpeta.squared_length();


	/// theta = pi - alpha
	float theta_pi_alpha = M_PI - alpha;
	float theta_pi_alpha_beta = M_PI - alpha - beta;

	float integer_pi_alpha = TriIntTheta(lambda, c, theta_pi_alpha);
	float integer_pi_alpha_beta = TriIntTheta(lambda, c, theta_pi_alpha_beta);

	//float r = (integer_pi_alpha - integer_pi_alpha_beta - sqrt(c) * beta * beta) / (-4) * M_PI;
	float r = fabs(integer_pi_alpha - integer_pi_alpha_beta - sqrt(c) * beta) / (-4) * M_PI;

	//printf("r = %f\n", r);

	return r;
}

float GreenCoordinate::TriIntTheta(float lambda, float c, float theta)
{
	float S = sin(theta);
	float C = cos(theta);

	int sign_S = sign_func(S);

	float atan_term = sqrt(c) * C / sqrt(lambda + S * S * c);
	float log_first_term = 2 * sqrt(lambda) * S * S / (1 - C) * (1 - C);
	float log_sec_term = 1 - 2 * c * C / (c * (1 + C) + lambda + sqrt(lambda * lambda + lambda * c * S * S));

	float  int_theta = sign_S * (-0.5) * (2 * sqrt(c) * atan(atan_term) + sqrt(lambda) * log(log_first_term * log_sec_term));

	return int_theta;
}

/*void GreenCoordinate::each_facet(
	Point_3 x, int point_index, int triangle_index, vector<int> triangle_vertices_index,
	vector<Point_3> triangle_vertices,
	Vector_3 triangle_norm,
	float* result_vertex, float* result_facet,
	int no_cage_vertices, int no_cage_facets)
{
	const float thres = 1e-5;

	vector<int> s(3);
	vector<float> I(3);
	vector<float> II(3);
	vector<Vector_3> q(3);
	vector<Vector_3> N(3);

	vector<Vector_3> triangle_vector(triangle_vertices.size());
	for (int l = 0; l < 3; l++)
		triangle_vector[l] = triangle_vertices[l] - x;

	Vector_3 p = (triangle_vector[0] * triangle_norm) * triangle_norm;

	for (int l = 0; l < 3; l++)
	{
		// sign function
		Vector_3 cross_vector = CGAL::cross_product(triangle_vector[l] - p, triangle_vector[(l+1)%3] - p);
		float value = cross_vector * triangle_norm;
		if(value >=0)
			s[l] = 1;
		else
			s[l] = -1;

		// I : GCTriInt
		Vector_3 zero(0, 0, 0);
		I[l] = TriInt(p, triangle_vector[l], triangle_vector[(l+1)%3], zero);

		// II : GCTriInt
		II[l] = TriInt(zero, triangle_vector[(l+1)%3], triangle_vector[l], zero);

		q[l] = CGAL::cross_product(triangle_vector[(l+1)%3], triangle_vector[l]);

		N[l] = q[l] / sqrt(q[l].squared_length());
	}

	float sum_I = 0;
	for (int k = 0; k < 3; k++)
		sum_I += s[k] * I[k];
	sum_I = (-1) * fabs(sum_I);

	//result_facet[triangle_index] = (-1) * sum;
	result_facet[point_index * no_cage_facets + triangle_index] = (-1) * sum_I;

	Vector_3 w1 = sum_I * triangle_norm;

	float x = 0, y = 0, z = 0;
	for (int k = 0; k < 3; k++)
	{
		x += II[k] * N[k].x();
		y += II[k] * N[k].y();
		z += II[k] * N[k].z();
	}
	Vector_3 w2(x, y, z);

	Vector_3 w = w1 + w2;

	float w_length = sqrt(w.squared_length());
	if(w_length > thres)
	{
		for (int l = 0; l < 3; l++)
		{
			int index = point_index * no_cage_vertices + triangle_vertices_index[l];

			if (triangle_vertices_index[l] > no_cage_vertices)
				printf("Error! Array out of range\n");

			float add_term = (N[(l+1)%3] * w) / (N[(l+1)%3] * triangle_vector[l]);

			//printf("add_term: %g\n", add_term);

			result_vertex[index] += add_term;
		}
	}

}*/

float* GreenCoordinate::gc(
	Point_3 x, int point_index, 
	vector<Point_3> mesh_coord, int no_cage_pts, 
	vector<vector<int> >mesh_trg, 
	vector<Vector_3> mesh_trg_norm, int no_cage_trgs)
{	
	const float thres = 1e-5;
	Point_3 origin(0, 0, 0);

	float* gc = new float[no_cage_pts + no_cage_trgs];

	vector<int> s(3);
	vector<float> I(3);
	vector<float> II(3);
	vector<Vector_3> q(3);
	vector<Vector_3> N(3);

	vector<Point_3> triangle_vertices(3);

	for (int j = 0; j < no_cage_trgs; j++)
	{
		vector<int> triangle_vertices_index = mesh_trg[j];
		Vector_3 triangle_norm = mesh_trg_norm[j];

		
		// move x to origin point
		for (int l = 0; l < 3; l++)	
			triangle_vertices[l] = Point_3_minus(mesh_coord[triangle_vertices_index[l]], x);

		float project_length = triangle_vertices[0].x() * triangle_norm.x() + triangle_vertices[0].y() * triangle_norm.y() + triangle_vertices[0].z() * triangle_norm.z();
		Point_3 p(project_length * triangle_norm.x(), project_length * triangle_norm.y(), project_length * triangle_norm.z());


		float sum_I = 0;
		float x = 0, y = 0, z = 0;
		Vector_3 w2(0, 0, 0);
		for (int l = 0; l < 3; l++)
		{
			// sign function
			Point_3 p1 = Point_3_minus(triangle_vertices[l], p);
			Point_3 p2 = Point_3_minus(triangle_vertices[(l+1)%3], p);

			Vector_3 v1(p1.x(), p1.y(), p1.z());
			Vector_3 v2(p2.x(), p2.y(), p2.z());

			Vector_3 cross_vector = CGAL::cross_product(v1, v2);
			float value = cross_vector * triangle_norm;
			s[l] = sign_func(value);


			// I : GCTriInt
			I[l] = TriInt(p, triangle_vertices[l], triangle_vertices[(l+1)%3], origin);


			// II : GCTriInt
			II[l] = TriInt(origin, triangle_vertices[(l+1)%3], triangle_vertices[l], origin);


			// get the normal vector of triangle (x, vi+1, vi)
			Vector_3 v3(triangle_vertices[(l+1)%3].x(), triangle_vertices[(l+1)%3].y(), triangle_vertices[(l+1)%3].z());
			Vector_3 v4(triangle_vertices[l].x(), triangle_vertices[l].y(), triangle_vertices[l].z());
			q[l] = CGAL::cross_product(v3, v4);


			// normalize the normal vector of triangle (x, vi+1, vi)
			N[l] = q[l] / sqrt(q[l].squared_length());


			sum_I += s[l] * I[l];
			w2 = w2 + II[l] * N[l];
		}
		
		sum_I = (-1) * fabs(sum_I);

		gc[no_cage_pts + j] = (-1) * sum_I;


		Vector_3 w1 = sum_I * triangle_norm;	

		Vector_3 w = w1 + w2;

		float w_length = sqrt(w.squared_length());

		if(w_length > thres)
		{
			// not co-planar
			//printf("not co-planar\n");

			for (int l = 0; l < 3; l++)
			{
				Vector_3 v(triangle_vertices[l].x(), triangle_vertices[l].y(), triangle_vertices[l].z());

				gc[triangle_vertices_index[l]] = gc[triangle_vertices_index[l]] + ((float) (N[(l+1)%3] * w) / (float) (N[(l+1)%3] * v));

				if (triangle_vertices_index[l] > no_cage_pts)
					printf("Error!!!\n");
			}

		}

	}
	
/*
	/// normalization
	float sum_v = 0;
	for (int j = 0; j < no_cage_pts; j++)
		sum_v += gc[j];
	for (int j = 0; j < no_cage_pts; j++)
		gc[j] = gc[j] / sum_v;

	float sum_f = 0;
	for (int j = 0; j < no_cage_trgs; j++)
	{
		sum_f = sum_f + gc[no_cage_pts + j];
		printf("%g\n", sum_f);
	}
	for (int j =0; j < no_cage_trgs; j++)
		gc[no_cage_pts + j] = gc[no_cage_pts + j]  / sum_f;
*/

	return gc;
}


float* GreenCoordinate::gc_all(
	vector<Point_3> interior_points, 
	vector<Point_3> mesh_coord, 
	vector<vector<int> > mesh_trg, 
	vector<Vector_3> mesh_trg_norm)
{
	/** initialization */
	int no_inter_pts = interior_points.size();
	int no_cage_pts = mesh_coord.size();
	int no_cage_trgs = mesh_trg_norm.size();

	int size_result_pts = no_inter_pts * no_cage_pts;
	int size_result_trgs = no_inter_pts * no_cage_trgs;
	int size_result = size_result_pts + size_result_trgs;

	float *result = new float[size_result];

	/** coordinate computation */
	float *gc_out;
	int point_index;
	for (int i = 0; i < no_inter_pts; i++)
	{
		point_index = i;

		gc_out = gc(interior_points[i], point_index, mesh_coord, no_cage_pts, mesh_trg, mesh_trg_norm, no_cage_trgs);

		for (int j = 0; j<no_cage_pts; j++)
			result[j + no_cage_pts*i] = gc_out[j];
		for (int j= 0; j < no_cage_trgs; j++)
			result[size_result_pts + j + no_cage_trgs*i] = gc_out[no_cage_pts + j];

		delete [] gc_out;
	}

	return result;
}

bool GreenCoordinate::go(AppObject* appObject, Mesh* origin_mesh, Mesh* cage_mesh, Mesh* deform_cage_mesh)
{
	int no_inter_pts =  origin_mesh->size_of_vertices();
	int no_cage_pts = cage_mesh->size_of_vertices();
	int no_cage_trg = cage_mesh->size_of_facets();

	cage_mesh->computeNormals();
	deform_cage_mesh->computeNormals();

	/// convert origin mesh vertices to interior_points array	
	vector<Point_3> interior_points = preprocess_vertex(origin_mesh);


	/// debug use
	//printf("interior points = %d\n\n", interior_points.size());*/
	//for (int i = 0; i < interior_points.size(); i++)
	//	printf("(%f, %f, %f)\n", interior_points[i].x(), interior_points[i].y(), interior_points[i].z());


	/// convert cage mesh vertices to cage_pcage_verticesoints array
	vector<Point_3> cage_vertices = preprocess_vertex(cage_mesh);


	/// debug use
	//printf("cage points = %d\n\n", cage_vertices.size());*/
	//for (int i = 0; i < cage_vertices.size(); i++)
	//	printf("(%f, %f, %f)\n", cage_vertices[i].x(), cage_vertices[i].y(), cage_vertices[i].z());


	/// convert cage mesh facets to cage_facets 2d array, size: no_cage_trg x 3
	vector<vector<int> > cage_facets = preprocess_facet(cage_mesh);


	///debug
	//	printf("cage facets index = %d\n\n", cage_facets.size());
	//	for (int i = 0; i < cage_facets.size(); i++)
	//		printf("(%d, %d, %d)\n", cage_facets[i][0], cage_facets[i][1], cage_facets[i][2]);
	

	/// convert cage mesh facets normal to cage_trg_norm array
	vector<Vector_3> cage_trg_norm = preprocess_facet_norm(cage_mesh);


	///debug
	//printf("cage triangle normal = %d\n\n", cage_trg_norm.size());
	//for (int i = 0; i < cage_trg_norm.size(); i++)
	//	printf("(%f, %f, %f)\n", cage_trg_norm[i].x(), cage_trg_norm[i].y(), cage_trg_norm[i].z());
	

   // allocate memory
	//float* result_vertices = new float[no_inter_pts * no_cage_pts];
	//float* result_facets = new float[no_inter_pts * no_cage_trg];


	/// calculate green coordinate, the result is store in result_vertices and result_facets
	float *gc_coord = gc_all(interior_points, cage_vertices, cage_facets, cage_trg_norm/*, result_vertices, result_facets*/);


	/// save the gc
	QString vertexGCTxtFileName = appObject->fileDir + "\\" + appObject->fileName + ".gc_v.txt"; // load cage mesh
	QString facetGCTxtFileName = appObject->fileDir + "\\" + appObject->fileName + ".gc_f.txt"; // load cage mesh

	saveGC(vertexGCTxtFileName, facetGCTxtFileName, gc_coord, no_inter_pts, no_cage_pts, no_cage_trg);


	/// convert deformed cage mesh vertices to deform_cage_vertices array
	vector<Point_3> deform_cage_vertices = preprocess_vertex(deform_cage_mesh);


	/// convert deformed cage mesh facets normal to deform_cage_trg_norm array
	vector<Vector_3> deform_cage_trg_norm = preprocess_facet_norm(deform_cage_mesh);

	vector<float> deform_scaling_factor = get_scaling_factor(cage_trg_norm, deform_cage_trg_norm);

	vector<Point_3> deform_result = postprocess(gc_coord, no_inter_pts, deform_cage_vertices, deform_cage_trg_norm, deform_scaling_factor);

	delete [] gc_coord;

	/// save the result back to the origin_mesh
	int ind = 0;
	for (Mesh::Vertex_iterator ivertex = origin_mesh->vertices_begin(); ivertex != origin_mesh->vertices_end(); ivertex++)
	{
		Mesh::Vertex_handle v = ivertex;
		v->point() = deform_result[ind];
		ind++;
	}


	return true;
}