/*
------------------------------------------------------------------------
Copyright LIRIS, M2DisCo 2009, Licence QPL
Author: Kai Wang
------------------------------------------------------------------------
*/

#ifndef ATTACKS
#define ATTACKS

#include "Mesh/enriched_polyhedron.h"
#include "Mesh/Polyhedron_Copy.h"

// for MSDM measurement
//#include "Distance3D.h"
//#include "Mesh/normal_cycle.h"

#include <CGAL/Subdivision_method_3.h>
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>`
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

#include <fstream>
#include <iostream>

//#include <Afx.h>
//#include <afxstr.h>

using std::ifstream;
using std::ofstream;

#include <windows.h>
#include <stdio.h>
#include <tchar.h>
#include <conio.h>

#include <QDir>

enum SimplificationType {LINDSTROMTURK, EDGELENGTHMIDPOINT};
enum SubdivisionType {CATMULLCLARK, DOOSABIN, LOOP, MIDPOINT, SQRT3};

template <class Polyhedron,class kernel>
class Attack
{
private:
	Mesh* m_MeshOriginal;
	Mesh* m_MeshAttacked;

	CCopyPoly<Polyhedron, kernel> copier;

	QString fileSeparator;
	QString attackFolderName;
	QString meshFormat;
	QString meshName;
	QString folderName;
	QString inFilename;

	QString attackFolderPath;

	double meshCentreX;
	double meshCentreY;
	double meshCentreZ;
	int previousTargetEdgeNumberSimp;
	int originalEdgeNumber;

public:

	Attack() {}

	Attack(Mesh* _mesh, QString _meshFormat, QString _meshName, QString _folderName)
	{
		MeshOriginal(_mesh);

		m_MeshAttacked = new Mesh;

		fileSeparator = "\\";
		
		meshFormat = _meshFormat;
		meshName = _meshName;
		folderName = _folderName;
		inFilename = _meshName + _meshFormat;

		attackFolderName = meshName + "_attacks";

		m_MeshOriginal->computeNormals();
		m_MeshOriginal->compute_type();

		meshCentreX = 0.0;
		meshCentreY = 0.0;
		meshCentreZ = 0.0;

		previousTargetEdgeNumberSimp = 0;
		originalEdgeNumber = (m_MeshOriginal->size_of_halfedges())/2.0;


		attackFolderPath = folderName + fileSeparator + attackFolderName;

		QDir dir(attackFolderPath);
		if (!dir.exists())
			dir.mkdir(attackFolderPath);
	}

	~Attack()
	{
		delete m_MeshAttacked;
		m_MeshAttacked = 0;
	}

	const Mesh* MeshOriginal() const { return m_MeshOriginal; }
	Mesh* MeshOriginal() { return m_MeshOriginal; }
	void MeshOriginal(Mesh* m) { m_MeshOriginal = m; }

	const Mesh* MeshAttacked() const { return m_MeshAttacked; }
	Mesh* MeshAttacked() { return m_MeshAttacked; }
	void MeshAttacked(Mesh* m) { m_MeshAttacked = m; }

	void PreProcessingAttackedMesh()
	{
		m_MeshAttacked->clear();
        copier.copy(m_MeshOriginal, m_MeshAttacked);
        m_MeshAttacked->computeNormals();
		m_MeshAttacked->compute_type();
	}

	void SaveAttackedObject(const QString& exportFileName)
	{
		//cout << "#Export file: " << qPrintable(exportFileName) << endl;

		/*if(meshFormat=="obj")*/
			m_MeshAttacked->write_obj(exportFileName);
	}

	void SaveAttackedObjectwithReordering(const QString& exportFileName, int seed)
	{
		//cout << "#Export file: " << qPrintable(exportFileName) << endl;

		/*if(meshFormat=="obj")*/
			m_MeshAttacked->writeObjRandom(exportFileName, seed);
	}

	void ElementReordering()
	{
		PreProcessingAttackedMesh();

		QString stringoutFilename = attackFolderPath + fileSeparator + meshName + "_elementreordering" + meshFormat;

		SaveAttackedObjectwithReordering(stringoutFilename, 100000);
	}

	void NoiseAdditionUniform(float noiseIntensity, bool preserveBoundaries=true)
	{
		PreProcessingAttackedMesh();

		int numVertex = (int)m_MeshAttacked->size_of_vertices();
		Vector_3 centroid = Point_3(0,0,0) - CGAL::ORIGIN;
		Mesh::Vertex_iterator pVertex;

		// mesh centre is calculated based on volume moments
		centroid = Point_3(meshCentreX,meshCentreY,meshCentreZ) - CGAL::ORIGIN;

		// calculate the average distance from the vertices to the mesh centre
		double distancetoCentroid = 0.0;
		for(pVertex=m_MeshAttacked->vertices_begin();pVertex!=m_MeshAttacked->vertices_end();pVertex++)
		{
			Vector_3 vectemp = pVertex->point() - CGAL::ORIGIN;
			distancetoCentroid = distancetoCentroid + (double)std::sqrt((vectemp - centroid) * (vectemp - centroid));
		}
		distancetoCentroid = distancetoCentroid/numVertex;

		// add pseudo-random uniformly distributed noise (between [-noiseLevel, +noiseLevel])
		unsigned int seed = (unsigned int)time(NULL) + 100000;
		srand(seed);
		double noisex, noisey, noisez;
		double noiseLevel = distancetoCentroid * noiseIntensity;
		for(pVertex=m_MeshAttacked->vertices_begin();pVertex!=m_MeshAttacked->vertices_end();pVertex++)
		{
			// keep boundaries unchanged if required by user
			bool is_border_vertex = false;
			bool stopFlag = false;
			Mesh::Halfedge_around_vertex_circulator hav = (*pVertex).vertex_begin();
			do
			{
				if(hav->is_border()==true)
				{
					is_border_vertex = true;
					stopFlag = true;
				}
				hav++;
			} while((hav!=(*pVertex).vertex_begin())&&(stopFlag==false));

			if((preserveBoundaries==true)&&(is_border_vertex==true))
				continue;

			noisex = noiseLevel * (1.0*rand()/RAND_MAX-0.5)*2;
			noisey = noiseLevel * (1.0*rand()/RAND_MAX-0.5)*2;
			noisez = noiseLevel * (1.0*rand()/RAND_MAX-0.5)*2;
			Vector_3 noise = Point_3(noisex, noisey, noisez) - CGAL::ORIGIN;
			pVertex->point() = pVertex->point() + noise;
		}

		// for correct rendering (if there is GUI), we need to update the mesh normals
		m_MeshAttacked->computeNormals();

		QString stringIntensity;
		//stringIntensity.setNum(noiseIntensity, 'e', 2);
		stringIntensity.setNum(noiseIntensity * 100 * 100);

		QString boundary = (preserveBoundaries==true) ? "keepboundary" : "notkeepboundary";

		QString stringoutFilename = attackFolderPath + fileSeparator + meshName + "_" + boundary + "_noise_" + stringIntensity + "%%" + meshFormat;

		SaveAttackedObject(stringoutFilename);
	}


	void LaplacianSmoothing(double deformFactor, int iteraNum, bool preserveBoundaries=true)
	{
		PreProcessingAttackedMesh();

		Mesh::Vertex_iterator pVertex;
		int numVertex = (int)m_MeshAttacked->size_of_vertices();
		Vector_3 * newPositions = new Vector_3[numVertex];

		for(int i=0;i<iteraNum;i++)
		{
			int n = 0;
			for(pVertex=m_MeshAttacked->vertices_begin();pVertex!=m_MeshAttacked->vertices_end();pVertex++)
			{
				Vector_3 currentVector_3 = pVertex->point() - CGAL::ORIGIN;

				// do not smooth the boundary vertices if required by user
				bool is_border_vertex = false;
				bool stopFlag = false;
				Mesh::Halfedge_around_vertex_circulator hav = (*pVertex).vertex_begin();
				do
				{
					if(hav->is_border()==true)
					{
						is_border_vertex = true;
						stopFlag = true;
					}
					hav++;
				} while((hav!=(*pVertex).vertex_begin())&&(stopFlag==false));

				if((preserveBoundaries==true)&&(is_border_vertex==true))
				{
					newPositions[n] = currentVector_3;
					n++;
					continue;
				}

				std::size_t degree = (*pVertex).vertex_degree();
				double alpha = 1.0/degree;
				Vector_3 vectemp = Point_3(0,0,0) - CGAL::ORIGIN;
				Mesh::Halfedge_around_vertex_circulator h = (*pVertex).vertex_begin();
				do
				{
					vectemp = vectemp+(h->opposite()->vertex()->point()-CGAL::ORIGIN-currentVector_3)*alpha;
					++h;
				} while(h != (*pVertex).vertex_begin());
				newPositions[n] = currentVector_3 + deformFactor*vectemp;
				n++;
			}

			n = 0;
			for(pVertex=m_MeshAttacked->vertices_begin();pVertex!=m_MeshAttacked->vertices_end();pVertex++)
			{
				pVertex->point() = Point_3(0,0,0) + newPositions[n];
				n++;
			}
		}
		
		delete [] newPositions;
		newPositions = 0;

		m_MeshAttacked->computeNormals();

		QString stringLambda, stringIteration;
		//stringLambda.setNum(deformFactor, 'e', 2);
		stringLambda.setNum(deformFactor * 100);
		stringIteration.setNum(iteraNum, 10);
		
		QString boundary = (preserveBoundaries ? "keepboundary" : "notkeepboundary");
		
		QString stringoutFilename = attackFolderPath + fileSeparator + meshName + "_smoothing_iteration_" + stringIteration + meshFormat;
		
		SaveAttackedObject(stringoutFilename);
	}

	// TODO: fixing after quantization the potentially introduced degeneracies, such as removing the null surface facet
	void CoordinateQuantization(int bitDepth)
	{
		PreProcessingAttackedMesh();

		Mesh::Vertex_iterator pVertex;
		double quantizationLevel = std::pow(2.0,bitDepth);
		pVertex = m_MeshAttacked->vertices_begin();
		Point_3 point = pVertex->point();
		double xmax = double(point.x());
		double xmin = xmax;
		double ymax = double(point.y());
		double ymin = ymax;
		double zmax = double(point.z());
		double zmin = zmax;
		pVertex++;

		for(;pVertex!=m_MeshAttacked->vertices_end();pVertex++)
		{
			point = pVertex->point();
			double x = double(point.x());
			double y = double(point.y());
			double z = double(point.z());
			if (x>xmax)
				xmax = x;
			if (x<xmin)
				xmin = x;
			if (y>ymax)
				ymax = y;
			if (y<ymin)
				ymin = y;
			if (z>zmax)
				zmax = z;
			if (z<zmin)
				zmin = z;
		}

		double xstep = (xmax-xmin)/quantizationLevel;
		double ystep = (ymax-ymin)/quantizationLevel;
		double zstep = (zmax-zmin)/quantizationLevel;

		for(pVertex=m_MeshAttacked->vertices_begin();pVertex!=m_MeshAttacked->vertices_end();pVertex++)
		{
			point = pVertex->point();
			double x = double(point.x());
			double y = double(point.y());
			double z = double(point.z());

			double xquantified, yquantified, zquantified;

			double xint = 1.0*std::floor((x-xmin)/xstep)*xstep + xmin;
			double xfrac = x - xint;
			if (xfrac<=(0.5*xstep))
				xquantified = xint;
			else
				xquantified = xint + xstep;

			double yint = 1.0*std::floor((y-ymin)/ystep)*ystep + ymin;
			double yfrac = y - yint;
			if (yfrac<=(0.5*ystep))
				yquantified = yint;
			else
				yquantified = yint +ystep;

			double zint = 1.0*std::floor((z-zmin)/zstep)*zstep + zmin;
			double zfrac = z - zint;
			if (zfrac<=(0.5*zstep))
				zquantified = zint;
			else
				zquantified = zint + zstep;

			pVertex->point() = Point_3(xquantified,yquantified,zquantified);
		}

		m_MeshAttacked->computeNormals();
		
		QString stringBitDepth;
		stringBitDepth.setNum(bitDepth, 10);

		QString stringoutFilename = attackFolderPath + fileSeparator + meshName + "_quantization_bitdepth_" + stringBitDepth + meshFormat;

		SaveAttackedObject(stringoutFilename);
	}

	void SimilarityTransformation()
	{
		PreProcessingAttackedMesh();

		unsigned int seed = (unsigned int)time(NULL) + 100000;

		int numVertex = (int)m_MeshAttacked->size_of_vertices();
		Vector_3 centroid = Point_3(0,0,0) - CGAL::ORIGIN;
		Mesh::Vertex_iterator pVertex;

		// mesh centre is calculated based on volume moments
		centroid = Point_3(meshCentreX,meshCentreY,meshCentreZ) - CGAL::ORIGIN;

		double distancetoCentroid = 0.0;
		for(pVertex=m_MeshAttacked->vertices_begin();pVertex!=m_MeshAttacked->vertices_end();pVertex++)
		{
			Vector_3 vectemp = pVertex->point() - CGAL::ORIGIN;
			distancetoCentroid = distancetoCentroid + (double)std::sqrt((vectemp - centroid) * (vectemp - centroid));
		}
		distancetoCentroid = distancetoCentroid/numVertex;

		double xTrans = (1.0*rand()/RAND_MAX-0.5)*10.0*distancetoCentroid;
		double yTrans = (1.0*rand()/RAND_MAX-0.5)*10.0*distancetoCentroid;
		double zTrans = (1.0*rand()/RAND_MAX-0.5)*10.0*distancetoCentroid;
		Translation(xTrans,yTrans,zTrans);

		double xRot = (1.0*rand()/RAND_MAX-0.5)*2.0;
		double yRot = (1.0*rand()/RAND_MAX-0.5)*2.0;
		double zRot = (1.0*rand()/RAND_MAX-0.5)*2.0;
		double angleRot = (1.0*rand()/RAND_MAX)*360.0;
		Rotation(xRot,yRot,zRot,angleRot);

		double factorScale = (1.0*rand()/RAND_MAX)*9.9 + 0.1;
		UniformScaling(factorScale);
		
		QString stringoutFilename = attackFolderPath + fileSeparator + meshName + "_similaritytransform" + meshFormat;

		SaveAttackedObject(stringoutFilename);
	}

	void Translation(double xTranslation,double yTranslation,double zTranslation)
	{
		Vector_3 translationVector_3(xTranslation,yTranslation,zTranslation);
		Aff_transformation_3 translation(CGAL::TRANSLATION,translationVector_3);
		std::transform(m_MeshAttacked->points_begin(),m_MeshAttacked->points_end(),m_MeshAttacked->points_begin(),translation);
		m_MeshAttacked->computeNormals();
	}

	void Translation_Save(double xTranslation, double yTranslation, double zTranslation)
	{
		PreProcessingAttackedMesh();

		Translation(xTranslation, yTranslation, zTranslation);

		QString stringoutFilename = folderName + fileSeparator + meshName + "_translation" + meshFormat;

		SaveAttackedObject(stringoutFilename);
	}

	void Rotation(double xAxis,double yAxis,double zAxis,double angle)
	{
		// normalize the rotation axis Vector_3
		double normAxis = sqrt(xAxis*xAxis + yAxis*yAxis + zAxis*zAxis);
		xAxis = xAxis / normAxis;
		yAxis = yAxis / normAxis;
		zAxis = zAxis / normAxis;

		// construction of the rotation matrix
		double c = cos(angle/180.0*M_PI);
		double s = sin(angle/180.0*M_PI);
		double m00 = xAxis*xAxis + (1.0-xAxis*xAxis)*c;
		double m01 = xAxis*yAxis*(1.0-c) - zAxis*s;
		double m02 = xAxis*zAxis*(1.0-c) + yAxis*s;
		double m10 = xAxis*yAxis*(1.0-c) + zAxis*s;
		double m11 = yAxis*yAxis + (1.0-yAxis*yAxis)*c;
		double m12 = yAxis*zAxis*(1.0-c) - xAxis*s;
		double m20 = xAxis*zAxis*(1.0-c) - yAxis*s;
		double m21 = yAxis*zAxis*(1.0-c) + xAxis*s;
		double m22 = zAxis*zAxis + (1.0-zAxis*zAxis)*c;

		// realize rotation by applying general affine transformation through matrix multiplication
		Aff_transformation_3 rotation(m00,m01,m02,m10,m11,m12,m20,m21,m22);
		std::transform(m_MeshAttacked->points_begin(),m_MeshAttacked->points_end(),m_MeshAttacked->points_begin(),rotation);
		m_MeshAttacked->computeNormals();
	}

	void Rotation_Save(double xAxis,double yAxis,double zAxis,double angle)
	{
		PreProcessingAttackedMesh();

		Rotation(xAxis, yAxis, zAxis, angle);

		QString stringoutFilename = folderName + fileSeparator + meshName + "_rotation" + meshFormat;

		SaveAttackedObject(stringoutFilename);
	}

	void UniformScaling(double scalingFactor)
	{
		Aff_transformation_3 uniformScaling(CGAL::SCALING, scalingFactor);
		std::transform(m_MeshAttacked->points_begin(), m_MeshAttacked->points_end() ,m_MeshAttacked->points_begin(), uniformScaling);
		m_MeshAttacked->computeNormals();
	}

	void UniformScaling_Save(double scalingFactor)
	{
		PreProcessingAttackedMesh();

		UniformScaling(scalingFactor);

		QString stringoutFilename = attackFolderPath + fileSeparator + meshName + "_uniformscaling" + meshFormat;

		SaveAttackedObject(stringoutFilename);
	}

	void NonUniformScaling(double scalingFactorX, double scalingFactorY, double scalingFactorZ)
	{
		for(Mesh::Vertex_iterator pVertex = m_MeshAttacked->vertices_begin(); pVertex != m_MeshAttacked->vertices_end(); pVertex++)
		{	
			double x = pVertex->point().x();
			double y = pVertex->point().y();
			double z = pVertex->point().z();

			double newX = x * scalingFactorX;
			double newY = y * scalingFactorY;
			double newZ = z * scalingFactorZ;

			Point_3 newPoint = Point_3(newX, newY, newZ); /*CGAL::ORIGIN;*/
			pVertex->point() = newPoint; /*pVertex->point() + noise;*/
		}

		m_MeshAttacked->computeNormals();
	}

	void NonUniformScaling_Save(double scalingFactorX, double scalingFactorY, double scalingFactorZ)
	{
		PreProcessingAttackedMesh();

		NonUniformScaling(scalingFactorX, scalingFactorY, scalingFactorZ);

		QString stringoutFilename = attackFolderPath + fileSeparator + meshName + "_nonuniformscaling" + meshFormat;

		SaveAttackedObject(stringoutFilename);
	}

	void Subdivision(SubdivisionType subdivisionType, int depth)
	{
		PreProcessingAttackedMesh();

		if((subdivisionType==LOOP)||(subdivisionType==SQRT3)||(subdivisionType==MIDPOINT))
		{
			if(!m_MeshAttacked->is_pure_triangle())
			{
				m_MeshAttacked->triangulate();
				m_MeshAttacked->computeNormals();
				m_MeshAttacked->compute_type();
			}
		}

		if(subdivisionType==CATMULLCLARK)
			SubdivisionCatmullClark(depth);
		else if(subdivisionType==LOOP)
			SubdivisionLoop(depth);
		else if(subdivisionType==DOOSABIN)
			SubdivisionDooSabin(depth);
		else if(subdivisionType==SQRT3)
			SubdivisionSqrt3(depth);
		else if(subdivisionType==MIDPOINT)
			SubdivisionMidpoint(depth);

		m_MeshAttacked->compute_type();
		m_MeshAttacked->computeNormals();

		QString stringType;
		if(subdivisionType==CATMULLCLARK)
			stringType = "catmullclark";
		else if(subdivisionType==LOOP)
			stringType = "loop";
		else if(subdivisionType==DOOSABIN)
			stringType = "doosabin";
		else if(subdivisionType==SQRT3)
			stringType = "sqrt3";
		else if(subdivisionType==MIDPOINT)
			stringType = "midpoint";
		stringType.toLower();

		QString stringDepth;
		stringDepth.setNum(depth, 10);

		QString stringoutFilename = attackFolderPath + fileSeparator + meshName + "_subdivision_" + stringType + "_iteration_"+ stringDepth + meshFormat;

		SaveAttackedObject(stringoutFilename);
	}

	void SubdivisionCatmullClark(int depth)
	{
		CGAL::Subdivision_method_3::CatmullClark_subdivision(*m_MeshAttacked,depth);
	}

	void SubdivisionLoop(int depth)
	{
		CGAL::Subdivision_method_3::Loop_subdivision(*m_MeshAttacked,depth);
	}

	void SubdivisionDooSabin(int depth)
	{
		CGAL::Subdivision_method_3::DooSabin_subdivision(*m_MeshAttacked,depth);
	}

	void SubdivisionSqrt3(int depth)
	{
		CGAL::Subdivision_method_3::Sqrt3_subdivision(*m_MeshAttacked,depth);
	}

	// The midPoint_3 subdivision does not introduce any geometry distortion, it simply adds vertices in the middle of edges.
	void SubdivisionMidpoint(int depth)
	{
		CGAL::Subdivision_method_3::PTQ(*m_MeshAttacked,CGAL::Linear_mask_3<Mesh>(),depth);
	}

	void Simplification(SimplificationType simpType, double removalRatio)
	{
		int targetEdgeNum = ceil((1.0-removalRatio/100.0) * 1.0 * originalEdgeNumber);

		if (targetEdgeNum>previousTargetEdgeNumberSimp)
		{
			PreProcessingAttackedMesh();
		}

		if (!m_MeshAttacked->is_pure_triangle())
		{
			m_MeshAttacked->triangulate();	
		}

		m_MeshAttacked->computeNormals();
		m_MeshAttacked->compute_type();

		CGAL::Surface_mesh_simplification::Count_stop_predicate<Surface> stop(targetEdgeNum);

		// upcast...
		Surface * surface = (Surface*)(m_MeshAttacked);

		if(simpType==LINDSTROMTURK)
		{
			int r = CGAL::Surface_mesh_simplification::edge_collapse
					(*surface
					,stop
					,CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index,*surface))
						  .edge_index_map  (boost::get(CGAL::edge_external_index  ,*surface))
					);
		}
		else if(simpType==EDGELENGTHMIDPOINT)
		{
			int r = CGAL::Surface_mesh_simplification::edge_collapse
					(*surface
					,stop
					,CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index,*surface))
						  .edge_index_map  (boost::get(CGAL::edge_external_index  ,*surface))
						  .get_cost(CGAL::Surface_mesh_simplification::Edge_length_cost<Surface>())
						  .get_placement(CGAL::Surface_mesh_simplification::Midpoint_placement<Surface>())
					);
		}

		m_MeshAttacked->compute_type();
		m_MeshAttacked->computeNormals();

		QString stringType;
		if(simpType==EDGELENGTHMIDPOINT)
			stringType = "edgelengthmidpoint";
		else if(simpType==LINDSTROMTURK)
			stringType = "lindstromturk";
		stringType.toLower();

		QString stringRatio;
		//stringRatio.setNum(removalRatio, 'e', 2);
		stringRatio.setNum(removalRatio);

		QString stringoutFilename = attackFolderPath + fileSeparator + meshName + "_simplification_" + stringType + "_ratio_" + stringRatio + "%" + meshFormat;

		SaveAttackedObject(stringoutFilename);

		previousTargetEdgeNumberSimp = targetEdgeNum;
	}

};


#endif // ATTACKS