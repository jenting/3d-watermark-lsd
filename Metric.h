#ifndef __METRIC_H
#define __METRIC_H

#include "stdafx.h"

#include <atlstr.h>

class MSDM {

private:
	bool sphereClipVector(Point_3 &O, double r, const Point_3 &P, Vector_3 &V);

	double processRoughnessPerVertexCurve(
		Mesh* mesh, Mesh::Vertex* pVertex, double radius, 
		std::vector<double> &TabDistance, std::vector<Point_3> &TabPoint, 
		double & meanRet, double dim, bool IsGauss = true);
	
	double processCovariance(
		Mesh::Vertex* pVertex, double curvature, double curvatureDeg,
		std::vector<double> TabDistance, std::vector<Point_3> TabPoint,
		std::vector<double> TabDistanceDeg, std::vector<Point_3> TabPointDeg, 
		double dim, bool IsGauss = true);

	void computeMSDMDistance(
		Mesh* originalMesh, Mesh* watermarkedMesh,
		double Param, double &L1, double &L2, double &L3, double &L4);

public:
	double calculateMSDM(Mesh* originalMesh, Mesh* watermarkedMesh, double radius);

	double processRoughnessCurveDual(
		Mesh* originalMesh, Mesh* watermarkedMesh,
		double radius, double maxdim, bool IsGauss = true);

};

class DeformationLevel {

private:
	double mean(double mean_p, double mean_q);
	double deviation(double dev_p, double dev_q);
	double covariance(double cov_p, double cov_q, double cov_pq);

public:
	void calculateVSIPointDeformationLevel(Mesh* originalMesh, Mesh* watermarkedMesh);
	void calculateSDFPointDeformationLevel(Mesh* originalMesh, Mesh* watermarkedMesh);
	void calculateNSDFPointDeformationLevel(Mesh* originalMesh, Mesh* watermarkedMesh);

	double calculateLocalWindowDeformationLevel(Mesh* originalMesh, Mesh* watermarkedMesh);
	double calculateGlobalWindowsDeformationLevel(Mesh* originalMesh, Mesh* watermarkedMesh);

};

class Quality {

private:
	void TrimString(CString & aString);
		
public:
	/** Signal-to-Noise ratio */
	double measureSNR(Mesh* originalMesh, Mesh* watermarkedMesh);

	/** Peak Signal-to-Noise ratio */
	double measurePSNR(Mesh* originalMesh, Mesh* watermarkedMesh);
	
	void calculateGeoDistances(
		const QString& originalMeshFileName, const QString& watermarkedMeshFileName, 
		double& mrms, double& mrmsBB, double& hausdorff, double& hausdorffwrtBB);

};

class Statistic {

private:
	
public:
	/** Bit Error Rate (BER) */
	double calculateBitErrorRate(const int* seqA, const int sizeA, const int* seqB, const int sizeB);
	double calculateBitErrorRate(const TNT::Array1D<int>& seqA, const TNT::Array1D<int>& seqB);

	/** Correlation Coefficient (CC) */
	double calcuateCorrelationCoefficient(const int* seqA, const int sizeA, const int* seqB, const int sizeB);
	double calcuateCorrelationCoefficient(const TNT::Array1D<int>& seqA, const TNT::Array1D<int>& seqB);
	
};

#endif