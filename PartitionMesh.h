#ifndef PARTITIONMESH_H
#define PARTITIONMESH_H

#include "stdafx.h"

class PartitionMesh
{
    public:
        PartitionMesh();
        virtual ~PartitionMesh();
        void partition(Mesh* mesh, const bool onVertices);
    protected:
    private:
        vector< vector<double>> evaluateResponsibility(vector<double> sdfArray, vector<double> alphaArray, vector<double> meanArray, vector<double> varianceArray);
        vector<double> sumOfResponsibility(vector< vector<double>> responsibilityMatrix);
        vector<double> updateAlpha(vector<double> sumOfResponsibilityArray, int numOfPoints);
        vector<double> updateMean(vector< vector<double>> responsibilityMatrix, vector<double> sumOfResponsibilityArray, vector<double> sdfArray);
        vector<double> updateVariance(vector< vector<double>> responsibilityMatrix, vector<double> sumOfResponsibilityArray, vector<double> sdfArray, vector<double> meanArray);
        double compute1dGaussianProbability(const double x, const double mean, const double variance);
        double computeLogLikelihood(const vector<double> sdfArray, const vector<double> alphaArray, const vector<double> meanArray, const vector<double> varianceArray);
        double computeLikelihood(const vector<double> sdfArray, const vector<double> alphaArray, const vector<double> meanArray, const vector<double> varianceArray);
		double computePointDistance(const Point_3 point1, const Point_3 point2);
		double computeDihedralAngle(Mesh::Facet_const_handle& f1, Mesh::Facet_const_handle& f2, Mesh::Halfedge_const_handle& edgeBetween);
		double computeEnergyTerm(const vector< vector<double>> estimate, Mesh *mesh);
		double computeEnergyTerm(const vector<int> labelMatrix, const vector< vector<double>> estimate, Mesh *mesh);

        /** list<int> getTwoFaces(const int pointIndexA, const int pointIndexB, const int numOfPoints, Mesh* mesh); */
        /** double useTwoPointToCalculateTwoFaceDihdedralAngle(const int pointIndexA, const int pointIndexB, const int numOfPoints, Mesh* mesh); */
        /** double calculateNormalLength(const double x, const double y, const double z); */
        /** double calculateTwoFaceCommonEdgeLength(int faceIndexA, int faceIndexB, Mesh *mesh); */
        /** double calculateInnerProduct(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2); */
        /** double calculateTwoFaceDihedralAngle(const int faceIndexA, const int faceIndexB, Mesh* mesh); */
        /** double getNLinkWeight(const int index1, const int index2, const Point_3 point1, const Point_3 point2, const vector<int> labelMatrix, const int numOfPoints, Mesh* mesh); */
		/**double computeMaxWeightK(const vector<int> labelMatrix, const int numOfPoints, Mesh* mesh); */
};

#endif // PARTITIONMESH_H
