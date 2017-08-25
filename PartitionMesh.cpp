#include "stdafx.h"

#include "PartitionMesh.h"
#include "graph.h"

#define smoothLamda 0.3
#define smoothDelta 0.001

PartitionMesh::PartitionMesh()
{
	//ctor
}

PartitionMesh::~PartitionMesh()
{
    //dtor
}

/**
void ParitionMesh :: partition(Mesh* mesh, const bool onVertices)
{
	//////////////////////////////////////////////////////////////////////////
	///// The greedy EM algorithm
	/// step 1 (initial parameters)
	int numOfPoints = (onVertices ? mesh->size_of_vertices() : mesh->size_of_facets());

	// calculate m=E[x] and S=Cov(x)=Var(x)=E[x^2]-(E[x])^2
	double mean = 0.0;
	double second_central_moment = 0.0;
	for(int i = 0 ; i < numOfPoints ; i++)
	{
		mean += results[i];
		second_central_moment += pow(results[i],2);
	}
	mean /= numOfPoints;
	second_central_moment /= numOfPoints;

	double covariance = second_central_moment - pow(mean,2);
	//

	//
	int numOfTrainingSet = 0.8 * numOfPoints;

	double beta = covariance / 2;
	double sigma = pow(0.8/numOfTrainingSet, 1/7) * beta;
	double variance = pow(sigma,2);

	vector< vector<double> > kernal_matrix;
	double constant = pow(2 * M_PI * variance, -1.5);
	for(int i = 0 ; i < numOfPoints ; i++)
	{
		for(int j = i ; j < numOfPoints ; j++)
		{
			kernal_matrix[i][j] = exp( (-0.5) * abs(results[i] - results[j]) / variance ) * constant;

			if(i!=j)
				kernal_matrix[j][i] = kernal_matrix[i][j];
		}
	}
	///

	/// step 2
	double constant = (-0.5) * log(covariance);
	double pre_log_likelihood = 0.0;
	for(int i = 0 ; i < numOfPoints ; i++)
	{
		if(covariance == 0)
			cout << "Error!! Covariance matrix is zero" << endl;

		pre_log_likelihood += constant * pow(results[i] - mean, 2) / covariance;
	}

	int clusters = 1;
	double thres = exp(-6.0);

	while(true)
	{
		if(clusters == 1)
		{

		}
		else
		{
			// E-Step (expectation)
			vector<double> adjust_weight;

			for(int i = 0; i < numOfPoints; i++)
			{
				compute1dGaussianProbability(results[i], mean, covariance);
			}


			// M-Step (maximization)
			if( abs(log_likelihood/pre_log_likelihood - 1) < thres )
				break;
			if( log_likelihood <= pre_log_likelihood)
				break;
		}

		pre_log_likelihood = log_likelihood;
	}
	///
	//////////////////////////////////////////////////////////////////////////
}
*/

/*
v*oid PartitionMesh::partition(Mesh* mesh, const bool onVertices)
{
	const int numValues = 50;
    const float interval = (float) 1/numValues;
    for ( Mesh::Facet_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++ ) {
        int intervalIndex = fit->volumeNSDF() / interval;
        if(intervalIndex == numValues) {
			fit->cluster(intervalIndex - 1);
            continue;
        }
        fit->cluster(intervalIndex);
    }
}
*/

void PartitionMesh::partition(Mesh* mesh, const bool onVertices)
{
	//////////////////////////////////////////////////////////////////////////
	///// The EM algorithm (SDF => 1D Gaussian)
	const int partitionLevels = 2;
	const int numOfPoints = (onVertices ? mesh->size_of_vertices() : mesh->size_of_facets());
	vector<double> sdfArray;
	sdfArray.resize(numOfPoints);
	int ind = 0;
	if(onVertices) {
		for(Mesh::Vertex_const_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end(); vit++) {
				//sdfArray[ind] = vit->volumeSDF();
				sdfArray[ind] = vit->volumeNSDF();
				ind++;
		}
	} else {
		for(Mesh::Facet_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) {
			//sdfArray[ind] = fit->volumeSDF();
			sdfArray[ind] = fit->volumeNSDF();
			ind++;
		}
	}

	cout << "\nPartitioning using EM algorithm, k=" << partitionLevels << endl;
	//cout << "Finding clusters within SDF values" << endl;

	// step 1. Guess Initial {alpha, mean, covariance}
	
	vector<double> alphaArray;
	vector<double> meanArray;
	vector<double> varianceArray;
	vector< vector<double>> estimateMatrix; // Matrix[partitionLevels x numOfPoints]

	alphaArray.resize(partitionLevels);
	meanArray.resize(partitionLevels);
	varianceArray.resize(partitionLevels);
	estimateMatrix.resize(partitionLevels, vector<double>(numOfPoints,0.0)); // Matrix[partitionLevels x numOfPoints]

	cout << "=====Soft partitioning=====" << endl;
    cout << "---Initial GMM parameters---" << endl;
	//cout << "Initial parameters using K-means, k = " << partitionLevels << endl;
	cout << "initial parameters -> proportion : equal, mean : pick group mean, variance : use overall variance" << endl;

    // overall mean
    double sum = 0;
    for(int i = 0; i < numOfPoints; i++)
        sum += sdfArray[i];
    double mean = (double) sum / numOfPoints;
    cout << "overall mean : " << mean << endl;

    // overall variance
    sum = 0;
    for(int i = 0 ; i < numOfPoints; i++)
        sum += pow(sdfArray[i] - mean, 2);
    double variance = (double) sum / numOfPoints;
    cout << "overall variance : " << variance << endl;

    // initial parameters
    for(int i = 0; i < partitionLevels; i++) {
        alphaArray[i] = (double) 1 / partitionLevels;
        meanArray[i] = sdfArray[rand() % numOfPoints];
        varianceArray[i] = variance;
    }

	meanArray[0] = 0.776509;
	meanArray[1] = 0.582864;
	varianceArray[0] = 0.008582;
	varianceArray[1] = 0.022687;

	cout << "Found Gaussians:" << endl;
	for(int i = 0; i < partitionLevels; i++)
		cout << "#" << i << " with alpha = " << alphaArray[i] << " , at [" << meanArray[i] << "] with variance " << varianceArray[i] << endl;

	double pre_log_likelihood = computeLogLikelihood(sdfArray, alphaArray, meanArray, varianceArray);

    vector< vector<double>> responsibilityMatrix;
    vector<double> newSumOfResponsibilityArray;
    vector<double> newMeanArray;
    vector<double> newVarianceArray;
    vector<double> newAlphaArray;
	while(true) {
	    responsibilityMatrix = evaluateResponsibility(sdfArray, alphaArray, meanArray, varianceArray);
        newSumOfResponsibilityArray = sumOfResponsibility(responsibilityMatrix);
        newMeanArray = updateMean(responsibilityMatrix, newSumOfResponsibilityArray, sdfArray);
        newVarianceArray = updateVariance(responsibilityMatrix, newSumOfResponsibilityArray, sdfArray, newMeanArray);
        newAlphaArray = updateAlpha(newSumOfResponsibilityArray, numOfPoints);

        double log_likelihood = computeLogLikelihood(sdfArray, newAlphaArray, newMeanArray, newVarianceArray);
		if( log_likelihood <= pre_log_likelihood )
            break;
		//if( abs(log_likelihood / pre_log_likelihood - 1) < exp(-6.0) )
			//break;
        alphaArray = newAlphaArray;
        meanArray = newMeanArray;
        varianceArray = newVarianceArray;
        pre_log_likelihood = log_likelihood;
	}

	//////////////////////////////////////////////////////////////////////////

	// use "Iteractive Graph Cuts for Optimal Boundary & Region Segmentation of Objects in N-D Images" //
	cout << "=====Hard partitioning=====" << endl;
	cout << "Performing K-way MinCut" << endl;
	for(int i = 0; i < numOfPoints; i++)
        for(int j = 0; j < partitionLevels; j++)
			estimateMatrix[j][i] = compute1dGaussianProbability(sdfArray[i], meanArray[j], varianceArray[j]);

	// step 1. start with an arbitrary labeling f
	cout << "Building Labels Matrix" << endl;

	vector<int> labelMatrix;
	labelMatrix.resize(numOfPoints);
	for(int i = 0; i < numOfPoints; i++) {
		double maxValue = FLT_MIN;
		double maxIndex = -1;
		for(int j = 0; j < partitionLevels; j++) {
			if(maxValue < estimateMatrix[j][i]) {
				maxValue = estimateMatrix[j][i];
				maxIndex = j;
			}
		}
		labelMatrix[i] = maxIndex;
	}

	cout << "original" << endl;
	vector<int> tmpLabelTotalVector;
	tmpLabelTotalVector.resize(partitionLevels, 0);
	for(int i = 0; i < numOfPoints; i++)
		tmpLabelTotalVector[labelMatrix[i]]++;
	for(int i = 0; i < partitionLevels; i++)
		cout << "Label " << i << ": " << tmpLabelTotalVector[i] << " facets" << endl;

	cout << "Calling K-Way Min-Cut algorithm" << endl;
	while(true) {
		// step 2. set success := 0
		int success = 0;
		// step 3. For each label alpha belongs to L
		for(int labelIndex = 0; labelIndex < partitionLevels; labelIndex++) {
			// step 3.1 Find f^ = arg min E(f`) among f` within one alpha-expansion of f
			int numOfGraphPoints = numOfPoints;

			// declare a graph and add node to this graph
			Graph *graph = new Graph();
			vector<Graph::node_id> graphNodes;
			graphNodes.resize(numOfPoints);
			for(int index = 0; index < numOfPoints; index++)
				graphNodes[index] = graph->add_node();

			//cout << "Assign weight of T-link" << endl;
			for(Mesh::Facet_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) {
				int faceIndex = fit->index();
				double souceWeight = -log10(estimateMatrix[labelIndex][faceIndex] + smoothDelta);
				double terminalWeight;
				if(labelMatrix[faceIndex] == labelIndex) // 3, 1
					terminalWeight = FLT_MAX;
				else // 3, 2
					terminalWeight = -log10(estimateMatrix[labelMatrix[faceIndex]][faceIndex] + smoothDelta);
				graph->set_tweights(graphNodes[faceIndex], souceWeight, terminalWeight);
			}

			//cout << "Create graph auxiliary edges and assign its T-link and N-link weight" << endl;
			for (Mesh::Edge_const_iterator eit = mesh->edges_begin(); eit != mesh->edges_end(); eit++) {
				Mesh::Facet_const_handle f1 = eit->facet();
				Mesh::Facet_const_handle f2 = eit->opposite()->facet();

				if (f1 == NULL || f2 == NULL) continue;

				Point_3 point1 = eit->prev()->vertex()->point();
				Point_3 point2 = eit->vertex()->point();
				int facetIndex1 = f1->index();
				int facetIndex2 = f2->index();
				if(labelMatrix[facetIndex1] == labelMatrix[facetIndex2]) { // 7
					double weight1, weight2;
					if(labelMatrix[facetIndex1] == labelIndex)
						weight1 = 0;
					else {
						double length = computePointDistance(point1, point2);
						double angle = computeDihedralAngle(f1, f2, eit);
						//weight1 = smoothLamda * length * (1 - log10(angle / M_PI));
						weight1 = smoothLamda * length * (1 - log10(angle));
					}
					graph->add_edge(graphNodes[facetIndex1], graphNodes[facetIndex2], weight1, weight1);

					if(labelMatrix[facetIndex2] == labelIndex)
						weight1 = 0;
					else {
						double length = computePointDistance(point1, point2);
						double angle = computeDihedralAngle(f1, f2, eit);
						//weight2 = smoothLamda * length * (1 - log10(angle / M_PI));
						weight2 = smoothLamda * length * (1 - log10(angle));
					}
					graph->add_edge(graphNodes[facetIndex1], graphNodes[facetIndex2], weight2, weight2);
				} else { // 4, 5, 6
					graphNodes.push_back(graph->add_node());
					numOfGraphPoints++;

					double weight1, weight2, weight3;
					if(labelMatrix[facetIndex1] == labelIndex) // 4, 5
						weight1 = 0;
					else {
						double length = computePointDistance(point1, point2);
						double angle = computeDihedralAngle(f1, f2, eit);
						//weight1 = smoothLamda * length * (1 - log10(angle / M_PI));
						weight1 = smoothLamda * length * (1 - log10(angle));
					}

					if(labelMatrix[facetIndex2] == labelIndex)
						weight2 = 0;
					else {
						double length = computePointDistance(point1, point2);
						double angle = computeDihedralAngle(f1, f2, eit);
						//weight2 = smoothLamda * length * (1 - log10(angle / M_PI));
						weight2 = smoothLamda * length * (1 - log10(angle));
					}
					graph->add_edge(graphNodes[facetIndex1], graphNodes[numOfGraphPoints-1], weight1, weight1);
					graph->add_edge(graphNodes[facetIndex2], graphNodes[numOfGraphPoints-1], weight2, weight2);

					// 6
					double length = computePointDistance(point1, point2);
					double angle = computeDihedralAngle(f1, f2, eit);
					//weight3 = smoothLamda * length * (1 - log10(angle / M_PI));
					weight3 = smoothLamda * length * (1 - log10(angle));
					graph->set_tweights(graphNodes[numOfGraphPoints - 1], 0, weight3);
				}
			}

			Graph::flowtype flow = graph->maxflow();
			cout << "Flow : "<< flow << endl;

			vector<int> tmpLabelMatrix;
			tmpLabelMatrix.resize(numOfPoints);
			for(int index = 0; index < numOfPoints; index++) {
				if (graph->what_segment(graphNodes[index]) == Graph::SOURCE)
					tmpLabelMatrix[index] = labelMatrix[index];
				else
					tmpLabelMatrix[index] = labelIndex;
			}

			delete graph;
			
			double oldEnergy = computeEnergyTerm(labelMatrix, estimateMatrix, mesh);
			double newEnergy = computeEnergyTerm(tmpLabelMatrix, estimateMatrix, mesh);

			cout << "old energy = " << oldEnergy << " , new energy = " << newEnergy << endl;

			if(newEnergy < oldEnergy) {
				for(int index = 0; index < numOfPoints; index++)
					labelMatrix[index] = tmpLabelMatrix[index];
				success = 1;
			}
		}
	
		if(success == 0)
			break;
	}

	cout << "Fixing Partition" << endl;
	vector<int> labelTotalVector;
	labelTotalVector.resize(partitionLevels, 0);
	int index = 0;
	for(Mesh::Facet_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) {
		fit->cluster(labelMatrix[index]+1);
		labelTotalVector[labelMatrix[index]]++;
		index++;
	}
	for(int i = 0; i <= partitionLevels; i++)
		cout << "Label " << i << ": " << labelTotalVector[i] << " facets" << endl;

/**
		     SOURCE
		       /         \
		   5/            \2
  	       /      3        \
	  node0 -----> node1
 	      |   <-----       |
		   |      4          |
		     \              /
		     1\          /6
		         \      /
		        SINK

	Graph::node_id nodes[2];
	Graph *g = new Graph();

	nodes[0] = g -> add_node();
	nodes[1] = g -> add_node();
	//g -> set_tweights(nodes[0], 1, 5);
	g -> set_tweights(nodes[0], 5, 1);
	g -> set_tweights(nodes[1], 2, 6);
	g -> add_edge(nodes[0], nodes[1], 3, 4);

	Graph::flowtype flow = g -> maxflow();

	printf("Flow = %d\n", flow);
	printf("Minimum cut:\n");
	if (g->what_segment(nodes[0]) == Graph::SOURCE)
		printf("node0 is in the SOURCE set\n");
	else
		printf("node0 is in the SINK set\n");
	if (g->what_segment(nodes[1]) == Graph::SOURCE)
		printf("node1 is in the SOURCE set\n");
	else
		printf("node1 is in the SINK set\n");

	delete g;
*/
}

vector< vector<double>> PartitionMesh::evaluateResponsibility(
    vector<double> sdfArray,
    vector<double> alphaArray,
    vector<double> meanArray,
    vector<double> varianceArray)
{
    int numOfPoints = sdfArray.size();
    int partitionLevels = alphaArray.size();

    vector< vector<double>> responsibilityMatrix; // Matrix [partitionLevels x numOfPoints]
    responsibilityMatrix.resize(partitionLevels, vector<double>(numOfPoints,0.0));
	//responsibilityMatrix.resize(numOfPoints, vector<double>(partitionLevels,0.0));

    for(int i = 0; i < numOfPoints; i++)
    {
        double sdf = sdfArray[i];
        double sum = 0;
        for(int j = 0; j < partitionLevels; j++)
        {
            double alpha = alphaArray[j];
            double mean = meanArray[j];
            double variance = varianceArray[j];

            double gaussian = alpha * compute1dGaussianProbability(sdf, mean, variance);
            responsibilityMatrix[j][i] = gaussian;
            sum += gaussian;
        }
        for(int j = 0; j < partitionLevels; j++)
            responsibilityMatrix[j][i] /= sum;
    }

    return responsibilityMatrix;
}

vector<double> PartitionMesh::sumOfResponsibility(vector< vector<double>> responsibilityMatrix)
{
    int partitionLevels = responsibilityMatrix.size();
    int numOfPoints = responsibilityMatrix[0].size();

    vector<double> responsibilitySum;
    responsibilitySum.resize(partitionLevels);
    for(int i = 0; i < partitionLevels; i++) {
        double sum = 0;
        for(int j = 0; j < numOfPoints; j++)
            sum += responsibilityMatrix[i][j];
        responsibilitySum[i] = sum;
    }

    return responsibilitySum;
}

vector<double> PartitionMesh::updateAlpha(
    vector<double> sumOfResponsibilityArray,
    int numOfPoints)
{
    int partitionLevels = sumOfResponsibilityArray.size();
    vector<double> alpha;
    alpha.resize(partitionLevels);
    for(int i = 0; i < partitionLevels; i++)
        alpha[i] = sumOfResponsibilityArray[i] / numOfPoints;
    return alpha;
}

vector<double> PartitionMesh::updateMean(
    vector< vector<double>> responsibilityMatrix,
    vector<double> sumOfResponsibilityArray,
    vector<double> sdfArray)
{
    int partitionLevels = responsibilityMatrix.size();
    int numOfPoints = responsibilityMatrix[0].size();

    vector<double> mean;
    mean.resize(partitionLevels);
    for(int i = 0; i < partitionLevels; i++) {
        double sum = 0;
        for(int j = 0; j < numOfPoints; j++)
            sum += responsibilityMatrix[i][j] * sdfArray[j];
        mean[i] = sum / sumOfResponsibilityArray[i];
    }

    return mean;
}

vector<double> PartitionMesh::updateVariance(
    vector< vector<double>> responsibilityMatrix,
    vector<double> sumOfResponsibilityArray,
    vector<double> sdfArray,
    vector<double> meanArray)
{
    int partitionLevels = responsibilityMatrix.size();
    int numOfPoints = responsibilityMatrix[0].size();

    vector<double> variance;
    variance.resize(partitionLevels);
    for(int i = 0; i < partitionLevels; i++)
    {
        double sum = 0;
        for(int j = 0; j < numOfPoints; j++)
            sum += responsibilityMatrix[i][j] * pow(sdfArray[j] - meanArray[i], 2);
        variance[i] = sum / sumOfResponsibilityArray[i];
    }

    return variance;
}

double PartitionMesh::compute1dGaussianProbability(
    const double sdf,
    const double mean,
    const double variance)
{
	/** 1-D Gaussian probability density function */
	double constant_term = pow( 2 * M_PI * variance, -0.5 );
	if(variance == 0)
        return 0;
	double exp_term = exp((-0.5) * pow(sdf - mean, 2) / variance);
	double gaussian = constant_term * exp_term;
	if(gaussian > 1)
		gaussian = 1;
	if(gaussian < 0)
		gaussian = 0;

	return gaussian;
}

double PartitionMesh::computeLogLikelihood(
    const vector<double> sdfArray,
    const vector<double> alphaArray,
    const vector<double> meanArray,
    const vector<double> varianceArray)
{
    /** log like-lihood */
    const int numOfPoints = sdfArray.size();
	const int numOfClusters = meanArray.size();

	double sum = 0.0;
	for(int i = 0; i < numOfPoints; i++) {
		double sdf = sdfArray[i];
		double innerSum = 0.0;
		for(int j = 0; j < numOfClusters; j++)
			innerSum += alphaArray[j] * compute1dGaussianProbability(sdf, meanArray[j], varianceArray[j]);

		sum += log(innerSum);
	}

	return sum;
}

double PartitionMesh::computeLikelihood(
    const vector<double> sdfArray,
    const vector<double> alphaArray,
    const vector<double> meanArray,
    const vector<double> varianceArray)
{
    const int numOfPoints = sdfArray.size();
	const int numOfClusters = meanArray.size();

	double sum = 0.0;
	for(int i = 0; i < numOfPoints; i++) {
		for(int j = 0; j < numOfClusters; j++)
			sum += alphaArray[j] * compute1dGaussianProbability(sdfArray[i], meanArray[j], varianceArray[j]);
    }

	return sum;
}

double PartitionMesh::computeEnergyTerm(const vector<int> labelMatrix, const vector< vector<double>> estimate, Mesh *mesh)
{
	int numOfPoints = estimate[0].size();

	double sumOfDataTerm = 0.0;
	double sumOfSmoothTerm = 0.0;
	for(int index = 0; index < numOfPoints; index++) // data term
		sumOfDataTerm += -log10( estimate[labelMatrix[index]][index] + smoothDelta );
	for (Mesh::Edge_const_iterator eit = mesh->edges_begin(); eit != mesh->edges_end(); eit++) { // smoothness term
        Mesh::Facet_const_handle f1 = eit->facet();
		Mesh::Facet_const_handle f2 = eit->opposite()->facet();

		if (f1 == NULL || f2 == NULL) continue;
		if(labelMatrix[f1->index()] == labelMatrix[f2->index()]) continue;

		Point_3 point1 = eit->prev()->vertex()->point();
		Point_3 point2 = eit->vertex()->point();

		double length = computePointDistance(point1, point2);
		double angle = computeDihedralAngle(f1, f2, eit);
		//sumOfSmoothTerm += length * (1 - log10(angle / M_PI));
		sumOfSmoothTerm += length * (1 - log10(angle));
	}

	return sumOfDataTerm + smoothLamda * sumOfSmoothTerm;
}

double PartitionMesh::computeEnergyTerm(const vector< vector<double>> estimate, Mesh *mesh)
{
	int numOfPoints = estimate[0].size();
	double sumOfDataTerm = 0.0;
	double sumOfSmoothTerm = 0.0;
	for(Mesh::Facet_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++) // data term
		sumOfDataTerm += -log10( estimate[fit->cluster()][fit->index()] + smoothDelta );
	for (Mesh::Edge_const_iterator eit = mesh->edges_begin(); eit != mesh->edges_end(); eit++) { // smoothness term
        Mesh::Facet_const_handle f1 = eit->facet();
		Mesh::Facet_const_handle f2 = eit->opposite()->facet();

		if (f1 == NULL || f2 == NULL) continue;
		if(f1->cluster() == f2->cluster()) continue;
	
		Point_3 point1 = eit->prev()->vertex()->point();
		Point_3 point2 = eit->vertex()->point();

		double length = computePointDistance(point1, point2);
		double angle = computeDihedralAngle(f1, f2, eit);
		//sumOfSmoothTerm += length * (1 - log10(angle / M_PI));
		sumOfSmoothTerm += length * (1 - log10(angle));
	}

	return sumOfDataTerm + smoothLamda * sumOfSmoothTerm;
}

double PartitionMesh::computePointDistance(const Point_3 point1, const Point_3 point2)
{
	double squareLength = pow(point1.x() - point2.x(), 2) + pow(point1.y() - point2.y(), 2) + pow(point1.z() - point2.z(), 2);
	return pow(squareLength, 1/2);
}

double PartitionMesh::computeDihedralAngle(Mesh::Facet_const_handle& f1, Mesh::Facet_const_handle& f2, Mesh::Halfedge_const_handle& edgeBetween)
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
	if (plane.oriented_side(p) == CGAL::ON_POSITIVE_SIDE) {
		concavityMultiplier = 1;
	}

	return angle * concavityMultiplier;
}

/**
template< class T >
bool isTheItemExistInSet(T item, set<T> s)
{
	if(s.count(item) == 1)
		return true;
	return false;
}
*/

/**
double calculateInnerProduct(Point_3 p1, Point_3 p2)
{
	return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}
*/

/**
double PartitionMesh::calculateInnerProduct(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	return x1 * x2 + y1 * y2 + z1 * z2;
}
*/

/**
double PartitionMesh::calculateNormalLength(const double x, const double y, const double z)
{
	double sum = pow(x,2)+pow(y,2)+pow(z,2);
	return pow(sum,1/2);
}
*/

/**
double PartitionMesh::calculateTwoFaceCommonEdgeLength(int faceIndexA, int faceIndexB, Mesh *mesh)
{
	Mesh::Facet_handle faceA = mesh->findFacet(faceIndexA);
	Mesh::Facet_handle faceB = mesh->findFacet(faceIndexB);

	Mesh::Halfedge_around_facet_const_circulator circulatorOfFaceA = faceA->facet_begin();
	Mesh::Halfedge_around_facet_const_circulator circulatorOfFaceB = faceB->facet_begin();

	Point_3 commonVertex1, commonVertex2;
	int commonCount = 0;
	do {
		Point_3 pointOfFaceA = circulatorOfFaceA->vertex()->point();

		double x1 = pointOfFaceA.x();
		double y1 = pointOfFaceA.y();
		double z1 = pointOfFaceA.z();

		do {
			Point_3 pointOfFaceB = circulatorOfFaceB->vertex()->point();

			double x2 = pointOfFaceB.x();
			double y2 = pointOfFaceB.y();
			double z2 = pointOfFaceB.z();

			if(x1==x2 && y1 == y2 && z1 == z2)
			{
				commonCount++;
				if(commonCount == 1)
					commonVertex1 = pointOfFaceA;
				if(commonCount == 2)
					commonVertex2 = pointOfFaceB;
			}
		} while (++circulatorOfFaceB != circulatorOfFaceB->facet_begin());
	} while (++circulatorOfFaceA != circulatorOfFaceA->facet_begin());

	double xA = commonVertex1.x();
	double yA = commonVertex1.y();
	double zA = commonVertex1.z();
	double xB = commonVertex2.x();
	double yB = commonVertex2.y();
	double zB = commonVertex2.z();

	return pow(pow(xA-yA, 2) + pow(yA-yB, 2) + pow(zA-zB,2), 1/2);
}
*/

/**
double PartitionMesh::calculateTwoFaceDihedralAngle(const int faceIndexA, const int faceIndexB, Mesh* mesh)
{
	// dihedral angle = 1 - inner_product(normal_face_index_A, normal_face_index_B) //
	Mesh::Facet_handle faceA = mesh->findFacet(faceIndexA);
	Mesh::Facet_handle faceB = mesh->findFacet(faceIndexB);

	double faceA_normal_x = faceA->normal().x();
	double faceA_normal_y = faceA->normal().y();
	double faceA_normal_z = faceA->normal().z();

	double faceB_normal_x = faceB->normal().x();
	double faceB_normal_y = faceB->normal().y();
	double faceB_normal_z = faceB->normal().z();

	double innerProduct = calculateInnerProduct(faceA_normal_x, faceA_normal_y, faceA_normal_z, faceB_normal_x, faceB_normal_y, faceB_normal_z);
	double faceALength = calculateNormalLength(faceA_normal_x, faceA_normal_y, faceA_normal_z);
	double faceBLength = calculateNormalLength(faceB_normal_x, faceB_normal_y, faceB_normal_z);

	double value = innerProduct / (faceALength * faceBLength);
	double angle = acos(value);

	if(angle < -1)
		angle = -1;
	if(angle > 1)
		angle = 1;

	return M_PI- angle;
}
*/

/**
list<int> PartitionMesh::getTwoFaces(const int pointIndexA, const int pointIndexB, const int numOfPoints, Mesh* mesh)
{
	list<int> faceIndexVector;
	for(int index = 0; index < numOfPoints; index++)
	{
		Mesh::Facet_handle face = mesh->findFacet(index);
		Mesh::Halfedge_around_facet_const_circulator faceIndexCirclator = face->facet_begin();

		int indexA = faceIndexCirclator->prev()->vertex()->index();
		int indexB = faceIndexCirclator->vertex()->index();
		int indexC = faceIndexCirclator->next()->vertex()->index();

		if(pointIndexA == indexA && pointIndexB == indexB)
			faceIndexVector.push_back(index);
		else if(pointIndexA == indexB && pointIndexB == indexA)
			faceIndexVector.push_back(index);
		else if (pointIndexA == indexB && pointIndexB == indexC)
			faceIndexVector.push_back(index);
		else if (pointIndexA == indexC && pointIndexB == indexB)
			faceIndexVector.push_back(index);
		else if (pointIndexA == indexC && pointIndexB == indexA)
			faceIndexVector.push_back(index);
		else if (pointIndexA == indexA && pointIndexB == indexC)
			faceIndexVector.push_back(index);

		if(faceIndexVector.size() == 2)
			return faceIndexVector;
	}

	return faceIndexVector;
}*/

/**
double PartitionMesh::useTwoPointToCalculateTwoFaceDihdedralAngle(const int pointIndexA, const int pointIndexB, const int numOfPoints, Mesh* mesh)
{
	list<int> faceList = getTwoFaces(pointIndexA, pointIndexB, numOfPoints, mesh);
	int faceIndex1 = faceList.front();
	faceList.pop_front();
	int faceIndex2 = faceList.front();
	faceList.pop_front();

	return calculateTwoFaceDihedralAngle(faceIndex1, faceIndex2, mesh);
}
*/

/**
double PartitionMesh::getNLinkWeight(const int index1, const int index2, const Point_3 point1, const Point_3 point2, const vector<int> labelMatrix, const int numOfPoints, Mesh* mesh)
{
	if(labelMatrix[index1] == labelMatrix[index2])
		return 0.0;
	double length = computePointDistance(point1, point2);
	double angle = useTwoPointToCalculateTwoFaceDihdedralAngle(index1, index2, numOfPoints, mesh);
	return length * (1 - log10(angle/M_PI));
}
*/

/**
double PartitionMesh::computeMaxWeightK(const vector<int> labelMatrix, const int numOfPoints, Mesh* mesh)
{
	// K = 1 +      max               summation          B{p,q} //
	//         p belongs to P   q:{p,q} belongs to N            //
	vector<double> sumOfFacesSmoothtermVector;
	sumOfFacesSmoothtermVector.resize(numOfPoints);

	for(int index = 0; index < numOfPoints; index++)
		sumOfFacesSmoothtermVector[index] = 0.0;
	for (Mesh::Edge_const_iterator eit = mesh->edges_begin(); eit != mesh->edges_end(); eit++) {
        Mesh::Facet_const_handle f1 = eit->facet();
		Mesh::Facet_const_handle f2 = eit->opposite()->facet();

		if (f1 == NULL || f2 == NULL) continue;

		Point_3 point1 = eit->prev()->vertex()->point();
		Point_3 point2 = eit->vertex()->point();
		int facetIndex1 = f1->index();
		int facetIndex2 = f2->index();

		if(labelMatrix[facetIndex1] != labelMatrix[facetIndex2]) {
			double length = computePointDistance(point1, point2);
			double angle = computeDihedralAngle(f1, f2, eit);
			double weight = length * (1 - log10(angle/M_PI));

			sumOfFacesSmoothtermVector[facetIndex1] = weight;
			sumOfFacesSmoothtermVector[facetIndex2] = weight;
		}
	}

	double max = 0.0;
	for(int index = 0; index < numOfPoints; index++)
		if(max < sumOfFacesSmoothtermVector[index])
			max = sumOfFacesSmoothtermVector[index];

	return 1 + max;
}
*/

/*
double calculateDataTerm(double assignProbability)
{
	double value = assignProbability + smoothDelta;
	return -log10(value);
}

double calculateSmoothnessTerm(int faceIndexA, int faceIndexB, int assignLabelIndex, vector<int> labelMatrix, Mesh* mesh)
{
	if(assignLabelIndex == labelMatrix[faceIndexB])
		return 0;

	double dihedralAngle = calculateTwoFaceDihedralAngle(faceIndexA, faceIndexB, mesh);
	double commonEdgeLength = calculateTwoFaceCommonEdgeLength(faceIndexA, faceIndexB, mesh);
	return commonEdgeLength * (1 - log10(dihedralAngle / M_PI));
}

double calculateSmoothnessTerm(int faceIndexA, int faceIndexB, int assignLabelIndexA, int assignLabelIndexB,  Mesh* mesh)
{
	if(assignLabelIndexA == assignLabelIndexB)
		return 0;

	double dihedralAngle = calculateTwoFaceDihedralAngle(faceIndexA, faceIndexB, mesh);
	double commonEdgeLength =
	(faceIndexA, faceIndexB, mesh);
	return commonEdgeLength * (1 - log10(dihedralAngle / M_PI));
}

double calculateSumOfSmoothnessTerm(int faceIndexA, int assignLabelIndex, vector<int> labelMatrix, Mesh* mesh)
{
	// find neighbors of faceIndex
	set<int> triangleIndexSet_t;


	double sumOfSmoothnessTerm = 0.0;


	//sumOfSmoothnessTerm += calculateSmoothnessTerm(faceIndexA, faceIndexB, assignLabelIndex, labelMatrix, mesh);

	return sumOfSmoothnessTerm;
}

double assignTLinkWeight(int faceIndexA, int faceIndexB, int assignLabelIndexA, int assignLabelIndexB, Mesh *mesh)
{
	return calculateSmoothnessTerm(faceIndexA, faceIndexB, assignLabelIndexA, assignLabelIndexB, mesh);
}

double assignNLinkWeight(int faceIndexA, int assignLabelIndex, vector<int> labelMatrix, vector< vector<double>> estimate, Mesh* mesh)
{
	double dataTerm = calculateDataTerm(estimate[faceIndexA][assignLabelIndex]);
	double smoothnessTerm = calculateSumOfSmoothnessTerm(faceIndexA, assignLabelIndex, labelMatrix, mesh);

	return dataTerm + smoothLamda * smoothnessTerm;
}

set<int> getFaceNeighborsIndexSet(int faceIndex, int numOfPoints, Mesh* mesh)
{
	set<int> faceIndexSet;

	Mesh::Facet_handle face = mesh->findFacet(faceIndex);
	Mesh::Halfedge_around_facet_const_circulator faceIndexCirclator = face->facet_begin();

	faceIndexSet.insert(faceIndexCirclator->vertex()->index());
	faceIndexSet.insert(faceIndexCirclator->next()->vertex()->index());
	faceIndexSet.insert(faceIndexCirclator->next()->next()->vertex()->index());

	//int numOfAuxNodes = 0;
	set<int> faceNeighborsIndexSet;
	for(int i = 0; i < numOfPoints; i++)
	{
		Mesh::Facet_handle fit = mesh->findFacet(i);
		Mesh::Halfedge_around_facet_const_circulator currentFaceCirclator = fit->facet_begin();

		int p1 = currentFaceCirclator->vertex()->index();
		int p2 = currentFaceCirclator->next()->vertex()->index();
		int p3 = currentFaceCirclator->next()->next()->vertex()->index();

		if(isTheItemExistInSet(p1, faceIndexSet) && isTheItemExistInSet(p2, faceIndexSet))
			faceNeighborsIndexSet.insert(i);
		else if(isTheItemExistInSet(p2, faceIndexSet) && isTheItemExistInSet(p3, faceIndexSet))
			faceNeighborsIndexSet.insert(i);
		else if(isTheItemExistInSet(p3, faceIndexSet) && isTheItemExistInSet(p1, faceIndexSet))
			faceNeighborsIndexSet.insert(i);
	}

	return faceNeighborsIndexSet;
}

int getNumOfFaceNeighbors(int faceIndex, int numOfPoints, Mesh* mesh)
{
	set<int> faceNeighborsSet = getFaceNeighborsIndexSet(faceIndex, numOfPoints, mesh);
	return faceNeighborsSet.size();
}

set<int> getDifferentLabelNeighbors(int faceIndex, int numOfPoints, set<int> alphaLabelSet, Mesh* mesh)
{
	set<int> faceIndexSet;

	Mesh::Facet_handle face = mesh->findFacet(faceIndex);
	Mesh::Halfedge_around_facet_const_circulator faceIndexCirclator = face->facet_begin();

	faceIndexSet.insert(faceIndexCirclator->vertex()->index());
	faceIndexSet.insert(faceIndexCirclator->next()->vertex()->index());
	faceIndexSet.insert(faceIndexCirclator->next()->next()->vertex()->index());

	//int numOfAuxNodes = 0;
	set<int> faceNeighbors;
	for(int i = 0; i < numOfPoints; i++)
	{
		Mesh::Facet_handle fit = mesh->findFacet(i);
		Mesh::Halfedge_around_facet_const_circulator currentFaceCirclator = fit->facet_begin();

		int p1 = currentFaceCirclator->vertex()->index();
		int p2 = currentFaceCirclator->next()->vertex()->index();
		int p3 = currentFaceCirclator->next()->next()->vertex()->index();

		if(isTheItemExistInSet(p1, faceIndexSet) && isTheItemExistInSet(p2, faceIndexSet))
			faceNeighbors.insert(i);
		else if(isTheItemExistInSet(p2, faceIndexSet) && isTheItemExistInSet(p3, faceIndexSet))
			faceNeighbors.insert(i);
		else if(isTheItemExistInSet(p3, faceIndexSet) && isTheItemExistInSet(p1, faceIndexSet))
			faceNeighbors.insert(i);
	}

	set<int> newFaceNeighbors;
	set<int>::iterator faceNeighborsIterator = faceNeighbors.begin();
	while(faceNeighborsIterator != faceNeighbors.end())
	{
		if(!isTheItemExistInSet(*faceNeighborsIterator, alphaLabelSet))
			newFaceNeighbors.insert(*faceNeighborsIterator);
		++faceNeighborsIterator;
	}

	return newFaceNeighbors;
}

int getNumOfDifferentLabelNeighbors(int faceIndex, int numOfPoints, set<int> alphaLabelSet, Mesh* mesh)
{
	return getDifferentLabelNeighbors(faceIndex, numOfPoints, alphaLabelSet, mesh).size();
}
*/
