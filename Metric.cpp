#include "stdafx.h"

#include <algorithm>
#include <atlstr.h>
#include <conio.h>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <map>
#include <set>
#include <stack>
#include <stdio.h>
#include <string>
#include <tchar.h>
#include <windows.h>

using std::ifstream;
using std::ofstream;

#include "Metric.h"
#include "Utility.h"

#ifndef MAX
#define MAX(a,b) (a>b?a:b)
#endif

bool MSDM::sphereClipVector(Point_3 &O, double r, const Point_3 &P, Vector_3 &V)
{
    Vector_3 W = P - O;
    double a = (V*V);
    double b = 2.0 * V * W ;
    double c = (W*W) - r*r ;
    double delta = b*b - 4*a*c ;

	if( a==0)
		return true ;

    if(delta < 0) // Should not happen, but happens sometimes (numerical precision)	
        return true ;

    double t = (- b + std::sqrt(delta)) / (2.0 * a) ;

    if(t < 0.0) // Should not happen, but happens sometimes (numerical precision)
		return true ;
    if(t >= 1.0) // Inside the sphere
        return false ;

	if(t==0)
		t=0.01;

    V=V*t;

    return true ;
}

double MSDM::calculateMSDM(Mesh* originalMesh, Mesh* watermarkedMesh, double radius)
{
	double maxDim = originalMesh->getMaxDim();
	double maxDim2 = watermarkedMesh->getMaxDim();
	
	printf("precess roughness curva dual\n");
	processRoughnessCurveDual(originalMesh, watermarkedMesh, radius * maxDim, maxDim);
	printf("precess roughness curva dual done\n");

	double param = 0.5;
	double L1, L2, L3, L4;

	printf("compute MSDM\n");
	computeMSDMDistance(originalMesh, watermarkedMesh, param, L1, L2, L3, L4);

	printf("L1 = %g\n", L1);
	printf("L2 = %g\n", L2);
	printf("L3 = %g\n", L3);
	printf("L4 = %g\n", L4);

	return L3;
}

void MSDM::computeMSDMDistance(
	Mesh* originalMesh,	Mesh* watermarkedMesh,
	double param, double &L1, double &L2, double &L3, double &L4)
{
	double SommeDist=0;
	double SommeDistCarre=0;
	double SommeDist3=0;
	double SommeDist4=0;
	double Max=0;
	int NbVertex=0;

	/*double PolyDegradRoughMin=m_PolyDegrad.MinNrmMeanCurvature();
	double PolyDegradRoughMax=m_PolyDegrad.MaxNrmMeanCurvature();

	double PolyOriginalRoughMin=m_PolyOriginal->MinNrmMeanCurvature();
	double PolyOriginalRoughMax=m_PolyOriginal->MaxNrmMeanCurvature();*/

	/*m_PolyOriginal->MinNrmGaussCurvature(10000);
	m_PolyOriginal->MaxNrmGaussCurvature(0);*/

	//int CMPT=0;

	Mesh::Vertex_iterator pVertex = originalMesh->vertices_begin();
	Mesh::Vertex_iterator pVertex_end = watermarkedMesh->vertices_end();
	Mesh::Vertex_iterator pVertexDef = watermarkedMesh->vertices_begin();	

	for (; pVertex != pVertex_end; pVertex++)
	{
		/////////////////////:avec covariance///////////////////////////////////::
		double MoyX=pVertex->curvatureMean;
		double MoyY=pVertexDef->curvatureMean;
		double SigX=pVertex->curvatureVariance;
		double SigY=pVertexDef->curvatureVariance;
		double SigXY=pVertex->curvatureCovariance;

		double angleRad;
		
		/// minkowsky
		double fact1=(fabs(MoyX-MoyY)) / (MAX(MoyX, MoyY) + 1); //(2*MoyX*MoyY)/(MoyX*MoyX+MoyY*MoyY);
		double fact2=(fabs(SigX-SigY)) / (MAX(SigX, SigY) + 1); //2*SigX*SigY/(SigX*SigX+SigY*SigY);
		double fact3=(fabs(SigX*SigY-SigXY)) / (SigX*SigY + 1); //fabs(SigXY/(SigX*SigY));

		angleRad=pow( (pow(fact1,3) + pow(fact2,3) + param * pow(fact3,3)) / (2.+param), 1./3.);
					
		if(MoyX<1)
			int h=2;

		//CMPT++;

		SommeDist+=angleRad;
		SommeDistCarre+=angleRad*angleRad;
		SommeDist3+=angleRad*angleRad*angleRad;
		SommeDist4+=angleRad*angleRad*angleRad*angleRad;
		
		/*if(Max<angleRad)
			Max=angleRad;

		double Indice=angleRad;

		pVertex->Kg(Indice);
		if(Indice<m_PolyOriginal->MinNrmGaussCurvature())
			m_PolyOriginal->MinNrmGaussCurvature(Indice);

		if(Indice>m_PolyOriginal->MaxNrmGaussCurvature())
			m_PolyOriginal->MaxNrmGaussCurvature(Indice);*/

		NbVertex++;

		pVertexDef++;
	}

	L1=SommeDist / (double)NbVertex;
	L2=sqrt(SommeDistCarre / (double)NbVertex);
	L3=pow(SommeDist3 / (double)NbVertex, 0.33333);
	L4=pow(SommeDist4 / (double)NbVertex, 0.25);
}

double MSDM::processRoughnessCurveDual(
	Mesh* originalMesh, Mesh* watermarkedMesh,
	double radius, double maxdim, bool IsGauss)
{
	Mesh::Vertex_iterator pVertex = originalMesh->vertices_begin();
	Mesh::Vertex_iterator pVertexDef=watermarkedMesh->vertices_begin();
	Mesh::Vertex_iterator pVertex_end = originalMesh->vertices_end();

	for (; pVertex != pVertex_end; pVertex++)
	{
		std::vector<double> TabDistance;
		std::vector<Point_3> TabPoint;

		std::vector<double> TabDistanceDeg;
		std::vector<Point_3> TabPointDeg;

		double mean, meanDeg;

		double var1 = processRoughnessPerVertexCurve(originalMesh, (&(*pVertex)), radius, TabDistance, TabPoint, mean, maxdim, IsGauss);

		double var2 = processRoughnessPerVertexCurve(watermarkedMesh, (&(*pVertexDef)), radius, TabDistanceDeg, TabPointDeg, meanDeg, maxdim, IsGauss);

		double cov = processCovariance((&(*pVertexDef)), mean, meanDeg, TabDistance, TabPoint, TabDistanceDeg, TabPointDeg, maxdim, IsGauss);

		//double cov=1;
		pVertexDef->curvatureCovariance = pVertex->curvatureCovariance = cov;

		pVertexDef++;
	}
	
	return 0;
}

double MSDM::processRoughnessPerVertexCurve(
	Mesh* mesh, Mesh::Vertex* pVertex, 
	double radius, 	std::vector<double> &TabDistance, std::vector<Point_3> &TabPoint, 
	double & meanRet, double dim, bool IsGauss)
{
	Point_3 O = pVertex->point();

	std::set<Mesh::Vertex*> vertices;
    std::stack<Mesh::Vertex*> S;

    vertices.insert(pVertex);
	S.push(pVertex);	

	int numInSphere=0;
	double sumDist=0;

    while(!S.empty())
	{
		Mesh::Vertex* v = S.top() ;
        S.pop() ;

        Point_3 P = v->point() ;

        Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
		Mesh::Halfedge_around_vertex_circulator pHalfedgeStart = h;
		CGAL_For_all(h, pHalfedgeStart)
		{
            Point_3 p1 = h->vertex()->point();
			Point_3 p2 = h->opposite()->vertex()->point();
			Vector_3 V = (p2-p1);

            if(v == pVertex || V * (P - O) > 0.0) 
			{
				double lenOld = std::sqrt(V*V);
				bool isect = sphereClipVector(O, radius, P, V) ;
				double lenEdge = std::sqrt(V*V);

				numInSphere++;

				///ici on prend en compte la distance map des sommets
				
				double weightedDist;
				Point_3 PPondere;

				if(lenOld != 0)
				{
					//weightedDist=(1 - lenEdge/lenOld) * h->vertex()->Kmax() + (lenEdge / lenOld )* h->opposite()->vertex()->Kmax();
					weightedDist=(1 - lenEdge/lenOld) * h->vertex()->volumeMinDist() + (lenEdge / lenOld )* h->opposite()->vertex()->volumeMinDist();
					//weightedDist=(1 - lenEdge/lenOld) * h->vertex()->volumeMaxDist() + (lenEdge / lenOld )* h->opposite()->vertex()->volumeMaxDist();
					//weightedDist=(1 - lenEdge/lenOld) * h->vertex()->volumeVSI() + (lenEdge / lenOld )* h->opposite()->vertex()->volumeVSI();
					//weightedDist=(1 - lenEdge/lenOld) * h->vertex()->volumeSDF() + (lenEdge / lenOld )* h->opposite()->vertex()->volumeSDF();
					//weightedDist=(1 - lenEdge/lenOld) * h->vertex()->volumeNSDF() + (lenEdge / lenOld )* h->opposite()->vertex()->volumeNSDF();


					PPondere=p1+V;
				}
				else
				{
					//weightedDist=h->opposite()->vertex()->Kmax();
					weightedDist=h->opposite()->vertex()->volumeMinDist();
					//weightedDist=h->opposite()->vertex()->volumeMaxDist();
					//weightedDist=h->opposite()->vertex()->volumeVSI();
					//weightedDist=h->opposite()->vertex()->volumeSDF();
					//weightedDist=h->opposite()->vertex()->volumeNSDF();

					PPondere=p2;
				}

				TabDistance.push_back(weightedDist);
				TabPoint.push_back(PPondere);

				sumDist+=weightedDist;

				if(!isect) 
				{
					Mesh::Vertex_iterator w=h->opposite()->vertex();
                    if(vertices.find(&(*w)) == vertices.end())
					{
                        vertices.insert(&(*w)) ;
                        S.push(&(*w)) ;
                    }
                }
				
			}
            
		}
		
	}

	
	/// calculate mean and variance in a local window
	double mean=0;
	double variance=0;

	if(!IsGauss)
	{
		mean=sumDist / (double)numInSphere;
		variance=0;
		
		if(numInSphere != 0)
		{
			for(int i=0;i<numInSphere;i++)
				variance += pow(TabDistance[i] - mean, 2);

			variance = variance / (double)numInSphere;
			variance = sqrt(variance);
		}
	}
	else
	{
		//variance = 0.008
		sumDist=0;
		double sumWi = 0;
		for(int i=0; i<numInSphere; i++)
		{
			Vector_3 DistancePt=TabPoint[i] - pVertex->point();
			double distPt = sqrt(DistancePt * DistancePt);
			double wi = 1/0.008/dim/sqrt(2*3.141592) * exp(-(distPt*distPt)/2/0.008/0.008/dim/dim);
			sumWi += wi;
			sumDist += TabDistance[i] * wi;
		}
		
		mean=sumDist / (double)sumWi;
		
		for(int i=0; i<numInSphere; i++)
		{
			Vector_3 DistancePt=TabPoint[i] - pVertex->point();
			double distPt = sqrt(DistancePt * DistancePt);
			double wi = 1/0.008/dim/sqrt(2*3.141592) * exp(-(distPt*distPt)/2/0.008/0.008/dim/dim);

			variance += wi * pow(TabDistance[i] - mean, 2);
		}

		variance = variance / sumWi;
		variance = sqrt(variance);
	}
	
	//variance=sqrt(variance);
	//variance=log(variance+10);
	
	//////////////////////
	/*pVertex->Kh(variance);
	pVertex->Kg(mean);*/
	
	pVertex->curvatureMean = mean;
	pVertex->curvatureVariance = variance;
	
	meanRet = mean;

	return variance;
}

double MSDM::processCovariance(
	Mesh::Vertex* pVertex, double curvature, double curvatureDeg,
	std::vector<double> TabDistance, std::vector<Point_3> TabPoint,
	std::vector<double> TabDistanceDeg, std::vector<Point_3> TabPointDeg, 
	double dim, bool IsGauss)
{
	double Cov1,Cov2;
	double somme1=0;
	double somme2=0;
	int Nb1,Nb2;
	double x1,x2;

	double buff,maxBuff;
	Vector_3 VBuff;

	int IndMax;
	int taille1=TabDistance.size();
	int taille2=TabDistanceDeg.size();
	if(taille1==0 || taille2==0)
		int h=3;

	int ecart=abs(taille1-taille2);

	double sommewi=0;
	for(int i=0;i<TabDistance.size();i++)
	{
		int IndMin=i-2*ecart;
		if(IndMin<0)
			IndMin=0;

		int IndMax=i+2*ecart;
		if(IndMax>TabPointDeg.size()-1)
			IndMax=TabPointDeg.size()-1;

		x1=TabDistance[i];
		x2=TabDistanceDeg[IndMin];
		VBuff=TabPointDeg[IndMin] - TabPoint[i];
		maxBuff=sqrt(VBuff*VBuff);
		
		for(int j=IndMin+1;j<=IndMax;j++)
		{
			VBuff=TabPointDeg[j] - TabPoint[i];
			buff=sqrt(VBuff*VBuff);
			if(buff<maxBuff)
			{		
				maxBuff=buff;
				x2=TabDistanceDeg[j];
			}
		}

		if(!IsGauss)
		{
			somme1+=(x1-curvature)*(x2-curvatureDeg);
			sommewi+=1;
		}
		else
		{
			Vector_3 DistancePt=TabPoint[i] - pVertex->point();
			double distPt=sqrt(DistancePt*DistancePt);
			double wi=1/0.008/dim/sqrt(2*3.141592)*exp(-(distPt*distPt)/2/0.008/0.008/dim/dim);

			sommewi+=wi;
			somme1+=wi*(x1-curvature)*(x2-curvatureDeg);
		}
	}

	Cov1=somme1/sommewi;

	sommewi=0;
	for(int i=0;i<TabDistanceDeg.size();i++)
	{
		int IndMin=i-2*ecart;
		if(IndMin<0)
			IndMin=0;

		int IndMax=i+2*ecart;
		if(IndMax>TabPoint.size()-1)
			IndMax=TabPoint.size()-1;

		x1=TabDistanceDeg[i];
		x2=TabDistance[IndMin];
		VBuff=TabPoint[IndMin]-TabPointDeg[i];
		maxBuff=sqrt(VBuff*VBuff);

		for(int j=IndMin+1;j<=IndMax;j++)
		{
			VBuff=TabPoint[j] - TabPointDeg[i];
			buff=sqrt(VBuff*VBuff);
			if(buff<maxBuff)
			{
				maxBuff=buff;
				x2=TabDistance[j];
			}
		}

		if(!IsGauss)
		{
			somme2+=(x1-curvatureDeg)*(x2-curvature);
			sommewi+=1;
		}
		else
		{
			Vector_3 DistancePt=TabPointDeg[i] - pVertex->point();
			double distPt=sqrt(DistancePt*DistancePt);
			double wi=1/0.008/dim/sqrt(2*3.141592)*exp(-(distPt*distPt)/2/0.008/0.008/dim/dim);

			sommewi+=wi;
			somme2+=wi*(x1-curvatureDeg)*(x2-curvature);
		}
	}

	Cov2=somme2/sommewi;

	if(Cov2>100000 || Cov2>100000)
		int klj=3;

	return (Cov2+Cov1)/2;
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


double DeformationLevel::mean(double mean_p, double mean_q)
{
	return abs(mean_p - mean_q) / max(mean_p, mean_q);
}

double DeformationLevel::deviation(double dev_p, double dev_q)
{
	return abs(dev_p - dev_q) / max(dev_p, dev_q);
}

double DeformationLevel::covariance(double cov_p, double cov_q, double cov_pq)
{
	return abs(cov_p*cov_q - cov_pq) / (cov_p * cov_q);
}


void DeformationLevel::calculateVSIPointDeformationLevel(Mesh* originalMesh, Mesh* watermarkedMesh)
{
	Mesh::Vertex_const_iterator ori_vit = originalMesh->vertices_begin();
	Mesh::Vertex_const_iterator def_vit = watermarkedMesh->vertices_begin();
	Mesh::Vertex_const_iterator ori_vit_end = originalMesh->vertices_end();
	
	assert(originalMesh->size_of_vertices() == watermarkedMesh->size_of_vertices());

	double minDiffVol = FLT_MAX;
	double maxDiffVol = FLT_MIN;
	double avgDiffVol =0;

	double absDiffVol;
	for (; ori_vit != ori_vit_end; ori_vit++, def_vit++)
	{
		absDiffVol = abs( ori_vit->volumeVSI() - def_vit->volumeVSI() );

		avgDiffVol += absDiffVol;
		if (absDiffVol < minDiffVol)
			minDiffVol = absDiffVol;
		if (absDiffVol > maxDiffVol)
			maxDiffVol = maxDiffVol;
	}

	printf("maxDiffVolVSI = %f\n", maxDiffVol);
	printf("minDiffVolVSI = %f\n", minDiffVol);
	printf("avgDiffVolVSI = %f\n", avgDiffVol);
}

void DeformationLevel::calculateSDFPointDeformationLevel(Mesh* originalMesh, Mesh* watermarkedMesh)
{
	Mesh::Vertex_const_iterator ori_vit = originalMesh->vertices_begin();
	Mesh::Vertex_const_iterator def_vit = watermarkedMesh->vertices_begin();
	Mesh::Vertex_const_iterator ori_vit_end = originalMesh->vertices_end();
	
	int ori_num_pts = originalMesh->size_of_vertices();
	int def_num_pts = watermarkedMesh->size_of_vertices();

	assert(ori_num_pts == def_num_pts);

	double minDiffVol = FLT_MAX;
	double maxDiffVol = FLT_MIN;
	double avgDiffVol =0;

	double absDiffVol;
	for (; ori_vit != ori_vit_end; ori_vit++, def_vit++)
	{
		absDiffVol = abs( ori_vit->volumeSDF() - def_vit->volumeSDF() );

		avgDiffVol += absDiffVol;
		if (absDiffVol > maxDiffVol)		maxDiffVol = absDiffVol;
		if (absDiffVol < minDiffVol)		minDiffVol = absDiffVol;
	}

	avgDiffVol = avgDiffVol / ori_num_pts;

	printf("#maxDiffVolSDF = %f\n", maxDiffVol);
	printf("#minDiffVolSDF = %f\n", minDiffVol);
	printf("#avgDiffVolSDF = %f\n", avgDiffVol);
}


void DeformationLevel::calculateNSDFPointDeformationLevel(Mesh* originalMesh, Mesh* watermarkedMesh)
{
	Mesh::Vertex_const_iterator ori_vit = originalMesh->vertices_begin();
	Mesh::Vertex_const_iterator def_vit = watermarkedMesh->vertices_begin();
	Mesh::Vertex_const_iterator ori_vit_end = originalMesh->vertices_end();
	
	int ori_num_pts = originalMesh->size_of_vertices();
	int def_num_pts = watermarkedMesh->size_of_vertices();

	assert(ori_num_pts == def_num_pts);

	double minDiffVol = FLT_MAX;
	double maxDiffVol = FLT_MIN;
	double avgDiffVol =0;

	double absDiffVol;
	for (; ori_vit != ori_vit_end; ori_vit++, def_vit++)
	{
		absDiffVol = abs( ori_vit->volumeNSDF() - def_vit->volumeNSDF() );

		avgDiffVol += absDiffVol;
		if (absDiffVol > maxDiffVol)		maxDiffVol = absDiffVol;
		if (absDiffVol < minDiffVol)		minDiffVol = absDiffVol;
	}

	avgDiffVol = avgDiffVol / ori_num_pts;

	printf("#maxDiffVolSDF = %f\n", maxDiffVol);
	printf("#minDiffVolSDF = %f\n", minDiffVol);
	printf("#avgDiffVolSDF = %f\n", avgDiffVol);
}


double DeformationLevel::calculateGlobalWindowsDeformationLevel(Mesh* originalMesh, Mesh* watermarkedMesh)
{
	return 0;	
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

double Quality::measureSNR(Mesh* originalMesh, Mesh* watermarkedMesh)
{
	double sum = 0.0, differSum = 0.0;

	Mesh::Vertex_const_iterator cvit = originalMesh->vertices_begin();
	Mesh::Vertex_const_iterator cvit_end = originalMesh->vertices_end();
	Mesh::Vertex_const_iterator svit = watermarkedMesh->vertices_begin();
	
	for (; cvit != cvit_end; cvit++, svit++)
	{
		Point_3 coverPoint = cvit->point();
		Point_3 stegoPoint = svit->point();

		sum += pow(coverPoint.x(), 2) + pow(coverPoint.y(), 2) + pow(coverPoint.z(), 2);
		//differSum += pow((stegoPoint.x() - coverPoint.x()), 2) + pow((stegoPoint.y() - coverPoint.y()), 2) + pow((stegoPoint.z() - coverPoint.z()), 2);
		differSum += Utility::pointSquaredDistance(coverPoint, stegoPoint);
	}

	double SNR = 10 * log10(safeDiv(sum, differSum));

	return SNR;
}

double Quality::measurePSNR(Mesh* originalMesh, Mesh* watermarkedMesh)
{
	originalMesh->computeBoundingBox();

	double diagonalLength = originalMesh->diagonalLength();
	int totalVertices = originalMesh->size_of_vertices();

	double differSum = 0.0;

	Mesh::Vertex_const_iterator cvit = originalMesh->vertices_begin();
	Mesh::Vertex_const_iterator cvit_end = originalMesh->vertices_end();
	Mesh::Vertex_const_iterator svit = watermarkedMesh->vertices_begin();
	for(; cvit != cvit_end; cvit++, svit++)
	{
		Point_3 coverPoint = cvit->point();
		Point_3 stegoPoint = svit->point();

		differSum += pow(Utility::pointDistance(coverPoint, stegoPoint), 2);
	}
	
	double RMSE = sqrt(differSum / totalVertices);

	double PSNR = 20 * log10(diagonalLength / RMSE);

	return PSNR;
}

void Quality::TrimString(CString & aString)
{
	aString.TrimLeft();
	aString.TrimRight();
}

void Quality::calculateGeoDistances(
	const QString& originalMeshFileName, const QString& watermarkedMeshFileName, 
	double& mrms, double& mrmsBB,
	double& hausdorff, double& hausdorffwrtBB)
{
	STARTUPINFO si;
	PROCESS_INFORMATION pi;

	SECURITY_ATTRIBUTES sa = {0};
	sa.nLength = sizeof(SECURITY_ATTRIBUTES);
	sa.bInheritHandle = TRUE;

	HANDLE fp = CreateFile(_T("tempfileformetro.txt"),GENERIC_WRITE,0,&sa,CREATE_ALWAYS,FILE_ATTRIBUTE_NORMAL,NULL);

	ZeroMemory(&si,sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi,sizeof(pi));

	si.dwFlags = STARTF_USESHOWWINDOW|STARTF_USESTDHANDLES;
	si.wShowWindow = SW_HIDE;
	si.hStdOutput = fp;
	si.hStdError = fp;

	char* originalMeshName = originalMeshFileName.toAscii().data();
	char* attackedMeshName = watermarkedMeshFileName.toAscii().data();

	CString inFilename = originalMeshName;
	CString outFilename = attackedMeshName;
	CString commandLine = " " + inFilename + " " + outFilename;
	LPTSTR buff = new TCHAR[512];
	_tcscpy(buff,commandLine);
	
	if(!CreateProcess(L"metro.exe", buff, NULL, NULL, TRUE, 0, NULL, NULL, &si ,&pi))
	{
		printf("CreateProcess failed for metro.exe (%d)\n", GetLastError());
	}

	WaitForSingleObject(pi.hProcess,INFINITE);

	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	CloseHandle(fp);

	delete [] buff;
	buff = 0;

	mrms = -1.0;
	FILE * tempMetroFile = fopen("tempfileformetro.txt","r");
	char * pLine = new char[512];
	while(fgets(pLine,512,tempMetroFile))
	{
		CString currentLine = pLine;
		TrimString(currentLine);

		if (currentLine.Left(3)=="RMS")
		{
			int stringLength = currentLine.GetLength();
			int commaPosition = currentLine.Find(':');
			CString stringTemp = currentLine.Mid(commaPosition+2,stringLength-commaPosition-2);
			double rmsTemp = wcstod(stringTemp, NULL);
			if (rmsTemp>mrms)
				mrms = rmsTemp;
		}

		if (currentLine.Left(9)=="Hausdorff")
		{
			int spacePosition1 = currentLine.Find(':') + 1;
			int spacePosition2 = currentLine.Find(' ',spacePosition1+1);
			int leftParenPosition = currentLine.Find('(');
			int spacePosition3 = currentLine.Find(' ',leftParenPosition+1);
			CString stringTemp1 = currentLine.Mid(spacePosition1+1,spacePosition2-spacePosition1-1);
			hausdorff = wcstod(stringTemp1, NULL);
			CString stringTemp2 = currentLine.Mid(leftParenPosition+1,spacePosition3-leftParenPosition-1);
			hausdorffwrtBB = wcstod(stringTemp2, NULL);
		}
	}

	fclose(tempMetroFile);

	delete [] pLine;
	pLine = 0;

	remove("tempfileformetro.txt");

	mrmsBB = hausdorffwrtBB / hausdorff * mrms;
	if (mrms<=0.00000000001)
		mrmsBB = 0.0;

	//fprintf(logFile,"Induced MRMS distance is %f (%f wrt bounding box diagonal)\n",mrms,mrmsBB);
	//fprintf(logFile,"Induced MRMS distance is %f\n",mrms);
	//fprintf(logFile,"Induced Hausdorff distance is %f (%f wrt bounding box diagonal)\n",hausdorff,hausdorffwrtBB);
	//fprintf(logFile,"Induced Hausdorff distance is %f\n",hausdorff);
	//fflush(logFile);

	/*
	// implementation with WinExec() function, but we cannot wait for the termination of the child process...
	tempFileNum++;
	char charIndex[10];
	sprintf(charIndex,"%d",tempFileNum);
	CString stringIndex = charIndex;
	CString stringTempFileName = "tempformetro.txt";
	char charTempFilename[128];
	CStringtochar(stringTempFileName,charTempFilename);
	FILE * tempFile = fopen(charTempFilename,"w");
	CString commandLine = "cmd /c metro rabbit.off rabbit.off > " + stringTempFileName;
	char charcommandLine[128];
	CStringtochar(commandLine,charcommandLine);
	int tempint = WinExec(charcommandLine,SW_HIDE);
	fclose(tempFile);
	//remove(charTempFilename);
	*/
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

double Statistic::calcuateCorrelationCoefficient(const int* seqA, const int sizeA, const int* seqB, const int sizeB)
{
	assert (sizeA == sizeB);

	int totalBits = sizeA;
		
	double meanA = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		meanA += seqA[bitIndex];
	meanA = (meanA / totalBits);
	
	double meanB = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		meanB += seqB[bitIndex];
	meanB = (meanB / totalBits);
	
	double varianceA = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		varianceA += pow((seqA[bitIndex] - meanA), 2);
	
	double varianceB = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		varianceB += pow((seqB[bitIndex] - meanB), 2);
		
	double correlation = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		correlation += (seqA[bitIndex] - meanA) * (seqB[bitIndex] - meanB);
	correlation = (correlation / pow((varianceA * varianceB), 0.5));

	return correlation;
}


double Statistic::calcuateCorrelationCoefficient(const TNT::Array1D<int>& seqA, const TNT::Array1D<int>& seqB)
{
	const int lenSeqA = seqA.dim();
	const int lenSeqB = seqB.dim();

	assert (lenSeqA == lenSeqB);

	int totalBits = lenSeqA;
		
	double meanA = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		meanA += seqA[bitIndex];
	meanA = (meanA / totalBits);
	
	double meanB = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		meanB += seqB[bitIndex];
	meanB = (meanB / totalBits);
	
	double varianceA = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		varianceA += pow((seqA[bitIndex] - meanA), 2);
	
	double varianceB = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		varianceB += pow((seqB[bitIndex] - meanB), 2);
		
	double correlation = 0;
	for (int bitIndex = 0; bitIndex < totalBits; bitIndex++)
		correlation += (seqA[bitIndex] - meanA) * (seqB[bitIndex] - meanB);
	correlation = (correlation / pow((varianceA * varianceB), 0.5));

	return correlation;
}


double Statistic::calculateBitErrorRate(const int* seqA, const int sizeA, const int* seqB, const int sizeB)
{
	assert(sizeA == sizeB);

	int numErrors = 0;
	for (int bitIndex = 0; bitIndex < sizeA; bitIndex++)
	{
		if (seqA[bitIndex] != seqB[bitIndex])
			numErrors++;
	}

	return (double) numErrors / sizeA;
}


double Statistic::calculateBitErrorRate(const TNT::Array1D<int>& seqA, const TNT::Array1D<int>& seqB)
{
	const int lenSeqA = seqA.dim();
	const int lenSeqB = seqB.dim();

	assert (lenSeqA == lenSeqB);

	int numErrors = 0;
	for (int bitIndex = 0; bitIndex < lenSeqA; bitIndex++)
	{
		if (seqA[bitIndex] != seqB[bitIndex])
			numErrors++;
	}

	return (double) numErrors / lenSeqA;
}
