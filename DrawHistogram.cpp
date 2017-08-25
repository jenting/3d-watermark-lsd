#include "stdafx.h"

#include <qwt_plot.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_marker.h>
#include <qwt_interval_data.h>
#include <qwt_text.h>

#include "DrawHistogram.h"
#include "Exporter.h"
#include "histogram_item.h"

DrawHistogram::DrawHistogram()
{
    //ctor
}

DrawHistogram::~DrawHistogram()
{
    //dtor
}

void DrawHistogram::drawSDFValue(AppObject* appObject, const bool onVertices, const bool normalized, const int watermarkBits)
{
	Mesh* mesh = appObject->mesh;

	const int totalVertices = appObject->mesh->size_of_vertices();
	const int numBins = 2 * watermarkBits + 2;

    ///create histogram
    QWidget* nw = new QWidget(0, Qt::Window);
    QwtPlot* plot = new QwtPlot(nw);
    plot->setCanvasBackground(QColor(Qt::white));

	if (onVertices)
	{
		if (normalized)
			plot->setTitle("Vertex NLSD Histogram");
		else 
			plot->setTitle("Vertex LSD Histogram");
	}
	else
	{
		if (normalized) 
			plot->setTitle("Facet NLSD Histogram");
		else 
			plot->setTitle("Facet LSD Histogram");
	}

    HistogramItem *histogram = new HistogramItem();
    histogram->setColor(Qt::darkCyan);

    const double interval = (double) 1 / (double) numBins;

    QwtArray<QwtDoubleInterval> intervals(numBins);
    QwtArray<double> counts(numBins);
    for ( int i = 0; i < numBins; i++ )
        counts[i] = 0;

	double maxSDF = 0.0;
	double minSDF = FLT_MAX;
	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
		{
			double volume = (normalized ? vit->volumeNSDF() : vit->volumeSDF());

			int intervalIndex = volume / interval;

			if (volume > maxSDF)
				maxSDF = volume;
			if (volume < minSDF)
				minSDF = volume;

			if (intervalIndex == numBins)
			{
				counts[intervalIndex - 1]++;
				continue;
			}
			counts[intervalIndex]++;
		}
	}
	else
	{
		Mesh::Facet_const_iterator fit = mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = mesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			double volume  = (normalized ? fit->volumeNSDF() : fit->volumeSDF());
			
			int intervalIndex = volume / interval;

			if (volume > maxSDF)
				maxSDF = volume;
			if (volume < minSDF)
				minSDF = volume;

			if(intervalIndex == numBins)
			{
				counts[intervalIndex - 1]++;
				continue;
			}
			counts[intervalIndex]++;
		}
	}
	

    double max = 0;
    for ( int i = 0; i < (int)intervals.size(); i++ )
        if(max < counts[i])
            max = counts[i];

    double pos = 0.0;
    for ( int i = 0; i < (int)intervals.size(); i++ )
	{
        intervals[i] = QwtDoubleInterval(pos, pos + interval);
        pos += interval;
    }

    histogram->setData(QwtIntervalData(intervals, counts));
    histogram->attach(plot);

    plot->setAxisScale(QwtPlot::yLeft, 0.0, max);
    plot->setAxisScale(QwtPlot::xBottom, minSDF, maxSDF, (maxSDF - minSDF) / numBins);
    plot->setAxisTitle(QwtPlot::yLeft, "Number of elements");

	if (onVertices)
	{
		if (normalized)
			plot->setAxisTitle(QwtPlot::xBottom, "Vertex NLSD Value");
		else
			plot->setAxisTitle(QwtPlot::xBottom, "Vertex LSD Value");
	}
	else
	{
		if (normalized)
			plot->setAxisTitle(QwtPlot::xBottom, "Facet NLSD Value");
		else
			plot->setAxisTitle(QwtPlot::xBottom, "Facet LSD Value");
	}

    plot->replot();
    plot->resize(600, 400);
    //plot->print()

	//nw->show();

	QString histogramFileName;
	QString qBins;
	qBins.setNum(numBins, 10);
	if (onVertices)
	{
		if (normalized)
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "sdf_n_v_" + qBins + "_" + appObject->fileName + ".png";
		else
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "sdf_v_" + qBins + "_" + appObject->fileName + ".png";
	}
	else
	{
		if (normalized)
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "sdf_n_f_" + qBins + "_" + appObject->fileName + ".png";
		else
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "sdf_f_" + qBins + "_" + appObject->fileName + ".png";
	}

	QPixmap pixmap(nw->size());
	nw->render(&pixmap);
	pixmap.resize(600, 400);
	pixmap.save(histogramFileName, "PNG", 100);

	cout << "export png file : " << qPrintable(histogramFileName) << endl;
}


///////////////////////////////////////////////////////////////////////////////


void DrawHistogram::drawSDFBins(
	AppObject* appObject,  const bool onVertices, const bool normalize,
	const int watermarkBits, const int logAlpha)
{
	Mesh* mesh = appObject->mesh;

	const int totalVertices = appObject->mesh->size_of_vertices();
	
	const int meanNumBins = watermarkBits + 2;
	const int relationNumBins = 2 * watermarkBits + 2;

	/*
	double minSDF = FLT_MAX;
	double maxSDF = 0.0;
	double minNSDF = FLT_MAX;
	double maxNSDF = 0.0;
	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
		{
			if (vit->volumeSDF() < minSDF)
				minSDF = vit->volumeSDF();
			if (vit->volumeSDF() > maxSDF)
				maxSDF = vit->volumeSDF();
			if (vit->volumeNSDF() < minNSDF)
				minNSDF = vit->volumeNSDF();
			if (vit->volumeNSDF() > maxNSDF)
				maxNSDF = vit->volumeNSDF();
		}
	}
	else
	{
		Mesh::Facet_const_iterator fit = mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = mesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			if (fit->volumeSDF() < minSDF)
				minSDF = fit->volumeSDF();
			if (fit->volumeSDF() > maxSDF)
				maxSDF = fit->volumeSDF();
			if (fit->volumeNSDF() < minNSDF)
				minNSDF = fit->volumeNSDF();
			if (fit->volumeNSDF() > maxNSDF)
				maxNSDF = fit->volumeNSDF();
		}
	}
	*/

	const double minSDF = mesh->getMinVolume();
	const double maxSDF = mesh->getMaxVolume();
	const double avgSDF = mesh->getAvgVolume();
	const double minNSDF = mesh->getMinNVolume();
	const double maxNSDF = mesh->getMaxNVolume();
	const double avgNSDF = mesh->getAvgNVolume();
	const double minStdDevSDF = mesh->getMinStdDevVolume();
	const double maxStdDevSDF = mesh->getMaxStdDevVolume();
	const double avgStdDevSDF = mesh->getAvgStdDevVolume();
	const double minStdDevNSDF = mesh->getMinStdDevNVolume();
	const double maxStdDevNSDF = mesh->getMaxStdDevNVolume();
	const double avgStdDevNSDF = mesh->getAvgStdDevNVolume();

	cout.precision(10);
	cout << "minNSDF = " << minNSDF << " , maxNSDF = " << maxNSDF << " , avgNSDF = " << avgNSDF << endl;
	cout << "minStdDevNSDF = " << minStdDevNSDF << " , maxStdDevNSDF = " << maxStdDevNSDF << " , avgStdDevNSDF = " << avgStdDevNSDF << endl;
	
	//const double range = maxNSDF - minNSDF;
	const double range = maxStdDevNSDF - minStdDevNSDF;
	const double meanInterval = range / meanNumBins;
	const double relationInterval = range / relationNumBins;

	vector<int> meanCounts(meanNumBins, 0);
	vector<int> relationCounts(relationNumBins, 0);

	vector<vector<double> > vertexNSDFBins(meanNumBins);
	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
		{
			//double value = vit->volumeNSDF() - minNSDF;
			double nsdf = vit->volumeNSDF();
			if (nsdf < minStdDevNSDF || nsdf > maxStdDevNSDF)
				continue;

			double value = nsdf - minStdDevNSDF;
			int meanIntervalIndex = value / meanInterval;
			int relationIntervalIndex = value / relationInterval;
			if (meanIntervalIndex == meanNumBins)
			{
				meanCounts[meanIntervalIndex-1] = meanCounts[meanIntervalIndex-1] + 1;
				vertexNSDFBins[meanIntervalIndex-1].push_back(vit->volumeNSDF());
			}
			else
			{
				meanCounts[meanIntervalIndex] = meanCounts[meanIntervalIndex] + 1;
				vertexNSDFBins[meanIntervalIndex].push_back(vit->volumeNSDF());
			}
			
			if (relationIntervalIndex == relationNumBins)
				relationCounts[relationIntervalIndex-1]++;
			else
				relationCounts[relationIntervalIndex]++;
		}
	}
	else
	{
		/*
		Mesh::Facet_const_iterator fit = mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = mesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			double value = (normalize ? (fit->volumeNSDF() - minNSDF) : (fit->volumeSDF() - minSDF));
			int intervalIndex = value / interval;

			if(intervalIndex == numBins)		
			{
				counts[intervalIndex-1] += 1.0;
				vertexNSDFBins[intervalIndex-1].push_back(fit->volumeNSDF());
			}
			else
			{
				counts[intervalIndex] += 1.0;
				vertexNSDFBins[intervalIndex].push_back(fit->volumeNSDF());
			}
		}
		*/
	}
	
	//////////////////////////////////////////////////////////////////////////
	// Check Valid Pair of Bins
	//////////////////////////////////////////////////////////////////////////
	/*
	/// save counts ///
	QString qNumBins, qBinThreshold, qEmbeddingThreshold;
	qNumBins.setNum(numBins, 10);
	qBinThreshold.setNum(binThreshold, 10);
	qEmbeddingThreshold.setNum(embeddingThreshold, 10);

	QString exportFileName;
	if (onVertices)
	{
		QString fileName;
		if (normalize)
			fileName = "counts_vsdf_bins_" + qNumBins + "_binThres_" + qBinThreshold + "_";
		else
			fileName = "counts_vnsdf_bins_" + qNumBins + "_binThres_" + qBinThreshold + "_";

		exportFileName = appObject->fileDir + appObject->fileSeparator + fileName + appObject->fileName + ".txt";
	}
	else
	{
		QString fileName;
		if (normalize)
			fileName = "counts_fsdf_bins_" + qNumBins + "_binThres_" + qBinThreshold + "_";
		else
			fileName = "counts_fnsdf_bins_" + qNumBins + "_binThres_" + qBinThreshold + "_";

		exportFileName = appObject->fileDir + appObject->fileSeparator + fileName + appObject->fileName + ".txt";
	}
	
	QFile file(exportFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);

		for ( int i = 0; i < numBins; i++)
		{
			ts << counts[i];
			if (i != (numBins - 1))
				ts << '\n';
		}

		file.close();
	}
	cout << "export file : " << qPrintable(exportFileName) << endl;
	*/

	TNT::Array1D<int> meanWatermarkSequence = TNT::Array1D<int>(watermarkBits);
	int meanWatermarkIndex = 0;
	for (int i = 1; i < (meanNumBins-1); i++)
	{
		double leftInterval = meanInterval * i;
		double rightInterval = meanInterval * (i+1);
		double medianValue = (leftInterval + rightInterval) / 2;

		vector<double> vertexNSDFBin = vertexNSDFBins[i];
		int vertexNSDFBinSize = vertexNSDFBin.size();
		double sum = 0.0;
		for (int j = 0; j < vertexNSDFBinSize; j++)
			sum += vertexNSDFBin[j];
		
		double avg = sum / vertexNSDFBinSize;

		if (avg >= medianValue)
			meanWatermarkSequence[meanWatermarkIndex] = 1;
		else
			meanWatermarkSequence[meanWatermarkIndex] = 0;

		meanWatermarkIndex++;
	}

	TNT::Array1D<int> relationWatermarkSequence = TNT::Array1D<int>(watermarkBits);
	int relationWatermarkIndex = 0;
	for (int i = 1; i < (relationNumBins-1); i+=2)
	{
		if ((i+1) == (relationNumBins-1))
			break;

		int firstBinIndex = i;
		int secondBinIndex = firstBinIndex + 1;

		int firstBinSize = relationCounts[firstBinIndex];
		int secondBinSize = relationCounts[secondBinIndex];

		if (firstBinSize <= secondBinSize)
			relationWatermarkSequence[relationWatermarkIndex] = 1;
		else
			relationWatermarkSequence[relationWatermarkIndex] = 0;

		relationWatermarkIndex++;
	}

	/*
	cout << "mean : [";
	for (int i = 0; i < watermarkBits; i++)
	{
		cout << meanWatermarkSequence[i];
		if (i != (watermarkBits - 1))
			cout << ",";
	}
	cout << "]\n";
	*/

	/*cout << "relation : [";
	for (int i = 0; i < watermarkBits; i++)
	{
		cout << relationWatermarkSequence[i];
		if (i != (watermarkBits - 1))
			cout << ",";
	}
	cout << "]\n";*/
	

	QString qBins;
	QString qAlpha;
	qBins.setNum(relationNumBins, 10);
	qAlpha.setNum(logAlpha, 10);

	//QString exportMeanFileName = appObject->fileDir + appObject->fileSeparator + "bits_mean" + appObject->fileName + ".txt";
	QString exportRelationFileName = appObject->fileDir + appObject->fileSeparator + "bits_relation_bin_" + qBins + "_" + qAlpha + "_" + appObject->fileName + ".txt";

	ExporterArray exporterArray;
	//exporterArray.exportArray1D(meanWatermarkSequence, exportMeanFileName);
	//exporterArray.exportArray1D(relationWatermarkSequence, exportRelationFileName);

	//////////////////////////////////////////////////////////////////////////
	// Draw Histogram
	//////////////////////////////////////////////////////////////////////////
  
	QwtArray<QwtDoubleInterval> intervals(relationNumBins);
    QwtArray<double> counts(relationNumBins);
    for ( int i = 0; i < relationNumBins; i++)
        counts[i] = relationCounts[i];

	QWidget* nw = new QWidget(0, Qt::Window);
    QwtPlot* plot = new QwtPlot(nw);
    plot->setCanvasBackground(QColor(Qt::white));

	if (onVertices)
		plot->setTitle("Histogram of LSD values (vertex)");
	else 
		plot->setTitle("Histogram of LSD values (facet)");

    HistogramItem *histogram = new HistogramItem();
    histogram->setColor(Qt::darkCyan);

    double max = 0.0;
    for ( int i = 0; i < relationNumBins; i++ )
        if(counts[i] > max)
            max = counts[i];

    double pos = 0.0;
    for ( int i = 0; i < relationNumBins; i++ )
	{
        //intervals[i] = QwtDoubleInterval(pos, pos + interval);
        //pos += interval;
		intervals[i] = QwtDoubleInterval(pos, pos + 1);
        pos += 1;
    }

    histogram->setData(QwtIntervalData(intervals, counts));
    histogram->attach(plot);

    plot->setAxisScale(QwtPlot::yLeft, 0.0, max);
    plot->setAxisScale(QwtPlot::xBottom, 0, relationNumBins, (int) (relationNumBins/10.));
    plot->setAxisTitle(QwtPlot::yLeft, "Total elements");

	if (onVertices)
		plot->setAxisTitle(QwtPlot::xBottom, "Index of each bin");
	else
		plot->setAxisTitle(QwtPlot::xBottom, "Index of each bin");

    plot->replot();
    plot->resize(600, 400);
    //plot->print()

	//nw->show();

	QString histogramFileName;
	if (onVertices)
	{
		if (normalize)
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "vnsdf_bin_" + qBins + "_alpha_" + qAlpha + "_" + appObject->fileName + ".png";
		else
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "sdf_v_bin_" + qBins + "_alpha_" + qAlpha + "_" + appObject->fileName + ".png";
	}
	else
	{
		if (normalize)
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "fnsdf_bin_" +qBins + "_alpha_" + qAlpha + "_" + appObject->fileName + ".png";
		else
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "fsdf_bin_" +qBins + "_alpha_" + qAlpha + "_" + appObject->fileName + ".png";
	}

	QPixmap pixmap(nw->size());
	nw->render(&pixmap);
	pixmap.resize(600, 400);
	pixmap.save(histogramFileName, "PNG", 100);

	cout << "export png file : " << qPrintable(histogramFileName) << endl;
}


///////////////////////////////////////////////////////////////////////////////


void DrawHistogram::drawSDFProbability(AppObject* appObject, const bool onVertices, const bool normalized, const int watermarkBits)
{
	Mesh* mesh = appObject->mesh;

	const int totalVertices = appObject->mesh->size_of_vertices();
	const int numBins = 2 * watermarkBits + 2;

	double maxSDF = 0.0;
	double minSDF = FLT_MAX;
	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
		{
			double volume = (normalized ? vit->volumeNSDF() : vit->volumeSDF());
			if (volume > maxSDF)
				maxSDF = volume;
			if (volume < minSDF)
				minSDF = volume;
		}
	}
	else
	{
		Mesh::Facet_const_iterator fit = mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = mesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			double volume = (normalized ? fit->volumeNSDF() : fit->volumeSDF());
			if (volume > maxSDF)
				maxSDF = volume;
			if (volume < minSDF)
				minSDF = volume;
		}
	}

	double minSDFInt = ((int) floor(minSDF * 100)) * 0.01;
	double maxSDFInt = ((int) ceil(maxSDF * 100)) * 0.01;	
	const double interval = (double) (maxSDFInt - minSDFInt) / (double) (numBins * numBins);

	//////////////////////////////////////////////////////////////////////////
	// Draw Histogram
	//////////////////////////////////////////////////////////////////////////
    QWidget* nw = new QWidget(0, Qt::Window);
    QwtPlot* plot = new QwtPlot(nw);
    plot->setCanvasBackground(QColor(Qt::white));

	if (onVertices)
	{
		if (normalized)
			plot->setTitle("Probability of NLSD (vertex)");
		else 
			plot->setTitle("Probability of LSD (vertex)");
	}
	else
	{
		if (normalized) 
			plot->setTitle("Probability of NLSD (facet)");
		else 
			plot->setTitle("Probability of LSD (facet)");
	}

	HistogramItem *histogram = new HistogramItem();
    histogram->setColor(Qt::darkCyan);

    QwtArray<QwtDoubleInterval> intervals(numBins * numBins);
    QwtArray<double> counts(numBins * numBins);
    for ( int i = 0; i < numBins * numBins; i++ )
        counts[i] = 0;

	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
		{
			double volume = (normalized ? vit->volumeNSDF() : vit->volumeSDF());
			int intervalIndex = volume / interval;
			if (intervalIndex == numBins * numBins)
			{
				counts[intervalIndex - 1]++;
				continue;
			}
			counts[intervalIndex]++;
		}
	}
	else
	{
		Mesh::Facet_const_iterator fit = mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = mesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			double volume = (normalized ? fit->volumeNSDF() : fit->volumeSDF());
			int intervalIndex = volume / interval;
			if(intervalIndex == numBins * numBins)
			{
				counts[intervalIndex - 1]++;
				continue;
			}
			counts[intervalIndex]++;
		}
	}
	
    double max = 0.0;
	double min = FLT_MAX;
	int totalCounts = 0;
    for ( int i = 0; i < numBins * numBins; i++ )
	{
        if (counts[i] > max)
            max = counts[i];
		if (counts[i] < min)
			min = counts[i];
		totalCounts += counts[i];
	}

	double maxProbability = 0.0;
	QwtArray<double> probability(numBins * numBins);
	for (int i = 0; i < numBins * numBins; i++)
	{
		probability[i] = (double) counts[i] / (double) totalCounts;

		if (probability[i] > maxProbability)
			maxProbability = probability[i];
	}

    double pos = 0.0;
    for ( int i = 0; i < numBins * numBins; i++ )
	{
        intervals[i] = QwtDoubleInterval(pos, pos + interval);
        pos += interval;
    }

    histogram->setData(QwtIntervalData(intervals, probability));
    histogram->attach(plot);

	plot->setAxisScale(QwtPlot::yLeft, 0.0, maxProbability);

	if (normalized)
		plot->setAxisScale(QwtPlot::xBottom, minSDFInt, maxSDFInt, 0.1);
	else	
		plot->setAxisScale(QwtPlot::xBottom, minSDFInt, maxSDFInt, 0.1);

    plot->setAxisTitle(QwtPlot::yLeft, "probability");

	if (onVertices)
	{
		if (normalized)
			plot->setAxisTitle(QwtPlot::xBottom, "NLSD Value");
		else
			plot->setAxisTitle(QwtPlot::xBottom, "LSD Value");
	}
	else
	{
		if (normalized)
			plot->setAxisTitle(QwtPlot::xBottom, "NLSD Value");
		else
			plot->setAxisTitle(QwtPlot::xBottom, "LSD Value");
	}

    plot->replot();
    plot->resize(600, 400);
    //plot->print()

	//nw->show();

	QString histogramFileName;
	QString qBins;
	qBins.setNum(numBins, 10);
	if (onVertices)
	{
		if (normalized)
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "sdf_n_v_" + "prob_" + qBins + "_" + appObject->fileName + ".png";
		else
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "sdf_v_" + "prob_" + qBins + "_" + appObject->fileName + ".png";
	}
	else
	{
		if (normalized)
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "sdf_n_f_" + "prob_" + qBins + "_" + appObject->fileName + ".png";
		else
			histogramFileName = appObject->fileDir + appObject->fileSeparator + "sdf_f_" + "prob_" + qBins + "_" + appObject->fileName + ".png";
	}

	QPixmap pixmap(nw->size());
	nw->render(&pixmap);
	pixmap.resize(600, 400);
	pixmap.save(histogramFileName, "PNG", 100);

	cout << "export png file : " << qPrintable(histogramFileName) << endl;
}