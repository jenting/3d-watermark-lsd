#include "stdafx.h"

#include "AppManager.h"
#include "AppParams.h"

#include <QTextStream>
#include <QFileDialog>

 void AppManager::OnAppParams()
{
	int normalizeSdf;
	switch(m_appParameters.lsd_normalize)
	{
		case None: normalizeSdf = 0; break;
		case MinMax: normalizeSdf = 1; break;
		case Log: normalizeSdf = 2; break;
	}

	int searchSpace;
	switch(m_appParameters.watermark_searchSpace)
	{
		case MIN_DIST: searchSpace = 0; break;
		case AVG_DIST: searchSpace = 1; break;
		case MAX_DIST: searchSpace = 2; break;
	}

	m_appParamsDialog->ui.debugmode->setChecked(m_appParameters.debugMode);
	m_appParamsDialog->ui.multithreaded->setChecked(m_appParameters.multithreaded);
	m_appParamsDialog->ui.dirRunner->setText(m_appParameters.dirRunner);
	m_appParamsDialog->ui.fileName1->setText(m_appParameters.fileName1);
	m_appParamsDialog->ui.fileName2->setText(m_appParameters.fileName2);
	m_appParamsDialog->ui.threadsNum->setValue(m_appParameters.threadsNum);
	
	m_appParamsDialog->ui.lsd_gridSize->setValue(m_appParameters.lsd_gridSize);
	m_appParamsDialog->ui.lsd_numberOfCones->setValue(m_appParameters.lsd_numberOfCones);
	m_appParamsDialog->ui.lsd_coneSeparation->setValue(m_appParameters.lsd_coneSeparation);
	m_appParamsDialog->ui.lsd_raysInCone->setValue(m_appParameters.lsd_raysInCone);
	m_appParamsDialog->ui.lsd_removeSameDirectionIntersection->setChecked(m_appParameters.lsd_removeSameDirectionIntersection);
	m_appParamsDialog->ui.lsd_gaussianWeights->setChecked(m_appParameters.lsd_gaussianWeights);
	m_appParamsDialog->ui.lsd_removeOutlier->setChecked(m_appParameters.lsd_removeOutlier);
	m_appParamsDialog->ui.lsd_normalizeSDF->setCurrentItem(normalizeSdf);
	m_appParamsDialog->ui.lsd_log_alpha->setValue(m_appParameters.lsd_log_alpha);
	m_appParamsDialog->ui.lsd_smoothing->setChecked(m_appParameters.lsd_smoothing);
	m_appParamsDialog->ui.lsd_smoothingIterations->setValue(m_appParameters.lsd_smoothing_iterations);

	m_appParamsDialog->ui.watermark_seed->setValue(m_appParameters.watermark_seed);
	m_appParamsDialog->ui.watermark_bits->setValue(m_appParameters.watermark_bits);
	m_appParamsDialog->ui.robustness_factor->setValue(m_appParameters.robustness_factor);
	m_appParamsDialog->ui.watermark_searchSpace->setCurrentItem(searchSpace);
	m_appParamsDialog->ui.watermark_maxCycles->setValue(m_appParameters.watermark_maxCycles);
	m_appParamsDialog->ui.watermark_particleNumbers->setValue(m_appParameters.watermark_particleNumbers);

	m_appParamsDialog->ui.normsRev->setChecked(m_appParameters.normalsReversed);

	m_appParamsDialog->show();
}

void AppManager::appParamsOk()
{
	appParamsApply();
	m_appParamsDialog->hide();
}

void AppManager::appParamsApply()
{
	switch(m_appParamsDialog->ui.lsd_normalizeSDF->currentItem())
	{
		case 0: m_appParameters.lsd_normalize = None; break;
		case 1: m_appParameters.lsd_normalize = MinMax; break;	
		case 2: m_appParameters.lsd_normalize = Log; break;
	}

	switch(m_appParamsDialog->ui.watermark_searchSpace->currentItem())
	{
		case 0: m_appParameters.watermark_searchSpace = MIN_DIST; break;
		case 1: m_appParameters.watermark_searchSpace = AVG_DIST; break;
		case 2: m_appParameters.watermark_searchSpace = MAX_DIST; break;
	}

	m_appParameters.debugMode = m_appParamsDialog->ui.debugmode->isChecked();
	m_appParameters.multithreaded = m_appParamsDialog->ui.multithreaded->isChecked();
	m_appParameters.dirRunner = m_appParamsDialog->ui.dirRunner->text();
	m_appParameters.fileName1 = m_appParamsDialog->ui.fileName1->text();
	m_appParameters.fileName2 = m_appParamsDialog->ui.fileName2->text();
	m_appParameters.threadsNum = m_appParamsDialog->ui.threadsNum->value();

	m_appParameters.lsd_gridSize = m_appParamsDialog->ui.lsd_gridSize->value();
	m_appParameters.lsd_numberOfCones = m_appParamsDialog->ui.lsd_numberOfCones->value();
	m_appParameters.lsd_coneSeparation = m_appParamsDialog->ui.lsd_coneSeparation->value();
	m_appParameters.lsd_raysInCone = m_appParamsDialog->ui.lsd_raysInCone->value();
	m_appParameters.lsd_removeSameDirectionIntersection = m_appParamsDialog->ui.lsd_removeSameDirectionIntersection->isChecked();
	m_appParameters.lsd_gaussianWeights = m_appParamsDialog->ui.lsd_gaussianWeights->isChecked();
	m_appParameters.lsd_removeOutlier = m_appParamsDialog->ui.lsd_removeOutlier->isChecked();
	m_appParameters.lsd_log_alpha = m_appParamsDialog->ui.lsd_log_alpha->value();
	m_appParameters.lsd_smoothing = m_appParamsDialog->ui.lsd_smoothing->isChecked();
	m_appParameters.lsd_smoothing_iterations = m_appParamsDialog->ui.lsd_smoothingIterations->value();

	m_appParameters.watermark_seed = m_appParamsDialog->ui.watermark_seed->value();
	m_appParameters.watermark_bits = m_appParamsDialog->ui.watermark_bits->value();
	m_appParameters.robustness_factor = m_appParamsDialog->ui.robustness_factor->value();
	m_appParameters.watermark_maxCycles = m_appParamsDialog->ui.watermark_maxCycles->value();
	m_appParameters.watermark_particleNumbers = m_appParamsDialog->ui.watermark_particleNumbers->value();

	m_appParameters.normalsReversed = m_appParamsDialog->ui.normsRev->isChecked();

	saveApplicationParameters();
}

void AppManager::saveApplicationParameters()
{
	ProgSettings ps;
	ps.beginGroup("AppParams");

	ps.writeEntry("debugMode", m_appParameters.debugMode);
	ps.writeEntry("multithreaded", m_appParameters.multithreaded);
	ps.writeEntry("threadsNum", m_appParameters.threadsNum);
	ps.writeEntry("directory_runner", m_appParameters.dirRunner);
	ps.writeEntry("file_name_1", m_appParameters.fileName1);
	ps.writeEntry("file_name_2", m_appParameters.fileName2);
	ps.writeEntry("norms_rev", m_appParameters.normalsReversed);

	ps.writeEntry("lsd_gridSize", m_appParameters.lsd_gridSize);
	ps.writeEntry("lsd_numberOfCones", m_appParameters.lsd_numberOfCones);
	ps.writeEntry("lsd_coneSeparation", m_appParameters.lsd_coneSeparation);
	ps.writeEntry("lsd_raysInCone", m_appParameters.lsd_raysInCone);
	ps.writeEntry("lsd_removeSameDirectionIntersection", m_appParameters.lsd_removeSameDirectionIntersection);
	ps.writeEntry("lsd_gaussianWeights", m_appParameters.lsd_gaussianWeights);
	ps.writeEntry("lsd_removeOutlier", m_appParameters.lsd_removeOutlier);
	ps.writeEntry("lsd_normalize", m_appParameters.lsd_normalize);
	ps.writeEntry("lsd_log_alpha", m_appParameters.lsd_log_alpha);
	ps.writeEntry("lsd_smoothing", m_appParameters.lsd_smoothing);
	ps.writeEntry("lsd_smoothingIterations", m_appParameters.lsd_smoothing_iterations);

	ps.writeEntry("watermark_seed", m_appParameters.watermark_seed);
	ps.writeEntry("watermark_bits", m_appParameters.watermark_bits);
	ps.writeEntry("robustness_factor", m_appParameters.robustness_factor);
	ps.writeEntry("watermark_searchSpace", m_appParameters.watermark_searchSpace);
	ps.writeEntry("watermark_maxCycles", m_appParameters.watermark_maxCycles);
	ps.writeEntry("watermark_particleNumbers", m_appParameters.watermark_particleNumbers);

	ps.endGroup();
}

void AppManager::loadApplicationParameters()
{
	ProgSettings ps;
	ps.beginGroup("AppParams");
	switch(ps.readNumEntry("lsd_normalize", 0))
	{
		case 0: m_appParameters.lsd_normalize = None; break;
		case 1: m_appParameters.lsd_normalize = MinMax; break;	
		case 2: m_appParameters.lsd_normalize = Log; break;
	}

	switch(ps.readNumEntry("watermark_searchSpace", 0))
	{
		case 0: m_appParameters.watermark_searchSpace = MIN_DIST; break;
		case 1: m_appParameters.watermark_searchSpace = AVG_DIST; break;	
		case 2: m_appParameters.watermark_searchSpace = MAX_DIST; break;
	}

	m_appParameters.debugMode = ps.readBoolEntry("debugMode", false);
	m_appParameters.multithreaded = ps.readBoolEntry("multithreaded", true);
	m_appParameters.threadsNum = ps.readNumEntry("threadsNum", 8);
	m_appParameters.dirRunner = ps.readEntry("directory_runner", "");
	m_appParameters.fileName1 = ps.readEntry("file_name_1", "");
	m_appParameters.fileName2 = ps.readEntry("file_name_2", "");
	m_appParameters.normalsReversed = ps.readBoolEntry("norms_rev", false);

	m_appParameters.lsd_gridSize = ps.readNumEntry("lsd_gridSize", 40);
	m_appParameters.lsd_numberOfCones = ps.readNumEntry("lsd_numberOfCones", 4);
	m_appParameters.lsd_coneSeparation = ps.readNumEntry("lsd_coneSeparation", 5);
	m_appParameters.lsd_raysInCone = ps.readNumEntry("lsd_raysInCone", 36);
	m_appParameters.lsd_removeSameDirectionIntersection = ps.readBoolEntry("lsd_removeSameDirectionIntersection", true);
	m_appParameters.lsd_gaussianWeights = ps.readBoolEntry("lsd_gaussianWeights", false);
	m_appParameters.lsd_removeOutlier = ps.readBoolEntry("lsd_removeOutlier", false);
	m_appParameters.lsd_log_alpha = ps.readNumEntry("lsd_log_alpha", 4);
	m_appParameters.lsd_smoothing = ps.readBoolEntry("lsd_smoothing", false);
	m_appParameters.lsd_smoothing_iterations = ps.readNumEntry("lsd_smoothingIterations", 1);

	m_appParameters.watermark_seed = ps.readNumEntry("watermark_seed", 123456789);
	m_appParameters.watermark_bits = ps.readNumEntry("watermark_bits", 16);
	m_appParameters.robustness_factor = ps.readNumEntry("robustness_factor", 3);
	m_appParameters.watermark_maxCycles = ps.readNumEntry("watermark_maxCycles", 10);
	m_appParameters.watermark_particleNumbers = ps.readNumEntry("watermark_particleNumbers", 50);

	ps.endGroup();
}