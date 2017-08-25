#ifndef DIRRUNNER_H
#define DIRRUNNER_H

#include "stdafx.h"

#include "AppObject.h"

#include <QWidget>

//#include "GeneratedFiles\ui_DirRunner.h"

/*namespace Ui {
    class DirRunner;
}*/

class DirRunner : public QWidget
{
    Q_OBJECT
public:
    //explicit DirRunner(QWidget *parent = 0);
	DirRunner(QWidget *parent = 0);
	~DirRunner() {};
	
	// deformation analysis //
	void dirRunnerAvgDDegree(AppObject* appObject, const QString& fileDir);
	void dirRunnerMaxDDegree(AppObject* appObject, const QString& fileDir);

	// attack //
	void dirRunnerNoiseAttack(const QString& fileDir);
	void dirRunnerQuantizationAttack(const QString& fileDir);
	void dirRunnerReorderingAttack(const QString& fileDir);
	void dirRunnerSimilarityTransformAttack(const QString& fileDir);
	void dirRunnerSimplificationAttack(const QString& fileDir);
	void dirRunnerSmoothingAttack(const QString& fileDir);
	void dirRunnerSubdivisionAttack(const QString& fileDir);

private:
    //Ui::dirRunner *ui;

public slots:

signals:

};

#endif // DIRRUNNER_H
