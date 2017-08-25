#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "GeneratedFiles\ui_MainWindow.h"

class QActionGroup;
class FileBrowser;

class MainWindow : public QMainWindow
{
	Q_OBJECT
public:
	MainWindow(QWidget *parent = 0);
	~MainWindow() {}

	Ui::MainWindowClass ui;

	struct CamParams
	{
		CamParams(float _a = 0.0, float _b = 0.0, float _c = 0.0, float _d = 0.0) :a(_a), b(_b), c(_c), d(_d) {}
		float a,b,c,d;
	};
	QVector<CamParams> m_cam;
	FileBrowser* m_browse;

private:
	QActionGroup *selModeActions;
	QActionGroup *displayModeActions;

protected:
	bool event(QEvent* event);
	virtual QMenu* createPopupMenu() { return NULL; }

public slots:
	void selectedMode(QAction* action);

	void ChangeSceneParameters(int tab, float a, float b, float c, float d);

signals:
	void functionKeyPressed(int key, Qt::ButtonState buttonState);
	void arrowKeyPressed(int key, Qt::ButtonState buttonState);

	void selectionMode(int);
};

#endif // MAINWINDOW_H
