#include "MainWindow.h"

#include <QKeyEvent>
#include <QTabBar>
#include <QDockWidget>

#include "FileBrowser.h"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	selModeActions = new QActionGroup(this);
	ui.selectObjectAction->setData(1);
	ui.selectObjectAction->setChecked(true);
	selModeActions->addAction(ui.selectObjectAction);
	ui.selectFacetAction->setData(2);
	selModeActions->addAction(ui.selectFacetAction);

	connect(selModeActions, SIGNAL( selected(QAction*) ), this, SLOT( selectedMode(QAction*) ) );

	m_browse = NULL;

	m_browse = new FileBrowser(parent);
	QDockWidget *dock = new QDockWidget("File Browser", this);
	dock->setWidget(m_browse);
	dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

	addDockWidget(Qt::LeftDockWidgetArea, dock); 
	dock->hide();
	//dock->show();

	ui.menuView->addAction(ui.mainToolBar->toggleViewAction());
	ui.menuView->addAction(dock->toggleViewAction());
}

/*
void MainWindow::tabChanged()
{
	int curTab = ui.tabBar->currentIndex();
	if (curTab < m_cam.size())
	{
		CamParams &c = m_cam[curTab];
		ui.world->ChangeSceneParameters(c.a, c.b, c.c, c.d);
	}
	emit changedTab(curTab);
}
*/

void MainWindow::ChangeSceneParameters(int tab, float a, float b, float c, float d)
{
	if (m_cam.size() < tab + 1)
		m_cam.resize(tab + 1);
	m_cam[tab] = CamParams(a, b, c, d);

	ui.world->ChangeSceneParameters(a, b, c, d);
}


bool MainWindow::event(QEvent* e)
{
	//QMessageBox::about(this, "Shape Diameter Function Project", "SDF Project by Lior Shapira v1.0");
	if (e->type() == QEvent::KeyPress) 
	{
		QKeyEvent* ke = (QKeyEvent*)e;
		if (ke->key() >= Qt::Key_0 && ke->key() < Qt::Key_9)
		{
			emit functionKeyPressed(ke->key() - Qt::Key_0, ke->state()); 
			return true;
		}
		switch (ke->key())
		{
		case Qt::Key_Left: case Qt::Key_Up: case Qt::Key_Right: case Qt::Key_Down:
			emit arrowKeyPressed(ke->key(), ke->state());
			break;

		}
	}

	return QMainWindow::event(e);
}

void MainWindow::selectedMode(QAction* action)
{
	emit selectionMode(action->data().toInt());
}


/*
void MainWindow::on_xBot_clicked()
{
	int ind = ui.tabBar->currentIndex();
	if (ind == 0)
		return;
	ui.tabBar->removeTab(ind);
	emit tabRemoved(ind);
}
*/
