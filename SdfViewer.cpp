#include "stdafx.h"

#include "SdfViewer.h" 
#include <QGLViewer/Vec.h>
#include <math.h>

#include "GlobalLog.h"

// Constructor must call the base class constructor.
SdfViewer::SdfViewer(QWidget *parent, PartSignature* part, const QGLWidget* shareWith)
: QGLViewer(parent, shareWith, Qt::Widget), m_part(part)
{	
	restoreStateFromFile();
  //help();

	connect(this, SIGNAL(drawNeeded()), this, SLOT(OnDraw()));
}

void SdfViewer::ChangeSceneParameters(float centerx,float centery, float centerz, float radius)
{
	//this->camera()->setZNearCoefficient(0.001);

	setSceneCenter(qglviewer::Vec(centerx, centery, centerz));
	setSceneRadius(radius);

	showEntireScene();

	/** change viewer position and direction */
	/*qglviewer::Vec a(1.31, 0, 0);
	qglviewer::Vec b(-1.0, 0, 0);
	camera()->setPosition(a);
	camera()->setViewDirection(b);*/

	qglviewer::Vec p = camera()->position();
	qglviewer::Vec d = camera()->viewDirection();

	QString temp;
	temp.sprintf("Camera position: [%f %f %f]", p.x, p.y, p.z);
	SDFLOG4(temp);
	temp.sprintf("Camera direction: [%f %f %f]", d.x, d.y, d.z);
	SDFLOG4(temp);
	temp.sprintf("Near plane at %f, Far plane at %f", camera()->zNear(), camera()->zFar());
	SDFLOG4(temp);

	/*printf("Camera: position [%.2f %.2f %.2f] direction [%.2f %.2f %.2f]\nNear plane %.2f, Far plane %.2f\nScene radius %.2f around [%.2f %.2f %.2f]",
		camera()->position().x, camera()->position().y, camera()->position().z,
		camera()->viewDirection().x, camera()->viewDirection().y, camera()->viewDirection().z,
		camera()->zNear(), camera()->zFar(),
		camera()->sceneRadius(), camera()->sceneCenter().x, camera()->sceneCenter().y, camera()->sceneCenter().z);*/
}

void SdfViewer::draw()
{
	QGLViewer::draw();
}

void SdfViewer::drawWithNames()
{
	emit drawWithNamesNeeded();
}

void SdfViewer::OnDraw()
{
	makeCurrent();
	if (m_part == NULL)
		emit mydrawNeeded();
	else
		emit mydrawNeeded(m_part);
}

void SdfViewer::init()
{
	restoreStateFromFile();

	ProgSettings ps;
	ps.beginGroup("default");
	QString snapshotdir = ps.readEntry("snapshotdir", "e:/three_dimension/watermarking/tmp/snapshot"); //E:\three_dimension\watermarking\tmp\snapshot
	ps.endGroup();

	setSnapshotFileName(snapshotdir.toAscii());
	setSnapshotFormat("PNG");

	setShortcut(SAVE_SCREENSHOT, Qt::CTRL + Qt::Key_S);
	setShortcut(DISPLAY_FPS, Qt::CTRL + Qt::Key_Z);

	setMouseBinding(Qt::LeftButton, FRAME, ROTATE);
	//setMouseBinding(Qt::LeftButton + Qt::CTRL, CAMERA, ROTATE);
	setMouseBinding(Qt::RightButton, CAMERA, TRANSLATE);
	setMouseBinding(Qt::LeftButton, SELECT, true);

	setWheelBinding(Qt::NoModifier, CAMERA, MOVE_FORWARD);

	/*
#if QT_VERSION < 0x040000
	setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::AltButton);
	setHandlerKeyboardModifiers(QGLViewer::FRAME,  Qt::NoButton);
	setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::ControlButton);
#else
	setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::AltModifier);
	setHandlerKeyboardModifiers(QGLViewer::FRAME,  Qt::NoModifier);
	setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::ControlModifier);
#endif
	*/

	//setAxisIsDrawn();

	setBackgroundColor(Qt::white);
	setForegroundColor(Qt::black);

	emit viewerInitialized();
}

QString SdfViewer::helpString() const
{
  QString text("<h2>I n t e r f a c e</h2>");
  text += "A GUI can be added to a QGLViewer widget using Qt's <i>designer</i>. ";
  text += "All the available QGLViewer signals and slots are listed in a <b>qglviewer.cw</b> file, ";
  text += "located in the QGLViewer <i>include</i> directory.";
  return text;
}

// overrided
void SdfViewer::select(const QMouseEvent* event)
{
	beginSelection(event->pos());
	drawWithNames();
	endSelection(event->pos());
	postSelection(event->pos());
	emit selected((const unsigned int)selectedName(), event->modifiers());
}

void SdfViewer::resizeGL(int width, int height)
{
	QGLViewer::resizeGL(width, height);
	emit windowResized(width, height);
}

void SdfViewer::OnShowWorldInfo()
{
	QString temp;
	temp.sprintf("Camera: position [%.2f %.2f %.2f] direction [%.2f %.2f %.2f]\nNear plane %.2f, Far plane %.2f\nScene radius %.2f around [%.2f %.2f %.2f]",
		camera()->position().x, camera()->position().y, camera()->position().z,
		camera()->viewDirection().x, camera()->viewDirection().y, camera()->viewDirection().z,
		camera()->zNear(), camera()->zFar(),
		camera()->sceneRadius(), camera()->sceneCenter().x, camera()->sceneCenter().y, camera()->sceneCenter().z);

	QMessageBox::information(g_main, "World Information", temp);
}


void SdfViewer::mousePressEvent(QMouseEvent *event)
{
	if (event->button() == Qt::RightButton)
	{
		m_pressPos = event->globalPos();
	}
	QGLViewer::mousePressEvent(event);
}

void SdfViewer::mouseReleaseEvent(QMouseEvent *event) 
{
	if (event->button() == Qt::RightButton)
	{
		QPoint mv = m_pressPos - event->globalPos();
		int mvs = qAbs(mv.x()) + qAbs(mv.y());
		if (mvs < 4)
		{
			select(event);
			updateGL();
			emit showContextMenu(event->globalPos());
		}
	}

	QGLViewer::mouseReleaseEvent(event);
}