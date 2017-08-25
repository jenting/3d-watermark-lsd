#pragma once

#include <QGLViewer/qglviewer.h>

class PartSignature;

class SdfViewer : public QGLViewer
{
Q_OBJECT

public slots:
	void ChangeSceneParameters(float centerx,float centery, float centerz, float radius);

	void OnDraw();

	void OnShowWorldInfo();
	void saveManualSnap() { saveSnapshot(false, false); }

public :
	SdfViewer(QWidget *parent, PartSignature* part = NULL, const QGLWidget* shareWith = NULL);
	void setPart(PartSignature* part) { m_part = part; }
  
protected :
	virtual void resizeGL(int width, int height);

	virtual void draw();
	virtual void drawWithNames();
	virtual void init();
	virtual QString helpString() const;

	virtual void select(const QMouseEvent* event);

	virtual void mousePressEvent(QMouseEvent *event);
	virtual void mouseReleaseEvent(QMouseEvent *event);

signals:

	void drawWithNamesNeeded();
	void mydrawNeeded();
	void mydrawNeeded(PartSignature* psig);
	void selected(const unsigned int i, Qt::KeyboardModifiers modif);
	void windowResized(int width, int height);

	void showContextMenu(QPoint pos);

private:
	PartSignature *m_part;
	QPoint m_pressPos;
};

