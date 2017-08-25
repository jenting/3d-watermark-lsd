#pragma once

#include <QDialog>
#include <QPixmap>
#include "GeneratedFiles\ui_glLights.h"

#include "gl_info_structs.h"

class ConfigureGlLights : public QDialog
{
	Q_OBJECT
public:
	ConfigureGlLights(QWidget* parent, Qt::WFlags flags = 0);
	~ConfigureGlLights();

	void init();
	void applyLight();

public slots:
	virtual void OnSelectLight( int lightNum );
	virtual void OnSelectSpecular();
	virtual void OnSelectDiffuse();
	virtual void OnApplyLightSettings();
	virtual void LoadConfiguration();
	virtual void LoadConfiguration( int num );
	virtual void SaveConfiguration( int num );
	virtual void SaveConfiguration();

signals:
	void lightChanged();

protected:
	GL_Light m_lights[8];

private:
	Ui::ConfigureGlLights ui;

	QPixmap image0;

private slots:
	void on_ambient_color_clicked();
};
