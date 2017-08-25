#pragma once

#include <QDialog>
#include "ui_glLights.h"
#include "gl_info_structs.h"

class ConfigureGlLights : public QDialog {
	Q_OBJECT

private:
	Ui::ConfigureGlLights ui;

	GL_Light m_lights[8];

public:
	ConfigureGlLights(QWidget* parent = NULL, const char* name = NULL, Qt::WFlags flags = 0);
	~ConfigureGlLights();

private:
	void SaveConfiguration(int num);
	void LoadConfiguration(int num);

	void init();

public slots:
	void OnSelectLight(int lightNum );
	void OnSelectSpecular();
	void OnSelectDiffuse();
	void OnApplyLightSettings();	
	void LoadConfiguration();	
	void SaveConfiguration();
signals:

};