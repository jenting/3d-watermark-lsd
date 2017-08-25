#include "stdafx.h"
#include "glLights.h"

#include <QColorDialog>


ConfigureGlLights::ConfigureGlLights(QWidget* parent, Qt::WFlags flags /* = 0 */)
: QDialog(parent, flags) 
{
	ui.setupUi(this);

/*	done in setupUi
	connect( ui.selected_light, SIGNAL( activated(int) ), this, SLOT( OnSelectLight(int) ) );
	connect( ui.ok_button, SIGNAL( clicked() ), this, SLOT( close() ) );
	connect( ui.diffuse_color, SIGNAL( clicked() ), this, SLOT( OnSelectDiffuse() ) );
	connect( ui.specular_color, SIGNAL( clicked() ), this, SLOT( OnSelectSpecular() ) );
	connect( ui.apply_button, SIGNAL( clicked() ), this, SLOT( OnApplyLightSettings() ) );
	connect( ui.load_config, SIGNAL( clicked() ), this, SLOT( LoadConfiguration() ) );
	connect( ui.save_config, SIGNAL( clicked() ), this, SLOT( SaveConfiguration() ) );
	connect( ui.saved_configs, SIGNAL( selected(int) ), this, SLOT( SaveConfiguration(int) ) );
*/
	init();
}

ConfigureGlLights::~ConfigureGlLights()
{

}

void ConfigureGlLights::OnSelectLight( int lightNum )
{
	//GL_Light light;
	//light.LoadFromGL(lightNum);

	//color
	ui.diffuse_color_label->setPaletteBackgroundColor(
		QColor(m_lights[lightNum].diffuse));
	ui.specular_color_label->setPaletteBackgroundColor(
		QColor(m_lights[lightNum].specular));
	ui.ambient_color_label->setPaletteBackgroundColor(
		QColor(m_lights[lightNum].ambient));

	//Position
	ui.posx->setText(QString::number(m_lights[lightNum].posx));
	ui.posy->setText(QString::number(m_lights[lightNum].posy));
	ui.posz->setText(QString::number(m_lights[lightNum].posz));

	//Enabled or not
	ui.light_enabled->setChecked(m_lights[lightNum].enabled);
}


void ConfigureGlLights::OnSelectSpecular()
{
	QColor sc = QColorDialog::getColor(ui.specular_color_label->paletteBackgroundColor(), this, "Select specular color");
	ui.specular_color_label->setPaletteBackgroundColor(sc);
}


void ConfigureGlLights::OnSelectDiffuse()
{
	QColor dc = QColorDialog::getColor(ui.diffuse_color_label->paletteBackgroundColor(), this, "Select diffuse color");
	ui.diffuse_color_label->setPaletteBackgroundColor(dc);
}


void ConfigureGlLights::on_ambient_color_clicked()
{
	QColor dc = QColorDialog::getColor(ui.ambient_color_label->paletteBackgroundColor(), this, "Select ambient color");
	ui.ambient_color_label->setPaletteBackgroundColor(dc);

}


void ConfigureGlLights::OnApplyLightSettings()
{
	applyLight();
	emit lightChanged();
}

void ConfigureGlLights::applyLight()
{
	int lightNum = ui.selected_light->currentIndex();

	m_lights[lightNum].Set(ui.light_enabled->isChecked(),
		ui.posx->text().toFloat(),ui.posy->text().toFloat(),ui.posz->text().toFloat(),
		ui.specular_color_label->paletteBackgroundColor().rgb(),
		ui.diffuse_color_label->paletteBackgroundColor().rgb(),
		ui.ambient_color_label->paletteBackgroundColor().rgb());

	m_lights[lightNum].ApplyToGL(lightNum);
}


void ConfigureGlLights::LoadConfiguration( int num )
{
	//work with registry to load configuration
	ProgSettings qs;	

	qs.beginGroup(QString("light_config_")+QString::number(num));
	
	if ((num == 0) && (qs.childKeys().size() == 0)) // empty default
	{ // create a default
		qs.writeEntry("light_0", "1;4285493103;4285624689;10.000000;1.000000;20.000000");
		qs.writeEntry("light_1", "1;4278190080;4278190080;0.000000;0.000000;0.000000");
	}

	bool configExists;	

	for (int i=0;i<8;i++) 
	{
		bool thisExist = false;
		QString lightString = qs.readEntry("light_"+QString::number(i), "", &thisExist);
		if (thisExist)
		{
			configExists = true;
			m_lights[i].LoadFromString(lightString);
			m_lights[i].ApplyToGL(i);
		}
	}

	if (configExists) {

		//now load values to current selected light
		int lightNum = ui.selected_light->currentItem();
		ui.light_enabled->setChecked(m_lights[lightNum].enabled);
		ui.posx->setText(QString::number(m_lights[lightNum].posx));
		ui.posy->setText(QString::number(m_lights[lightNum].posy));
		ui.posz->setText(QString::number(m_lights[lightNum].posz));
		ui.specular_color_label->setPaletteBackgroundColor(m_lights[lightNum].specular);
		ui.diffuse_color_label->setPaletteBackgroundColor(m_lights[lightNum].diffuse);
	} else {
		QMessageBox::warning(g_main,
			"Config does not exist",
			"Configuration does not exist and cannot be loaded");
	}

	qs.endGroup();
}


void ConfigureGlLights::LoadConfiguration()
{
	int selectedConfig = ui.saved_configs->currentItem();

	if (selectedConfig>=0)
		LoadConfiguration(selectedConfig);
	emit lightChanged();
}


void ConfigureGlLights::SaveConfiguration( int num )
{
	//first save current settings if they changed...
	OnApplyLightSettings();

	//work with registry to save configuration
	ProgSettings qs;

	qs.beginGroup(QString("light_config_")+QString::number(num));

	for (int i=0;i<8;i++) {		
		qs.writeEntry("light_"+QString::number(i), m_lights[i].ToString());
	}

	qs.endGroup();
}


void ConfigureGlLights::SaveConfiguration()
{
	int selectedConfig = ui.saved_configs->currentItem();

	if (selectedConfig>=0)
		SaveConfiguration(selectedConfig);
}


void ConfigureGlLights::init()
{
	OnSelectLight(0);
}

