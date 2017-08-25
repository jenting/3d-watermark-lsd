#include "stdafx.h"
#include "BasicRenderingParams.h"

#include <QColorDialog>
#include <QCheckBox>
#include <QComboBox>
#include <Q3ListBox>
#include <QSpinBox>

#include "ColorSchemeManager.h"
#include "ColorSchemeEditor.h"

RenderingParamsDialog::RenderingParamsDialog(QWidget* parent, const char* name, Qt::WFlags flags):
QDialog(parent,name,false,flags)
{
	ui.setupUi(this);
	// signals and slots connections - in setupUI
	init();
}

RenderingParamsDialog::~RenderingParamsDialog()
{

}

void RenderingParamsDialog::changeMeshColor()
{	
	m_meshColor = QColorDialog::getColor(m_meshColor, this);
	ui.meshColorLabel->setPaletteBackgroundColor(m_meshColor);
}

void RenderingParamsDialog::changeVertexColor()
{
	m_vertexColor = QColorDialog::getColor(m_vertexColor, this);
	ui.vertexColorLabel->setPaletteBackgroundColor(m_vertexColor);
}

void RenderingParamsDialog::changeEdgeColor()
{
	m_edgeColor = QColorDialog::getColor(m_edgeColor, this);
	ui.edgeColorLabel->setPaletteBackgroundColor(m_edgeColor);
}


void RenderingParamsDialog::init()
{       
	ui.polygonMode->insertItem("Point");
	ui.polygonMode->insertItem("Line");
	ui.polygonMode->insertItem("Fill");

	ui.edgePainters->insertItem("None");        

	ui.vertexPainters->insertItem("None");    

	ui.colorSchemes->insertStringList(ColorSchemeManager::GetNames());		
}


void RenderingParamsDialog::OnAlphaLevelChanged( int val )
{
	//alphaLevelText->setText(QString::number(((float)val/255)*100, 'g', 1) + "%");
}


void RenderingParamsDialog::SetAlphaLevel( int val )
{
	//alphaLevelText->setText(QString::number(val) + "%");
	ui.alphaLevel->setValue(val);
}

void RenderingParamsDialog::onFillFacetPainters( const QStringList & painterNames )
{
	ui.facetPainters->clear();

	for (int i=0; i<painterNames.size(); i++) {
		QString name = painterNames[i];
		ui.facetPainters->insertItem(name);
	}
}


void RenderingParamsDialog::OnRenderingParams( RenderingParams * renderingParams )
{
	//set pointer to rendering params and parse them to text box
	m_renderingParams = renderingParams;

	//fill fields
	ReadRenderingParams();
}


void RenderingParamsDialog::OnApply()
{
	//update rendering params structure
	WriteRenderingParams();

	//emit signal that rendering params have changed

	emit renderingParamsChanged(m_renderingParams);
}


void RenderingParamsDialog::ReadRenderingParams()
{
	ui.cbSmoothShading->setChecked(m_renderingParams->m_smoothShading);
	ui.cbCulling->setChecked(m_renderingParams->m_culling);
	ui.cbImposeVertices->setChecked(m_renderingParams->m_superimposeVertices);
	ui.cbImposeEdges->setChecked(m_renderingParams->m_superimposeEdges);
	ui.cbLighting->setChecked(m_renderingParams->m_lighting);
	ui.cbAntialiasing->setChecked(m_renderingParams->m_antialiasing);
	ui.cbUseNormals->setChecked(m_renderingParams->m_useNormals);
	ui.cbSmoothNormals->setChecked(m_renderingParams->m_smoothNormals);
	ui.cbOverlay->setChecked(m_renderingParams->m_overlay);

	m_meshColor = m_renderingParams->m_facetColor;
	ui.meshColorLabel->setPaletteBackgroundColor(m_meshColor);
	m_edgeColor = m_renderingParams->m_edgeColor;
	ui.edgeColorLabel->setPaletteBackgroundColor(m_edgeColor);
	m_vertexColor = m_renderingParams->m_vertexColor;
	ui.vertexColorLabel->setPaletteBackgroundColor(m_vertexColor);	

	ui.alphaLevel->setValue(m_renderingParams->m_alpha);

	int itemIndex = 0;
	switch(m_renderingParams->m_polygonMode)
	{
		case 0x1B00 :
			itemIndex = 0;
			break;
		case 0x1B01 :
			itemIndex = 1;
			break;
		case 0x1B02 :
			itemIndex = 2;
			break;
	}
	ui.polygonMode->setCurrentItem(itemIndex);

	ui.facetPainters->setCurrentItem(m_renderingParams->m_renderModeFacets);
	ui.edgePainters->setCurrentItem(m_renderingParams->m_renderModeEdges);
	ui.vertexPainters->setCurrentItem(0);

	ui.debugPainters->setCurrentItem(m_renderingParams->m_debugPainter);

	ui.vertexThickness->setValue(m_renderingParams->m_vertexThickness);
	ui.edgeThickness->setValue(m_renderingParams->m_edgeThickness);

	int qlbi = ui.colorSchemes->findText(m_renderingParams->m_colorSchemeName);
	if (qlbi >= 0) {
		ui.colorSchemes->setCurrentItem(qlbi);		
	}

}


void RenderingParamsDialog::WriteRenderingParams()
{
	m_renderingParams->m_smoothShading = ui.cbSmoothShading->isChecked();
	m_renderingParams->m_culling = ui.cbCulling->isChecked();
	m_renderingParams->m_superimposeVertices = ui.cbImposeVertices->isChecked();
	m_renderingParams->m_superimposeEdges = ui.cbImposeEdges->isChecked();
	m_renderingParams->m_lighting = ui.cbLighting->isChecked();
	m_renderingParams->m_antialiasing = ui.cbAntialiasing->isChecked();
	m_renderingParams->m_useNormals = ui.cbUseNormals->isChecked();
	m_renderingParams->m_smoothNormals = ui.cbSmoothNormals->isChecked();
	m_renderingParams->m_overlay = ui.cbOverlay->isChecked();

	m_renderingParams->m_facetColor = m_meshColor;	
	m_renderingParams->m_edgeColor = m_edgeColor;
	m_renderingParams->m_vertexColor = m_vertexColor;

	m_renderingParams->m_alpha = (unsigned char)ui.alphaLevel->value();

	switch(ui.polygonMode->currentItem())
	{
		case 0:
			m_renderingParams->m_polygonMode = 0x1B00;
			break;
		case 1:
			m_renderingParams->m_polygonMode = 0x1B01;
			break;
		case 2:
			m_renderingParams->m_polygonMode = 0x1B02;
			break;
	}

	m_renderingParams->m_renderModeFacets = ui.facetPainters->currentItem();
	m_renderingParams->m_renderModeEdges = ui.edgePainters->currentItem();

	m_renderingParams->m_debugPainter = ui.debugPainters->currentItem();

	m_renderingParams->m_vertexThickness = ui.vertexThickness->value();
	m_renderingParams->m_edgeThickness = ui.edgeThickness->value();	

	m_renderingParams->m_colorSchemeName = ui.colorSchemes->currentText();
}


void RenderingParamsDialog::OnEditColorSchemes()
{
	ColorSchemeEditor* cse = new ColorSchemeEditor(this);

	if (cse->exec() == QDialog::Accepted) {
		ui.colorSchemes->clear();
		ui.colorSchemes->insertStringList(ColorSchemeManager::GetNames());		
	}
}


void RenderingParamsDialog::OnFillDebugPainters( const QStringList & painterNames )
{
	ui.debugPainters->clear();

	for (int i=0; i<painterNames.size(); i++) {
		QString name = painterNames[i];
		ui.debugPainters->insertItem(name);
	}
}
