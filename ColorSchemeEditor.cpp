#include "stdafx.h"
#include <QStringList>
#include <QColor>
#include <QInputDialog>
#include <QLayout>
#include <QLabel>
#include <QCheckBox>
#include <QPushButton>

#include "ColorSchemeEditor.h"
#include "ColorSchemeManager.h"

ColorSchemeEditor::ColorSchemeEditor(QWidget* parent, const char* name, Qt::WFlags flags):
QDialog(parent,name,false,flags)
{
	ui.setupUi(this);

	init();

/*	done in setupUI
	connect( ui.buttonOk, SIGNAL( clicked() ), this, SLOT( accept() ) );
	connect( ui.buttonCancel, SIGNAL( clicked() ), this, SLOT( reject() ) );
	connect( ui.newButton, SIGNAL( clicked() ), this, SLOT( OnNewScheme() ) );
	connect( ui.deleteButton, SIGNAL( clicked() ), this, SLOT( OnDeleteScheme() ) );
	connect( ui.renameButton, SIGNAL( clicked() ), this, SLOT( OnRenameScheme() ) );
	connect( ui.csNames, SIGNAL( selected(const QString&) ), this, SLOT( OnSelectScheme(const QString&) ) );
	connect( ui.csNames, SIGNAL( selectionChanged() ), this, SLOT( OnSchemeSelectionChanged() ) );
	*/
}

ColorSchemeEditor::~ColorSchemeEditor()
{

}

void ColorSchemeEditor::loadColorSchemes()
{
	ColorSchemeManager::Load(m_colorSchemes);
}


void ColorSchemeEditor::saveColorSchemes()
{
	ColorSchemeManager::Save(m_colorSchemes);
}


void ColorSchemeEditor::fillNamesList()
{
	ui.csNames->clear();

	colorSchemeHash_t::const_iterator it = m_colorSchemes.begin();
	colorSchemeHash_t::const_iterator it_end = m_colorSchemes.end();

	for (;it != it_end; it++) {
		QString name = it->first.c_str();
		ui.csNames->insertItem(name);
	}
}


ColorScheme* ColorSchemeEditor::getColorSchemeByName( const QString & name )
{
	colorSchemeHash_t::iterator it = m_colorSchemes.find(std::string(name.toAscii()));

	if (it != m_colorSchemes.end()) {
		return it->second;
	} else {
		return NULL;
	}
}


void ColorSchemeEditor::writeDetails( ColorScheme * cs )
{
	//put values in the form
	for (int i=0;i<5;i++) {
		m_colorCheckBoxes[i]->setChecked(i<cs->colors.size());
		if (i<cs->colors.size()) {
			m_colorLabels[i]->OnSetColor(cs->colors[i]);
		} else {
			m_colorLabels[i]->OnSetColor(Qt::white);
		}
	}
}


void ColorSchemeEditor::readDetails( ColorScheme * cs )
{
	//read back values from the form
	cs->colors.clear();
	for (int i=0;i<5;i++) {
		if (m_colorCheckBoxes[i]->isChecked()) {
			cs->colors.push_back(m_colorLabels[i]->color());
		}
	}
}


void ColorSchemeEditor::OnDeleteScheme()
{
	if (m_selectedColorScheme == QString::null) return;

	colorSchemeHash_t::iterator it = m_colorSchemes.find(std::string(m_selectedColorScheme.toAscii()));

	if (it != m_colorSchemes.end()) {
		delete it->second;
		m_colorSchemes.erase(it);
	}

	ui.csNames->removeItem(ui.csNames->currentItem());
}


void ColorSchemeEditor::OnRenameScheme()
{
	if (m_selectedColorScheme == QString::null) return;

	bool ok;
	QString newName = QInputDialog::getText("Rename color scheme", "New name", QLineEdit::Normal, m_selectedColorScheme, &ok, this);

	if (ok) {
		//change name in list
		ui.csNames->changeItem(newName, ui.csNames->currentItem());
		//change in hash map
		colorSchemeHash_t::iterator it = m_colorSchemes.find(std::string(m_selectedColorScheme.toAscii()));

		if (it != m_colorSchemes.end()) {
			ColorScheme* cs = it->second;
			m_colorSchemes.erase(it);
			m_colorSchemes.insert(std::pair<std::string, ColorScheme*>(std::string(newName.toAscii()), cs));
		}

		m_selectedColorScheme = newName;
	}
}


void ColorSchemeEditor::OnNewScheme()
{
	bool ok;
	QString newName = QInputDialog::getText("New color scheme", "Name", QLineEdit::Normal, "<new scheme>", &ok, this);

	if (ok) {
		ColorScheme* cs = new ColorScheme;
		cs->colors.push_back(Qt::black);
		cs->colors.push_back(Qt::white);
		m_colorSchemes.insert(std::pair<std::string, ColorScheme*>(std::string(newName.toAscii()), cs));

		ui.csNames->insertItem(newName);
		ui.csNames->setCurrentItem(ui.csNames->count()-1);
	}
}


void ColorSchemeEditor::OnSelectScheme( const QString & csName )
{
	m_selectedColorScheme = csName;

	writeDetails(getColorSchemeByName(csName));
}


void ColorSchemeEditor::init()
{
	ui.csNames->clear();

	loadColorSchemes();

	fillNamesList();

	//////////////////////////////////////////////////////////////////////////
	//create the place for editing color schemes
	const int rowCount = 5;
	QGridLayout* grid = new QGridLayout(ui.csDetails, rowCount, 4);
	for (int i=0; i<rowCount; i++) {
		QCheckBox* cb = new QCheckBox(ui.csDetails);
		cb->setMaximumSize(QSize(24,24));
		m_colorCheckBoxes[i] = cb;
		connect(cb, SIGNAL(clicked()), SLOT(OnSchemeChanged()));

		ColorLabel* label = new ColorLabel("#" + QString::number(i), ui.csDetails);
		label->setPaletteBackgroundColor(Qt::white);
		label->setMinimumSize(QSize(24,24));
		label->setMaximumSize(QSize(24,24));
		m_colorLabels[i] = label;

		ColorButton* button = new ColorButton("...", Qt::white, ui.csDetails);
		button->setMaximumSize(QSize(24, 24));
		connect(button, SIGNAL(changed()), SLOT(OnSchemeChanged()));

		//changes in button reflect on label
		label->connect(button, SIGNAL(colorSelected(const QColor&)), SLOT(OnSetColor(const QColor&)));

		QSpacerItem* spacer = new QSpacerItem(50, 24);

		grid->addWidget(cb,		i,	0);
		grid->addWidget(label,	i,	1);
		grid->addWidget(button, i,	2);
		grid->addMultiCell(spacer, i,  i, 3, 3);
	}

	ui.csNames->setSelected(0, true);
}


void ColorSchemeEditor::accept()
{
	saveColorSchemes();
	done(QDialog::Accepted);
}


void ColorSchemeEditor::OnSchemeChanged()
{
	//find current color scheme and update it
	if (m_selectedColorScheme == QString::null) return;

	colorSchemeHash_t::iterator it = m_colorSchemes.find(std::string(m_selectedColorScheme.toAscii()));

	if (it != m_colorSchemes.end()) {
		ColorScheme* cs = it->second;
		readDetails(cs);
	}
}


void ColorSchemeEditor::OnSchemeSelectionChanged()
{
	m_selectedColorScheme = ui.csNames->currentText();

	writeDetails(getColorSchemeByName(m_selectedColorScheme));
}
