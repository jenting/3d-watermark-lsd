#pragma once

#include <QDialog>
#include <QCheckBox>
#include <QString>

#include "ColorButton.h"
#include "GeneratedFiles\ui_ColorSchemeEditor.h"

class ColorSchemeEditor : public QDialog {
	Q_OBJECT

private:
	Ui::ColorSchemeEditor ui;

	ColorLabel* m_colorLabels[5];
	stdext::hash_map<std::string,ColorScheme*> m_colorSchemes;
	QString m_selectedColorScheme;
	QCheckBox* m_colorCheckBoxes[5];

protected:
	void loadColorSchemes();
	void saveColorSchemes();
	void fillNamesList();
	ColorScheme* getColorSchemeByName( const QString & name );
	void writeDetails( ColorScheme * cs );
	void readDetails( ColorScheme * cs );

	void init();

public:
	ColorSchemeEditor(QWidget* parent = NULL, const char* name = NULL, Qt::WFlags flags = 0);
	~ColorSchemeEditor();

public slots:
	void OnDeleteScheme();
	void OnRenameScheme();
	void OnNewScheme();
	void OnSelectScheme(const QString& csName);
	void accept();
	void OnSchemeChanged();
	void OnSchemeSelectionChanged();

};
