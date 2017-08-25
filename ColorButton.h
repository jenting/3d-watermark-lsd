#ifndef __COLOR_BUTTON_H_
#define __COLOR_BUTTON_H_

#include <qobject.h>
#include <qpushbutton.h>
#include <qcolor.h>
#include <qcolordialog.h>
#include <qlabel.h>
#include <qstring.h>

class ColorButton : public QPushButton {
	Q_OBJECT
private:
	QColor m_color;
protected:

public:
	ColorButton(const QString& caption, const QColor& initialColor, QWidget* parent, const char* name = 0) 
		:m_color(initialColor), QPushButton(caption, parent, name)
	{
		connect(this,SIGNAL(clicked()), SLOT(OnClicked()));	
	}

public slots:

	/** show color selection dialog */
	void OnClicked() {
		QColor newColor = QColorDialog::getColor(m_color, (QWidget*)parent());
		if (newColor.isValid()) {
			m_color = newColor;
			emit colorSelected(m_color);
			emit changed();
		}
	}

signals:
	void colorSelected(const QColor& Color);
	void changed();
};

class ColorLabel : public QLabel {
	Q_OBJECT
private:

public:
	ColorLabel(const QString& text, QWidget* parent, const char* name = 0)
		: QLabel(text, parent, name)
	{
	}

	const QColor& color() const { return paletteBackgroundColor();}

public slots:
	void OnSetColor(const QColor& color) {
		setPaletteBackgroundColor(color);
	}
};

#endif