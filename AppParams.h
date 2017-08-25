#ifndef __APPPARAMS_H_INC__
#define __APPPARAMS_H_INC__

#include <QDialog>
#include "GeneratedFiles\ui_AppParams.h"

class AppParams : public QDialog
{
	Q_OBJECT

public:
	AppParams(QWidget* parent = NULL, const char* name = NULL, Qt::WFlags flags = 0);
	~AppParams();

	Ui::AppParams ui;

signals:
	void pressdOk();
	void pressdApply();

private:
	void init();

private slots:
	void on_buttonOk_clicked();
	void on_applyBot_clicked();
};

#endif // __APPPARAMS_H_INC__
