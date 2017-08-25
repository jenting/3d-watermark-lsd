#include "stdafx.h"
#include "AppParams.h"

#include <deque>
#include <QToolTip>

AppParams::AppParams(QWidget* parent, const char* name, Qt::WFlags flags):
QDialog(parent, name, false, flags)
{
	ui.setupUi(this);
	init();
}

AppParams::~AppParams()
{

}

void AppParams::init()
{
	ui.lsd_normalizeSDF->insertItem("None");
	ui.lsd_normalizeSDF->insertItem("Min-Max");
	ui.lsd_normalizeSDF->insertItem("Log");

	ui.watermark_searchSpace->insertItem("Min Dist");
	ui.watermark_searchSpace->insertItem("Avg Dist");
	ui.watermark_searchSpace->insertItem("Max Dist");
}


void AppParams::on_applyBot_clicked()
{
	emit pressdApply();
}

void AppParams::on_buttonOk_clicked()
{
	emit pressdOk();
}