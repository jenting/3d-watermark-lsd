#ifndef FILEBROWSER_H
#define FILEBROWSER_H

#include <QDialog>
#include <QTreeView>

class MyTreeView : public QTreeView
{
	Q_OBJECT
public:
	MyTreeView(QWidget *parent = 0) : QTreeView(parent) {}
	virtual void selectionChanged(const QItemSelection & selected, const QItemSelection & deselected)
	{
		emit itemSelectionChanged();
	}
signals:
	void itemSelectionChanged();
};


#include "GeneratedFiles\ui_FileBrowser.h"

class QDirModel;



class FileBrowser : public QWidget
{
	Q_OBJECT

public:
	FileBrowser(QWidget *parent = 0);
	~FileBrowser() {}

private:
	Ui::FileBrowserClass ui;
	QDirModel *m_treemodel;
	QDirModel *m_listmodel;

protected:
	virtual void resizeEvent(QResizeEvent *event);

private slots:
	void changedDir();
	void selectFile(const QModelIndex& index);

signals:
	void openFile(QString name);

};

#endif // FILEBROWSER_H
