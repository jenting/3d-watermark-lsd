#include <QDialog>

#include "GeneratedFiles\ui_BasicRenderingParams.h"

#include "RenderingParams.h"

class RenderingParamsDialog : public QDialog {
	Q_OBJECT

private:
	Ui::RenderingParamsDialog ui;

	QColor m_edgeColor;
	QColor m_vertexColor;
	QColor m_meshColor;
	ColorScheme m_colorScheme;

	RenderingParams* m_renderingParams;

private:
	virtual void SetAlphaLevel( int val );
	virtual void ReadRenderingParams();
	virtual void WriteRenderingParams();

public:
	RenderingParamsDialog(QWidget* parent = NULL, const char* name = NULL, Qt::WFlags flags = 0);
	~RenderingParamsDialog();

public slots:
	virtual void changeMeshColor();
	virtual void changeVertexColor();
	virtual void changeEdgeColor();
	virtual void OnAlphaLevelChanged( int val );
	virtual void onFillFacetPainters( const QStringList & painterNames );
	virtual void OnRenderingParams( RenderingParams * renderingParams );
	virtual void OnApply();
	virtual void OnFillDebugPainters( const QStringList & painterNames );

protected slots:
	virtual void init();
	virtual void OnEditColorSchemes();

signals:
	void renderingParamsChanged(RenderingParams* renderingParams);
};
