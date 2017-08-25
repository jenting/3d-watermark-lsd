#ifndef __RENDERING_PARAMS_H_
#define __RENDERING_PARAMS_H_

#include <qcolor.h>
#include <vector>
/**
*	structure holding a series of colors, if one then solid, if 2 or 3 colors then 
*  it is a gradient color scheme
*/
struct ColorScheme {
	/** vector of selected colors in scheme */
	std::vector<QColor> colors;
	unsigned char alpha;

	QColor def() const {
		if (colors.size()>0)
			return colors[0];
		else
			return Qt::lightGray;
	}

	QColor gradient(const float val) const;	
};

/**
 *	Rendering parameters which may be shared between objects or be one per object
 */
struct RenderingParams {
	/** perform smooth shading */
	bool m_smoothShading;
	/** perform culling */
	bool m_culling;
	/** what polygon mode to draw in (full, lines...) */
	int m_polygonMode;		
	/** perform anti aliasing */
	bool m_antialiasing;
	/** draw visible vertices */
	bool m_superimposeVertices;
	/** draw visible edges */
	bool m_superimposeEdges;
	/** use normals when rendering */
	bool m_useNormals;
	/** use lighting model */
	bool m_lighting;

	/** in non-smooth painters how to do normals */
	bool m_smoothNormals;
	
	/** color for facets */
	QColor m_facetColor;
	/** color for edges */
	QColor m_edgeColor;
	/** color for vertices */
	QColor m_vertexColor;
	/** color scheme used by painters */
	QString m_colorSchemeName;

	/** what painter to use for facets */
	int m_renderModeFacets;
	/** what painter to use for vertices */
	int m_renderModeEdges;

	int m_debugPainter;

	bool m_overlay;

	unsigned char m_alpha;

	/** thickness of points */
	int m_vertexThickness;
	/** thickness of lines */
	int m_edgeThickness;

	bool m_vboSupported;
};

#endif