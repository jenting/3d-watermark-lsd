#ifndef DRAWHISTOGRAM_H
#define DRAWHISTOGRAM_H

#include "stdafx.h"
#include "AppObject.h"

class DrawHistogram
{
	private:

	protected:

    public:
        DrawHistogram();
        virtual ~DrawHistogram();

		///  minSDF/interval ///
		void drawSDFValue(AppObject* appObject, const bool onVertices, const bool normalized, const int watermarkBits);

		///  (sdf - minSDF)/interval ///
		void drawSDFBins(
			AppObject* appObject, const bool onVertices, const bool normalize,
			const int watermarkBits, const int logAlpha);

		///  ///
		void drawSDFProbability(AppObject* appObject, const bool onVertices, const bool normalized, const int watermarkBits);
};

#endif // DRAWHISTOGRAM_H
