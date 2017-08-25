#include "stdafx.h"
#include "RenderingParams.h"

#include <assert.h>

QColor ColorScheme::gradient(const float val) const
{
	if(val < 0 || val > 1)
		printf("X \"val < 0 or val > 1\"\n");

	assert(val>=0 && val<=1.0);

	if (colors.size()>0)
	{
		float intervalSize = (float) 1 / colors.size();
		int colorIndex = (float) val / intervalSize;
		if (colorIndex >= colors.size() - 1)
			return colors[colors.size()-1];
		else {
			QColor c1 = colors[colorIndex];
			QColor c2 = colors[colorIndex+1];

			float perc = (val - intervalSize * colorIndex) / intervalSize;

			int red = ((float) c1.red() * (1-perc) + (float) c2.red() * perc);
			int green = ((float) c1.green() * (1-perc) + (float) c2.green() * perc);
			int blue = ((float) c1.blue() * (1-perc) + (float) c2.blue() * perc);

			/*qDebug("Gradient value %f (%f\%) between [%d,%d,%d] and [%d,%d,%d] = [%d,%d,%d]",
				val, perc, c1.red(), c1.green(), c1.blue(),
				c2.red(), c2.green(), c2.blue(), red, green, blue);*/

			return QColor(red, green, blue);
		}		
	}
	else
	{
		return Qt::lightGray;
	}
}