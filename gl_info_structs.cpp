#include "stdafx.h"

#include "gl_info_structs.h"

#include <qstringlist.h>

GLenum getLightEnum(const int lightNum)
{
	switch (lightNum) {
		case 0: return GL_LIGHT0;
		case 1: return GL_LIGHT1;
		case 2: return GL_LIGHT2;
		case 3: return GL_LIGHT3;
		case 4: return GL_LIGHT4;
		case 5: return GL_LIGHT5;
		case 6: return GL_LIGHT6;
		case 7: return GL_LIGHT7;
	}
	return GL_LIGHT0;
}

int convertFloatToInt(const GLfloat fvalue)
{
	int result = (float) fvalue * 255;
	return result;
}

void GL_Light::LoadFromGL(const int lightNum)
{
	GLfloat diffuse_array[4];
	GLfloat specular_array[4];
	GLfloat ambient_array[4];

	GLfloat position[4];

	GLenum lightId = getLightEnum(lightNum);

	//Color
	glGetLightfv(lightId, GL_DIFFUSE,  diffuse_array);
	glGetLightfv(lightId, GL_SPECULAR, specular_array);
	glGetLightfv(lightId, GL_AMBIENT, ambient_array);

	QColor diffuseColor(convertFloatToInt(diffuse_array[0]),
		convertFloatToInt(diffuse_array[1]),
		convertFloatToInt(diffuse_array[2]));

	QColor specularColor(convertFloatToInt(specular_array[0]),
		convertFloatToInt(specular_array[1]),
		convertFloatToInt(specular_array[2]));

	QColor ambientColor(convertFloatToInt(ambient_array[0]),
		convertFloatToInt(ambient_array[1]),
		convertFloatToInt(ambient_array[2]));


	specular = specularColor.rgb();
	diffuse = diffuseColor.rgb();
	ambient = ambientColor.rgb();

	//Position
	glGetLightfv(lightId, GL_POSITION, position);
	posx = position[0];
	posy = position[1];
	posz = position[2];
	
	//Enabled or not
	GLboolean lightEnabled = glIsEnabled(lightId);	
	enabled = lightEnabled;	
}

void GL_Light::ApplyToGL(const int lightNum)
{
	GLenum lightId = getLightEnum(lightNum);

	if (enabled)
		glEnable(lightId);
	else
		glDisable(lightId);

	QColor specularColor = QColor(specular);
	QColor diffuseColor = QColor(diffuse);
	QColor ambientColor = QColor(ambient);

	GLfloat fspecular[4];
	fspecular[0] = (float) specularColor.red() / 255;
	fspecular[1] = (float) specularColor.green() / 255;
	fspecular[2] = (float) specularColor.blue() / 255;
	fspecular[3] = 0.0;
	glLightfv(lightId, GL_SPECULAR, fspecular);
	
	GLfloat fdiffuse[4];
	fdiffuse[0] = (float) diffuseColor.red() / 255;
	fdiffuse[1] = (float) diffuseColor.green() / 255;
	fdiffuse[2] = (float) diffuseColor.blue() / 255;
	fdiffuse[3] = 0.0;
	glLightfv(lightId, GL_DIFFUSE, fdiffuse);

	GLfloat fambient[4];
	fambient[0] = (float) ambientColor.red() / 255;
	fambient[1] = (float) ambientColor.green() / 255;
	fambient[2] = (float) ambientColor.blue() / 255;
	fambient[3] = 0.0;
	glLightfv(lightId, GL_AMBIENT, fambient);

	
	GLfloat position[4];
	position[0] = (float) posx;
	position[1] = (float) posy;
	position[2] = (float) posz;
	position[3] = 1.0;
	glLightfv(lightId, GL_POSITION, position);
}

void GL_Light::LoadFromString(const QString lightConfig)
{
	QStringList sl = QStringList::split(";", lightConfig);

	assert(sl.size()==6);

	enabled = (bool)sl.first().toInt; sl.pop_front();
	diffuse = sl.first().toUInt(); sl.pop_front();
	specular = sl.first().toUInt(); sl.pop_front();
	ambient = QColor(0, 0, 0).rgb();

	posx = sl.first().toFloat(); sl.pop_front();
	posy = sl.first().toFloat(); sl.pop_front();
	posz = sl.first().toFloat();
}

QString GL_Light::ToString()
{
	QString result;
	result.sprintf("%d;%lu;%lu;%f;%f;%f",enabled, diffuse, specular, posx, posy, posz);

	return result;
}