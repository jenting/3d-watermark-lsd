#ifndef __GL_INFO_STRUCTS_H_
#define __GL_INFO_STRUCTS_H_

#include <qstring.h>
#include <qcolor.h>

struct GL_Light {
	bool enabled;
	float posx,posy,posz;
	unsigned int specular;
	unsigned int diffuse;
	unsigned int ambient;

	GL_Light() :
	enabled(false), posx(0), posy(0), posz(0), specular(0), diffuse(0), ambient(0)
	{
	}

	GL_Light(
		const bool ienabled,
		const float iposx, 
		const float iposy,
		const float iposz,
		const unsigned int ispecular,
		const unsigned int idiffuse,
		const unsigned int iambient) : 
	enabled(ienabled), 
		posx(iposx), posy(iposy), posz(iposz),
		specular(ispecular), diffuse(idiffuse), ambient(iambient)
	{}

	void Set(const bool ienabled,
		const float iposx, 
		const float iposy,
		const float iposz,
		const unsigned int ispecular,
		const unsigned int idiffuse,
		const unsigned int iambient)
	{
		enabled = ienabled;
		posx = iposx;
		posy = iposy;
		posz = iposz;
		diffuse = idiffuse;
		specular = ispecular;
		ambient = iambient;
	}

	void LoadFromGL(const int lightNum);
	void ApplyToGL(const int lightNum);

	void LoadFromString(const QString lightConfig);
	QString ToString();
};

#endif