#ifndef _COLOR_SCHEME_MANAGER_H_
#define _COLOR_SCHEME_MANAGER_H_

#include "RenderingParams.h"
#include <qstring.h>
#include <qstringlist.h>
#include <qcolor.h>

#include <hash_map>

typedef stdext::hash_map<std::string,ColorScheme*> colorSchemeHash_t;

/**
 *	Utility class to load/save color schemes from registry
 */
class ColorSchemeManager {
public:
	static QStringList GetNames() {
		ProgSettings ps;
		return ps.entryList("ColorSchemes");
	}

	static ColorScheme GetScheme(const QString& name) {
		ProgSettings ps;
		ps.beginGroup("ColorSchemes");
		bool entryok = true;
		QString csValue;
		if (name != "default")
			csValue = ps.readEntry(name, QString::null, &entryok);
		else // default
			csValue = "4278522367;4278255594;4280745728;4294967040;4294901760;";

		if (entryok) {
			//fill the gradient in the list
			QStringList csColors = QStringList::split(";", csValue);

			ColorScheme cs;
			cs.colors.reserve(csColors.size());

			int numOfColors = csColors.size();
			for (int i=0;i<numOfColors;i++) {
				unsigned int cint = csColors.first().toUInt();
				QColor c = cint;
				cs.colors.push_back(c);
				if (i<numOfColors-1) 
					csColors.pop_front();
			}

			return cs;
		} else {
			return ColorScheme();
		}
	}

	static void Load(colorSchemeHash_t& colorSchemes) {
		ProgSettings ps;

		QStringList csList = ps.entryList("ColorSchemes");

		ps.beginGroup("ColorSchemes");

		QStringList::const_iterator it = csList.begin();
		bool entryok;
		while (it != csList.end()) {
			QString csName = *it;
			QString csValue = ps.readEntry(csName, QString::null, &entryok);
			if (!entryok) continue;

			//fill the gradient in the list
			QStringList csColors = QStringList::split(";", csValue);

			ColorScheme* cs = new ColorScheme;
			cs->colors.reserve(csColors.size());

			int numOfColors = csColors.size();
			for (int i=0;i<numOfColors;i++) {
				unsigned int cint = csColors.first().toUInt();
				QColor c = cint;
				cs->colors.push_back(c);
				if (i<numOfColors-1) 
					csColors.pop_front();
			}

			colorSchemes.insert(std::pair<std::string,ColorScheme*>(std::string(csName.toAscii()), cs));

			it++;
		}

		ps.endGroup();
	}

	static void Save(const colorSchemeHash_t& colorSchemes) {
		ProgSettings ps;
		ps.beginGroup("ColorSchemes");

		//save color schemes to registry
		colorSchemeHash_t::const_iterator it = colorSchemes.begin();
		colorSchemeHash_t::const_iterator it_end = colorSchemes.end();

		for (;it != it_end; it++) {
			ColorScheme* cs = it->second;
			QString colorValues;
			for (int i=0; i<cs->colors.size();i++)
				colorValues += QString::number(cs->colors[i].rgb()) + ";";
			ps.writeEntry(it->first.c_str(), colorValues);
		}

		ps.endGroup();
	}
};

#endif