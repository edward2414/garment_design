//
//  Filename         : config_os_dependant.h
//  Author           : Emmanuel Turquin
//  Purpose          : Things changing from one O.S. to another.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qstring.h>

#ifndef  CONFIG_OS_DEPENDANT_H
# define CONFIG_OS_DEPENDANT_H

namespace config {

  // Directory separators.
# ifdef WIN32
  static const QString DIR_SEP("\\");
  static const QString PATH_SEP(";");
# else
  static const QString DIR_SEP("/");
  static const QString PATH_SEP(":");
# endif // WIN32

} // End of namespace config.

#endif // CONFIG_OS_DEPENDANT_H
