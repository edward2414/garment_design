//
//  Filename         : config_app_info.h
//  Author           : Emmanuel Turquin
//  Purpose          : Information about the application.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  CONFIG_APP_INFO_H
# define CONFIG_APP_INFO_H

# include <qstring.h>

namespace config {

  static const QString	APP_NAME(QMAKE_APP_NAME);
  static const QString	APP_VERSION(QMAKE_APP_VERSION);
  static const QString	APP_ABOUT
    (
     "<H2>" + APP_NAME + " " + APP_VERSION + "</H2>"
     "<P>A sketching interface for garment design.</P>"
     "<UL>"
     "<LI>Emmanuel Turquin"
     "<LI>Marie-Paule Cani"
     "<LI>John Hughes"
     "</UL>"
     "<B>(C) 2004 Emmanuel Turquin & Evasion.</B>"
     );

  static const QString	APP_INTERFACE_FILE("interface.xml");

} // End of namespace config.

#endif // CONFIG_APP_INFO_H
