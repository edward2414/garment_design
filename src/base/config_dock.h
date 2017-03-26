//
//  Filename         : config_dock.h
//  Author           : Emmanuel Turquin
//  Purpose          : Information about the Dock windows.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  CONFIG_DOCK_H
# define CONFIG_DOCK_H

# include <qstring.h>

namespace config {

  static const QString	DOCK_DEFAULT_CAPTION("DockWindow");
  static const int	DOCK_MIN_WIDTH(180);
  static const int	DOCK_MIN_HEIGHT(120);

} // End of namespace config.

#endif // CONFIG_DOCK_H
