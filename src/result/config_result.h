//
//  Filename         : config_result.h
//  Author           : Emmanuel Turquin
//  Purpose          : Information about the Result.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  CONFIG_RESULT_H
# define CONFIG_RESULT_H

# include <qstring.h>

namespace config {

  static const QString	RESULT_CONFIG_FILE("result.xml");
  static const bool	RESULT_DISP_MODEL = true;
  static const bool	RESULT_DISP_SURFACE = true;
  static const bool	RESULT_DISP_CONTOUR = false;
  static const bool	RESULT_DISP_TYPE = false;

} // End of namespace config.

#endif // CONFIG_RESULT_H
