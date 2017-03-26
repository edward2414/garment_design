//
//  Filename         : config_canvas.h
//  Author           : Emmanuel Turquin
//  Purpose          : Information about the Canvas.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  CONFIG_CANVAS_H
# define CONFIG_CANVAS_H

# include <qstring.h>

namespace config {

  static const QString	CANVAS_CONFIG_FILE("canvas.xml");
  static const int	CANVAS_SAMPLING_INT(5);
  static const int	CANVAS_SNAPPING_INT(10);
  static const float	CANVAS_SPLITTING_THRESHOLD(0.5);
  static const bool	CANVAS_CHAIN_PTS(true);
  static const bool	CANVAS_STROKE_PTS(true);
  static const bool	CANVAS_TEXTURE(true);
  static const bool	CANVAS_BBOX(true);
  static const bool	CANVAS_TYPE(false);
  static const bool	CANVAS_ANTIALIASING(true);
  static const int	GRID_WIDTH(20);
  static const int	GRID_HEIGHT(20);
  static const int	CANVAS_GM_X(128);
  static const int	CANVAS_GM_Y(128);
  static const float	CANVAS_BORDER_X(1.0);
  static const float	CANVAS_BORDER_Y(1.0);

} // End of namespace config.

#endif // CONFIG_CANVAS_H
