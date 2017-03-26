//
//  Filename         : config_pretreatment.h
//  Author           : Emmanuel Turquin
//  Purpose          : Information about the Pretreatment.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  CONFIG_PRETREATMENT_H
# define CONFIG_PRETREATMENT_H

# include <qstring.h>

namespace config {

  static const QString	PRETREATMENT_CONFIG_FILE("pretreatment.xml");
  static const bool	PRETREATMENT_CAMERA_PERSPECTIVE(false);
  static const bool	PRETREATMENT_DISP_AXIS(true);
  static const bool	PRETREATMENT_DISP_BBOX(true);
  static const bool	PRETREATMENT_DISP_MODEL(true);
  static const bool	PRETREATMENT_DISP_OCTREE(false);
  static const bool	PRETREATMENT_DISP_DF(true);
  static const int	PRETREATMENT_OCTREE_LEVEL(5);
  static const double	PRETREATMENT_BBOX_X(0);
  static const double	PRETREATMENT_BBOX_Y(0);
  static const double	PRETREATMENT_BBOX_Z(0);
  static const double	PRETREATMENT_BBOX_X_STEP(1);
  static const double	PRETREATMENT_BBOX_Y_STEP(1);
  static const double	PRETREATMENT_BBOX_Z_STEP(1);
  static const unsigned	PRETREATMENT_TEXTURE_SIZE(512);
  static const unsigned	PRETREATMENT_TEXTURE_MAX_SIZE(1024);
  static const unsigned	PRETREATMENT_DF_X(32);
  static const unsigned	PRETREATMENT_DF_Y(32);
  static const unsigned	PRETREATMENT_DF_Z(32);

} // End of namespace config.

#endif // CONFIG_PRETREATMENT_H
