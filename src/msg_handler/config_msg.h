//
//  Filename         : config_msg.h
//  Author           : Emmanuel Turquin
//  Purpose          : Information about the MsgHandler class.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  CONFIG_MSG_H
# define CONFIG_MSG_H

# include <qstring.h>
# include <qcolor.h>

namespace config {

  static const QString	MSG_CONFIG_FILE("msg.xml");

  static const QString	MSG_NORMAL_STRING("Message");
  static const QString	MSG_WARNING_STRING("Warning");
  static const QString	MSG_ERROR_STRING("Error");
  static const QString	MSG_DEBUG_STRING("Debug");

  static const QColor	MSG_NORMAL_COLOR(QColor(0, 150, 0));
  static const QColor	MSG_WARNING_COLOR(QColor(200, 100, 0));
  static const QColor	MSG_ERROR_COLOR(QColor(200, 0, 0));
  static const QColor	MSG_DEBUG_COLOR(QColor(180, 180, 0));
  static const QColor	MSG_TAB_NAME_COLOR(QColor(0, 0, 200));

} // End of namespace config.

#endif // CONFIG_MSG_H
