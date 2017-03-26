//
//  Filename         : msg_handler.h
//  Author           : Emmanuel Turquin
//  Purpose          : A class to handle messages.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  MSG_HANDLER_H
# define MSG_HANDLER_H

# include <iostream>
# include <qstring.h>

class MsgHandler
{
public:

  typedef enum {
    MSG_NORMAL,
    MSG_WARNING,
    MSG_ERROR,
    MSG_DEBUG
  } MsgType;

  MsgHandler(const QString& prefix = "");
  virtual ~MsgHandler();

  virtual void message(const MsgType& type,
		       const QString& from,
		       const QString& msg);

protected:

  QString _prefix;
};

#endif // MSG_HANDLER_H
