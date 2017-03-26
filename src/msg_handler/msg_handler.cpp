//
//  Filename         : msg_handler.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : A class to handle messages.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "msg_handler.h"
#include "config_msg.h"

MsgHandler::MsgHandler(const QString& prefix) : _prefix(prefix)
{
  if (!prefix.isEmpty())
    _prefix += "> ";
}

MsgHandler::~MsgHandler()
{
}

void MsgHandler::message(const MsgType& type,
			 const QString& from,
			 const QString& msg) {
  switch (type) {
  case(MSG_DEBUG):
    std::cerr << _prefix << config::MSG_DEBUG_STRING<< " (from " << from << "): " << msg << std::endl;
    break;
  case(MSG_ERROR):
    std::cerr << _prefix << config::MSG_ERROR_STRING<< " (from " << from << "): " << msg << std::endl;
    break;
  case(MSG_WARNING):
    std::cerr << _prefix << config::MSG_WARNING_STRING<< " (from " << from << "): " << msg << std::endl;
    break;
  case(MSG_NORMAL):
  default:
    std::cout << _prefix << config::MSG_NORMAL_STRING << " (from " << from << "): " << msg << std::endl;
  }
}
