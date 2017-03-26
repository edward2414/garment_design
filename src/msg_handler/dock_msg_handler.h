//
//  Filename         : dock_msg_handler.h
//  Author           : Emmanuel Turquin
//  Purpose          : A graphical message handler.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_MSG_HANDLER_H
# define DOCK_MSG_HANDLER_H

# include <qsyntaxhighlighter.h>
# include "dock_window.h"
# include "msg_handler.h"

class QTabWidget;
class Config;

class DockMsgHandler : public DockWindow, public MsgHandler
{
  Q_OBJECT

public:

  DockMsgHandler(MainWindow* parent, const char* name = 0, const QString& prefix = "", Config *conf = 0);
  ~DockMsgHandler();

  void updateTabs();

public slots:

  void message(const MsgHandler::MsgType& type,
	       const QString& from,
	       const QString& msg);
private:

  class SyntaxHighlighter : public QSyntaxHighlighter
  {
  public:

    SyntaxHighlighter(QTextEdit*, const QStringList&, const QValueList<QColor>&);
    virtual ~SyntaxHighlighter();
  
    virtual int highlightParagraph (const QString & text, int endStateOfLastPara);

  private:

    const QStringList&		_keywords;
    const QValueList<QColor>&	_colors;
  };

  QStringList			_keywords;
  QValueList<QColor>		_colors;
  QPtrList<SyntaxHighlighter>	_highlighters;
  QTabWidget			*_tabs;
};

#endif // DOCK_MSG_HANDLER_H
