//
//  Filename         : dock_msg_handler.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : A graphical message handler.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qworkspace.h>
#include <qtabwidget.h>
#include <qtextedit.h>
#include "main_window.h"
#include "mdi_window.h"
#include "dock_msg_handler.h"
#include "config.h"
#include "config_msg.h"

DockMsgHandler::DockMsgHandler(MainWindow* parent, const char* name, const QString& prefix, Config *conf)
  : DockWindow(parent, name, 0, conf), MsgHandler(prefix)
{
  // Initialize keywords for syntax highlighting.
  _keywords.push_back(config::MSG_NORMAL_STRING);
  _colors.push_back(config::MSG_NORMAL_COLOR);
  _keywords.push_back(config::MSG_WARNING_STRING);
  _colors.push_back(config::MSG_WARNING_COLOR);
  _keywords.push_back(config::MSG_ERROR_STRING);
  _colors.push_back(config::MSG_ERROR_COLOR);
  _keywords.push_back(config::MSG_DEBUG_STRING);
  _colors.push_back(config::MSG_DEBUG_COLOR);

  _tabs = new QTabWidget(this);
  setWidget(_tabs);
  QTextEdit *te =  new QTextEdit();
  SyntaxHighlighter *sh = new SyntaxHighlighter(te, _keywords, _colors);
  _highlighters.append(sh);
  te->setReadOnly(true);
  _tabs->addTab(te, "All");

  if (parent) {
    parent->setDockEnabled(this, Qt::DockLeft, false);
    parent->setDockEnabled(this, Qt::DockRight, false);
    if (_place == UNDEFINED) {
      parent->moveDockWindow(this, Qt::DockBottom);
      if (_conf)
	_conf->io()->setValue("dock_window/place", DOCK_BOTTOM);
    }
    connect(parent, SIGNAL(message(const MsgHandler::MsgType&, const QString&, const QString&)),
	    this, SLOT(message(const MsgHandler::MsgType&, const QString&, const QString&)));
  }
  _tabs->show();
  updateTabs();
}

DockMsgHandler::~DockMsgHandler()
{
  for (QPtrList<SyntaxHighlighter>::iterator it = _highlighters.begin(), it_end = _highlighters.end();
       it != it_end;
       ++it) {
    delete(*it);
  }
  for (int i = 0; i < _tabs->count(); ++i)
    delete _tabs->page(i);
  delete _tabs;
}

void DockMsgHandler::updateTabs()
{
  MainWindow *mw = dynamic_cast<MainWindow*>(parent());
  if (!mw)
    mw = dynamic_cast<MainWindow*>(topLevelWidget());
  if (!mw)
    return;
  QWorkspace *ws = mw->workspace();
  if (!ws)
    return;
  QWidgetList mdi_windows = ws->windowList(QWorkspace::CreationOrder);
  for (int i = _tabs->count() - 1; i < int(mdi_windows.count()); ++i) {
    MDIWindow *m = (MDIWindow*)mdi_windows.at(i);
    if (m) {
      // update keywords for syntax highlighting.
      _keywords.push_back(m->caption());
      _colors.push_back(config::MSG_TAB_NAME_COLOR);

      QTextEdit *te =  new QTextEdit();
      SyntaxHighlighter *sh = new SyntaxHighlighter(te, _keywords, _colors);
      _highlighters.append(sh);
      te->setReadOnly(true);
      _tabs->addTab(te, m->caption());
    }
  }
}

void DockMsgHandler::message(const MsgHandler::MsgType& type,
			     const QString& from,
			     const QString& msg) {
  int found = -1;
  for (int i = 0; i < _tabs->count(); ++i) {
    if (from == _tabs->label(i)) {
      found = i;
      break;
    }
  }
  switch (type) {
  case(MSG_DEBUG):
    ((QTextEdit*)(_tabs->page(0)))->append(_prefix + config::MSG_DEBUG_STRING +
					   " (from " + from + "): " + msg + "\n");
    if (found >= 0)
      ((QTextEdit*)(_tabs->page(found)))->append(_prefix + config::MSG_DEBUG_STRING +": " + msg + "\n");
    break;
  case(MSG_ERROR):
    ((QTextEdit*)(_tabs->page(0)))->append(_prefix + config::MSG_ERROR_STRING +
					   " (from " + from + "): " + msg + "\n");
    if (found >= 0)
      ((QTextEdit*)(_tabs->page(found)))->append(_prefix + config::MSG_ERROR_STRING + ": " + msg + "\n");
    break;
  case(MSG_WARNING):
    ((QTextEdit*)(_tabs->page(0)))->append(_prefix + config::MSG_WARNING_STRING +
					   " (from " + from + "): " + msg + "\n");
    if (found >= 0)
      ((QTextEdit*)(_tabs->page(found)))->append(_prefix + config::MSG_WARNING_STRING + ": " + msg + "\n");
    break;
  case(MSG_NORMAL):
  default:
    ((QTextEdit*)(_tabs->page(0)))->append(_prefix + config::MSG_NORMAL_STRING +
					   " (from " + from + "): " + msg + "\n");
    if (found >= 0)
      ((QTextEdit*)(_tabs->page(found)))->append(_prefix + config::MSG_NORMAL_STRING + ": " + msg + "\n");
  }
}

DockMsgHandler::SyntaxHighlighter::SyntaxHighlighter(QTextEdit *te,
						     const QStringList& keywords,
						     const QValueList<QColor>& colors) :
  QSyntaxHighlighter(te), _keywords(keywords), _colors(colors)
{
}

DockMsgHandler::SyntaxHighlighter::~SyntaxHighlighter()
{
}

int DockMsgHandler::SyntaxHighlighter::highlightParagraph(const QString & text, int endStateOfLastPara) {
  int pos = 0;
  for (unsigned i = 0; i < _keywords.count(); ++i) {
    pos = 0;
    QString word = _keywords[i];
    while ((pos = text.find(word, pos)) != -1) {
      setFormat(pos, word.length(), _colors[i]);
      pos += text.length() + 1;
    }
  }
  return endStateOfLastPara;
}
