//
//  Filename         : base_qglviewer.h
//  Author           : Emmanuel Turquin
//  Purpose          : An intermediate class between QGLViewer and
//                     the application's viewers.
//  Date of creation : 04/21/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  BASE_QGLVIEWER_H
# define BASE_QGLVIEWER_H

# include <QGLViewer/qglviewer.h>
# include "msg_handler.h"

class BaseQGLViewer : public QGLViewer
{
  Q_OBJECT

public:

  BaseQGLViewer(QWidget* parent=NULL, const char* name=0,
		const QGLWidget* shareWidget=0, int wflags=0);

  void setSceneBoundingBox(const double min[3], const double max[3]) {
    //const float fmin[3] = {min[0], min[1], min[2]};
    //const float fmax[3] = {max[0], max[1], max[2]};

    // FIXME JDW method below should accept the float arrays above and 
    // convert them to qglviewer::Vec types according to the 
    // QGLViewer documentation for the Vec class. 
    //Unfortunately On this development machine
    //it causes a compiler error and I don't have time to 
    //work out why. So for now convert explicitly
 
    const qglviewer::Vec fmin( min[0], min[1], min[2] );
    const qglviewer::Vec fmax( max[0], max[1], max[2] );
    QGLViewer::setSceneBoundingBox(fmin, fmax);
  }

signals:

  void statusBarMessage(const QString&, int);
  void message(const MsgHandler::MsgType&, const QString&, const QString&);
};

#endif // BASE_QGLVIEWER_H
