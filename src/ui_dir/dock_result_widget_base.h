/****************************************************************************
** Form interface generated from reading ui file 'result/dock_result_widget_base.ui'
**
** Created: Tue Nov 8 16:41:01 2005
**      by: The User Interface Compiler ($Id: dock_result_widget_base.h,v 1.1.1.1 2005/11/09 08:01:28 wither Exp $)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#ifndef DOCKRESULTWIDGETBASE_H
#define DOCKRESULTWIDGETBASE_H

#include <qvariant.h>
#include <qwidget.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class QToolBox;
class QCheckBox;
class QGroupBox;

class DockResultWidgetBase : public QWidget
{
    Q_OBJECT

public:
    DockResultWidgetBase( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
    ~DockResultWidgetBase();

    QToolBox* toolBox2;
    QWidget* Display;
    QCheckBox* modelCheckBox;
    QGroupBox* groupBox1;
    QCheckBox* contourCheckBox;
    QCheckBox* typeCheckBox;
    QCheckBox* surfaceCheckBox;

public slots:
    virtual void modelCheckBox_toggled(bool);
    virtual void surfaceCheckBox_toggled(bool);
    virtual void contourCheckBox_toggled(bool);
    virtual void typeCheckBox_toggled(bool);

protected:
    QGridLayout* DockResultWidgetBaseLayout;
    QGridLayout* DisplayLayout;
    QSpacerItem* spacer4;
    QGridLayout* groupBox1Layout;

protected slots:
    virtual void languageChange();

};

#endif // DOCKRESULTWIDGETBASE_H
