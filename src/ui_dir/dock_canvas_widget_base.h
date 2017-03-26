/****************************************************************************
** Form interface generated from reading ui file 'canvas/dock_canvas_widget_base.ui'
**
** Created: Tue Nov 8 16:41:01 2005
**      by: The User Interface Compiler ($Id: dock_canvas_widget_base.h,v 1.1.1.1 2005/11/09 08:01:28 wither Exp $)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#ifndef DOCKCANVASWIDGETBASE_H
#define DOCKCANVASWIDGETBASE_H

#include <qvariant.h>
#include <qpixmap.h>
#include <qwidget.h>
#include "config.h"

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class FloatSpinBox;
class QToolBox;
class QPushButton;
class QGroupBox;
class QLabel;
class QSpinBox;
class QButtonGroup;
class QCheckBox;

class DockCanvasWidgetBase : public QWidget
{
    Q_OBJECT

public:
    DockCanvasWidgetBase( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
    ~DockCanvasWidgetBase();

    QToolBox* toolBox2;
    QWidget* Computation;
    QPushButton* computeButton;
    QGroupBox* groupBox1;
    QLabel* textLabel1_3_2;
    QSpinBox* gmXSpinBox;
    QLabel* textLabel1_3;
    QSpinBox* gmYSpinBox;
    QGroupBox* groupBox1_2;
    QLabel* textLabel1_3_2_2;
    FloatSpinBox* borderXSpinBox;
    QLabel* textLabel1_3_3;
    FloatSpinBox* borderYSpinBox;
    QWidget* DrawingOptions;
    QButtonGroup* buttonGroup2;
    QLabel* textLabel1;
    QSpinBox* samplingSpinBox;
    QLabel* textLabel2;
    QSpinBox* snappingSpinBox;
    QButtonGroup* buttonGroup4_2;
    FloatSpinBox* splittingSpinBox;
    QLabel* textLabel1_2;
    QWidget* Display;
    QButtonGroup* buttonGroup4;
    QCheckBox* strokeCheckBox;
    QCheckBox* chainCheckBox;
    QCheckBox* antialiasingCheckBox;
    QCheckBox* typeCheckBox;
    QCheckBox* textureCheckBox;
    QCheckBox* bbCheckBox;

public slots:
    virtual void splittingSpinBox_valueChanged(float);
    virtual void samplingSpinBox_valueChanged(int);
    virtual void snappingSpinBox_valueChanged(int);
    virtual void chainCheckBox_toggled(bool);
    virtual void strokeCheckBox_toggled(bool);
    virtual void typeCheckBox_toggled(bool);
    virtual void antialiasingCheckBox_toggled(bool);
    virtual void textureCheckBox_toggled(bool);
    virtual void bbCheckBox_toggled(bool);
    virtual void gmXSpinBox_valueChanged(int);
    virtual void gmYSpinBox_valueChanged(int);
    virtual void borderXSpinBox_valueChanged(float);
    virtual void borderYSpinBox_valueChanged(float);
    virtual void computeButton_clicked();

protected:
    QGridLayout* DockCanvasWidgetBaseLayout;
    QGridLayout* ComputationLayout;
    QSpacerItem* spacer5_2;
    QGridLayout* groupBox1Layout;
    QGridLayout* groupBox1_2Layout;
    QGridLayout* DrawingOptionsLayout;
    QSpacerItem* spacer5;
    QGridLayout* buttonGroup2Layout;
    QGridLayout* buttonGroup4_2Layout;
    QGridLayout* DisplayLayout;
    QSpacerItem* spacer4;
    QGridLayout* buttonGroup4Layout;

protected slots:
    virtual void languageChange();

private:
    QPixmap image0;

};

#endif // DOCKCANVASWIDGETBASE_H
