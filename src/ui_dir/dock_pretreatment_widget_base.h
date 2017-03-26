/****************************************************************************
** Form interface generated from reading ui file 'pretreatment/dock_pretreatment_widget_base.ui'
**
** Created: Tue Nov 8 16:41:01 2005
**      by: The User Interface Compiler ($Id: dock_pretreatment_widget_base.h,v 1.1.1.1 2005/11/09 08:01:28 wither Exp $)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#ifndef DOCKPRETREATMENTWIDGETBASE_H
#define DOCKPRETREATMENTWIDGETBASE_H

#include <qvariant.h>
#include <qpixmap.h>
#include <qwidget.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class FloatSpinBox;
class QToolBox;
class QGroupBox;
class QLabel;
class QPushButton;
class QSpinBox;
class QCheckBox;
class QButtonGroup;
class QRadioButton;

class DockPretreatmentWidgetBase : public QWidget
{
    Q_OBJECT

public:
    DockPretreatmentWidgetBase( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
    ~DockPretreatmentWidgetBase();

    QToolBox* toolBox;
    QWidget* Computation;
    QGroupBox* groupBox2;
    FloatSpinBox* bbXSpinBox;
    FloatSpinBox* bbYSpinBox;
    FloatSpinBox* bbZSpinBox;
    QLabel* textLabel1_4;
    QLabel* textLabel1_2_2;
    QLabel* textLabel1_3_2;
    QPushButton* bbResetButton;
    QGroupBox* groupBox3_2;
    QLabel* textLabel1_7;
    QSpinBox* textureSizeSpinBox;
    QGroupBox* groupBox3;
    QSpinBox* octreeLevelSpinBox;
    QLabel* textLabel1_5;
    QGroupBox* groupBox4;
    QLabel* textLabel2_2;
    QSpinBox* dfXSpinBox;
    QSpinBox* dfYSpinBox;
    QSpinBox* dfZSpinBox;
    QLabel* textLabel5;
    QLabel* textLabel4;
    QLabel* textLabel3;
    QPushButton* computeButton;
    QWidget* Display;
    QCheckBox* modelCheckBox;
    QCheckBox* bbCheckBox;
    QCheckBox* bbAxisCheckBox;
    QCheckBox* octreeCheckBox;
    QButtonGroup* buttonGroup1;
    QCheckBox* dfCheckBox;
    QLabel* textLabel1;
    QLabel* textLabel1_2;
    QLabel* textLabel1_3;
    QLabel* textLabel1_6;
    QSpinBox* dfXMinSpinBox;
    QSpinBox* dfXMaxSpinBox;
    QSpinBox* dfYMinSpinBox;
    QSpinBox* dfYMaxSpinBox;
    QSpinBox* dfZMinSpinBox;
    QSpinBox* dfZMaxSpinBox;
    QLabel* textLabel2;
    QButtonGroup* buttonGroup2;
    QRadioButton* perspectiveButton;
    QRadioButton* orthoButton;

public slots:
    virtual void perspectiveButton_toggled(bool);
    virtual void orthoButton_toggled(bool);
    virtual void bbCheckBox_toggled(bool);
    virtual void modelCheckBox_toggled(bool);
    virtual void bbZSpinBox_valueChanged(float);
    virtual void bbYSpinBox_valueChanged(float);
    virtual void bbXSpinBox_valueChanged(float);
    virtual void octreeCheckBox_toggled(bool);
    virtual void octreeLevelSpinBox_valueChanged(int);
    virtual void bbAxisCheckBox_toggled(bool);
    virtual void computeButton_clicked();
    virtual void bbResetButton_clicked();
    virtual void textureSizeSpinBox_valueChanged(int);
    virtual void dfXSpinBox_valueChanged(int);
    virtual void dfYSpinBox_valueChanged(int);
    virtual void dfZSpinBox_valueChanged(int);
    virtual void dfCheckBox_toggled(bool);
    virtual void dfXMinSpinBox_valueChanged(int);
    virtual void dfXMaxSpinBox_valueChanged(int);
    virtual void dfYMinSpinBox_valueChanged(int);
    virtual void dfYMaxSpinBox_valueChanged(int);
    virtual void dfZMinSpinBox_valueChanged(int);
    virtual void dfZMaxSpinBox_valueChanged(int);

protected:
    QGridLayout* DockPretreatmentWidgetBaseLayout;
    QGridLayout* ComputationLayout;
    QSpacerItem* spacer6;
    QGridLayout* groupBox2Layout;
    QGridLayout* groupBox3_2Layout;
    QGridLayout* groupBox3Layout;
    QGridLayout* groupBox4Layout;
    QGridLayout* DisplayLayout;
    QSpacerItem* spacer7;
    QGridLayout* buttonGroup1Layout;
    QGridLayout* buttonGroup2Layout;

protected slots:
    virtual void languageChange();

private:
    QPixmap image0;

};

#endif // DOCKPRETREATMENTWIDGETBASE_H
