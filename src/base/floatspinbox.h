//
//  Filename         : floatspinbox.h
//  Author           : Emmanuel Turquin
//  Purpose          : A spin box that works with floats.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  FLOATSPINBOX_H
# define FLOATSPINBOX_H

# include <qspinbox.h>
# include <qvalidator.h>

class FloatSpinBox : public QSpinBox
{
  Q_OBJECT

 public:

  FloatSpinBox(QWidget* parent = 0, const char* name = 0);
  ~FloatSpinBox();

  float value() const;
  float minValue() const;
  float maxValue() const;

  float lineStep () const;
  void setLineStep (float);

  void setRange(float minValue, float maxValue);

 public slots:

  void setValue(float value);

 signals:

  void valueChanged(float value);

 protected:

  virtual QString mapValueToText (int value);
  virtual int mapTextToValue (bool * ok = 0);
  virtual void valueChange ();

 private:

  QDoubleValidator* dval;

  int dec;
};

#endif // FLOATSPINBOX_H
