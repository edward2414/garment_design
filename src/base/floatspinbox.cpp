//
//  Filename         : floatspinbox.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : A spin box that works with floats.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qspinbox.h>
#include <limits>
#include <cmath>
#include "floatspinbox.h"

FloatSpinBox::FloatSpinBox( QWidget* parent, const char* name )
  : QSpinBox(0, std::numeric_limits<int>::max(), 1000, parent, name), dec(3)
{
  dval = new QDoubleValidator(std::numeric_limits<float>::min(), std::numeric_limits<float>::max(), 2, this, "dval" );
  setValidator( dval );
  setValue(1.0); // Here to force displaying 0.000 just after...
  setValue(0.0);
}

FloatSpinBox::~FloatSpinBox()
{
  delete dval;
}  

QString FloatSpinBox::mapValueToText ( int value )
{
  QString s;
  s.setNum( float( value )/ pow( 10, dec ), 'f', dec );
  return s;
}

int FloatSpinBox::mapTextToValue ( bool * ok )
{
  return int( cleanText().toFloat( ok ) * pow( 10, dec ) );
}

float FloatSpinBox::value() const
{
  return float( QRangeControl::value() ) / pow( 10, dec );
}

float FloatSpinBox::minValue() const
{
  return float( QRangeControl::minValue() ) / pow( 10, dec );
}

float FloatSpinBox::maxValue() const
{
  return float( QRangeControl::maxValue() ) / pow( 10, dec );
}

float FloatSpinBox::lineStep() const
{
  return float(QSpinBox::lineStep())/ pow(10, dec);
}

void FloatSpinBox::setLineStep (float value)
{
  QSpinBox::setLineStep(int(value *  pow(10, dec)));
}

void FloatSpinBox::setValue( float value )
{
  QRangeControl::setValue( int( value *  pow( 10, dec ) ) );
}

void FloatSpinBox::setRange( float minValue, float maxValue )
{
  QRangeControl::setRange( int( minValue *  pow( 10, dec ) ), 
			   int( maxValue *  pow( 10, dec ) ) );
  dval->setRange( minValue, maxValue, 2 );
}

void FloatSpinBox::valueChange()
{
  QSpinBox::valueChange();
  emit valueChanged( value() );
}
