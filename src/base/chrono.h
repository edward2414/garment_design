//
//  Filename         : chrono.h
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Class to measure ellapsed time in seconds.
//  Date of creation : 05/18/04
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  CHRONO_H
# define CHRONO_H

#include <time.h>

class Chronometer
{
 public:

  inline Chronometer() {}
  inline ~Chronometer() {}

  inline clock_t start() {
    _start = clock();
    return _start;
  }

  inline double stop() {
    clock_t stop = clock();
    return (double)(stop - _start) / CLOCKS_PER_SEC ;
  }

 private:

  clock_t _start;
};

#endif // CHRONO_H
