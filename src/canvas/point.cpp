#include "point.h"
#include <fstream>

ostream& operator<<(ostream& out, const Point& point) {
  out << "(" << point._coord[0] << ", " << point._coord[1] << ")";
  //out << point._coord[0];
  //out << " ";
  //out << point._coord[1];
  return(out); 
}

istream& operator>>(istream& in, Point& point) {
  char c; // use to strip the extra formatting to make it more human readable
  in >> c;
  in >> point._coord[0];
  in >> c;
  in >> point._coord[1];
  in >> c;
  return(in); 
}

/*
int main(int argc, char **argv )
{
  Point p(0.5,0.6);
  cout << p << endl;
  
  ofstream out("/tmp/point.txt");
  if(!out)
  {
    cerr << "Unable to open output file.\n";
    exit(EXIT_FAILURE); 
  }
  out << p;
  out.close();
  
  ifstream in("/tmp/point.txt");
  if(!in)
  {
    cerr << "Unable to open input file.\n";
    exit(EXIT_FAILURE); 
  }
  
  Point q(0,0);
  in >> q;
  in.close();
  cout << q << endl;
}
*/