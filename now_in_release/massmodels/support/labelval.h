// Handy for reading label = value pairs

#ifndef LABELVAL_H
#define LABELVAL_H 1

#include <iostream>
#include <string>


// For parameter values
class Labelval {
private:
  std::string label;
  double value;
    
public:
  friend std::istream& operator>> (std::istream& sce, Labelval& res);
  friend std::ostream& operator<< (std::ostream& dest, Labelval& x);
  bool match (const char* s) const {return label.find (s) == 0;}
  double val () const {return value;}
  double operator() () const {return value;}
};

#endif
