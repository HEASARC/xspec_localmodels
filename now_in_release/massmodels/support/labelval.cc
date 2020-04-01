// Handling of label = value pairs for reading parameters

#include "labelval.h"

// For reading values from the file of constants
std::istream& operator>> (std::istream& sce, Labelval& res) {
  // No format checks - label includes trailing space
  getline (sce, res.label, '=');
  sce >> res.value;
  // Ready for next entry
  sce.ignore (256, '\n');
  return sce;
}

std::ostream& operator<< (std::ostream& dest, Labelval& x) {
  dest << x.label << "= " << x.value;
  return dest;
}
