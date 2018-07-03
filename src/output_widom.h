#ifndef OUTPUT_H
#define OUTPUT_H

#include "stdlib.h"
#include <iostream>

class OutputWidom {
  public:
    OutputWidom();

    void header( std::ostream& );
    void add_line( std::ostream& , int, double, double, double );

};
#endif
