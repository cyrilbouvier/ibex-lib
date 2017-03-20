#include <cmath>
#include <limits>
#include <cassert>
#include <stdio.h>

int main(void)
{
  double a = -std::numeric_limits<double>::infinity();
  double b = nextafter(-std::numeric_limits<double>::max(),
                                    -std::numeric_limits<double>::infinity());
  printf ("%a\n%a\n%d\n", a, b, (a == b));
  assert(a == b);
}

