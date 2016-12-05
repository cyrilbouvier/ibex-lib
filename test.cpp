#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

struct test_s
{
  int a;
  double b;
  char c;
  int d;
};
typedef struct test_s test_t[1];
typedef struct test_s * test_ptr;

int
main (int argc, char *argv[])
{
  for (int i = 0; i < argc; i++)
    printf ("%s\n", argv[i]);

  double d = 1.0;
  printf ("%a %a %a\n", DBL_MAX, INFINITY, HUGE_VAL);
  printf ("%a %a %a %a %a\n", d, nextafter(d, 0), nextafter (d, DBL_MAX),
                          nextafter (d, INFINITY), nextafter (d, HUGE_VAL));

  test_t t;
  test_ptr ptr = t;

  printf ("%zu %zu %zu\n", sizeof (t), sizeof (ptr), sizeof (struct test_s));
  memset (t, 0xff, sizeof (struct test_s));
  printf ("%x %a %x %x\n", t->a, t->b, t->c, t->d);
  memset (ptr, 0x00, sizeof (struct test_s));
  printf ("%x %a %c %x\n", t->a, t->b, t->c, t->d);

  printf ("%a\n", log (2.));
  d = strtod ("0x1.62e42fefa39efp-1", NULL);
  printf ("%.20f %a => %f %a\n", d, d, pow (exp(1.0), d), pow (exp(1.0), d));
  d = strtod ("0x1.62e42fefa39f0p-1", NULL);
  printf ("%.20f %a => %f %a\n", d, d, pow (exp(1.0), d), pow (exp(1.0), d));

  printf ("%a\n", log10 (3.));
  d = strtod ("0x1.e8927964fd5fcp-2", NULL);
  printf ("%.20f %a => %f %a\n", d, d, pow (10., d), pow (10., d));
  d = strtod ("0x1.e8927964fd5fdp-2", NULL);
  printf ("%.20f %a => %f %a\n", d, d, pow (10., d), pow (10., d));
  d = strtod ("0x1.e8927964fd5fep-2", NULL);
  printf ("%.20f %a => %f %a\n", d, d, pow (10., d), pow (10., d));
  

  return EXIT_SUCCESS;
}
