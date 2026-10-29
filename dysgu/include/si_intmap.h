#ifndef DYSGU_SI_INTMAP_H
#define DYSGU_SI_INTMAP_H
// Convenience typedef for the int/int specialization used by the fragment
// genotyper (dysgu/extra_metrics.pyx). Cython cannot name a C++ template
// specialization directly, so bind this typedef instead.
#include "superintervals.h"
typedef si::IntervalMap<int, int> SiIntMap;
#endif
