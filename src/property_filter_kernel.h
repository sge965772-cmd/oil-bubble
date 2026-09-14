#ifndef PROPERTY_FILTER_KERNEL_H
#define PROPERTY_FILTER_KERNEL_H

static inline double property_filter_2d (
  double center,
  double north, double south, double east, double west,
  double northeast, double northwest, double southeast, double southwest)
{
  return (4.*center + 2.*(north + south + east + west) +
          northeast + northwest + southeast + southwest)/16.;
}

#endif
