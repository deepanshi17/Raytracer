#include "raytracer.h"

v3 v3_expr(double x, double y, double z)
{
  v3 result;
  result.x = x;
  result.y = y;
  result.z = z;
  return result;
}

v3 v3_add (v3 v1, v3 v2)
{
  return v3_expr(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

v3 v3_sub (v3 v1, v3 v2)
{

  return v3_expr(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

v3 v3_neg (v3 v)
{
  v3 result;
  result.x = -(v.x);
  result.y = -(v.y);
  result.z = -(v.z);
  return result;
}

v3 v3_scale (double s, v3 v)
{
  v3 result;
  result.x = s*(v.x);
  result.y = s*(v.y);
  result.z = s*(v.z);
  return result;
}

double v3_dot (v3 v1, v3 v2)
{
  return ((v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z));
}

double v3_mag (v3 v)
{
  return sqrt(v3_dot(v, v));
}

v3 v3_norm (v3 v)
{
  return v3_scale((1 / v3_mag(v)), v);
}

char *v3_tos(v3 v)
{
  char buf[128];
  sprintf(buf,"<%.2f,%.2f,%.2f>",v.x,v.y,v.z);
  return strdup(buf);
}

void v3_show(FILE *f, v3 v)
{
  fprintf(f,"<%.2f,%.2f,%.2f>",v.x,v.y,v.z);
}
