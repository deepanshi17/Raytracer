#include "raytracer.h"

v3 ray_position (ray r, double t)
{
  return v3_add(r.origin, v3_scale(t, r.direction));
}

char *ray_tos(ray r)
{
  char *o = v3_tos(r.origin);
  char *d = v3_tos(r.direction);
  char buf[128];
  sprintf(buf,"ray(%s,%s)",o,d);
  free(o);
  free(d);
  return strdup(buf);
}

void ray_show(FILE *f, ray r)
{
  char *s = ray_tos(r);
  fprintf(f,"%s",s);
  free(s);
}

