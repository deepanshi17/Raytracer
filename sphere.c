#include "raytracer.h"

sphere sphere_expr (v3 v, double r, rgb (*sc) (v3, v3), rgb sh)
{
  sphere result;
  result.center = v;
  result.radius = r;
  result.surface_color = sc;
  result.shine = sh;
  return result;
}

char *sphere_tos(sphere s)
{
  char buf[256];
  memset(buf,0,256);
  char *s1 = v3_tos(s.center);
  char *s2 = rgb_tos(s.shine);
  snprintf(buf,256,"sph(%s,%2.2f,<FN>,%s",s1,s.radius,s2);
  free(s1);
  free(s2);
  return strdup(buf);
}

void sphere_show(FILE *f, sphere s)
{
  char *z = sphere_tos(s);
  fprintf(f,"%s",z);
  free(z);
}

