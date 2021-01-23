#include "raytracer.h"

char *light_tos(light l)
{
  char buf[256];
  memset(buf,0,256);
  char *s1 = v3_tos(l.direction);
  char *s2 = rgb_tos(l.color);
  snprintf(buf,256,"light(dir=%s,color=%s)",s1,s2);
  free(s1);
  free(s2);
  return strdup(buf);
}

void light_show(FILE *f, light l)
{
  char *s = light_tos(l);
  fprintf(f,"%s",s);
  free(s);
}

