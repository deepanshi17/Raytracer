#include "raytracer.h"

poster poster_expr (v3 ult, double w, double h, rgb (*sc) (v3, v3), rgb sh)
{
  poster result;
  result.upper_left = ult;
  result.w = w;
  result.h = h;
  result.surface_color = sc;
  result.shine = sh;
  return result;
}

char *poster_tos(poster p)
{
  char buf[256];
  memset(buf,'\0',256);
  char *s1 = v3_tos(p.upper_left);
  char *s2 = rgb_tos(p.shine);
  snprintf(buf,256,"poster(%s,%.2fx%.2f,<FN>,%s)",
	   s1,p.w,p.h,s2);
  free(s1);
  free(s2);
  return strdup(buf);
}

void poster_show(FILE *f, poster p)
{
  char *z = poster_tos(p);
  fprintf(f,"%s",z);
  free(z);
}

