#include "raytracer.h"

char *camera_tos(camera cm)
{
  char buf[256];
  memset(buf,0,256);
  char *s = v3_tos(cm.loc);
  sprintf(buf,"camera(loc=%s,w=%u,h=%u)",s,cm.w,cm.h);
  free(s);
  return strdup(buf);
}

void camera_show(FILE *f, camera cm)
{
  char *s = camera_tos(cm);
  fprintf(f,"%s",s);
  free(s);
}
