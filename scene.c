#include "raytracer.h"

#define SCENE_FMT_STRING "--- begin scene\nbg:    %s\nlight: %s\namb:   %s\nobjs:  %s\n--- end scene\n"

char *scene_tos(scene sc)
{
  char buf[1024];
  memset(buf,0,1024);
  char *bgs  = rgb_tos(sc.bg);
  char *ls   = light_tos(sc.light);
  char *ambs = rgb_tos(sc.amb);
  char *os   = ol_tos(sc.objects);
  snprintf(buf,1024,SCENE_FMT_STRING,bgs,ls,ambs,os);
  char *result = strdup(buf);
  free(bgs);
  free(ls);
  free(ambs);
  free(os);
  return result;
}

void scene_show(FILE *f, scene sc)
{
  char *s = scene_tos(sc);
  fprintf(f,"%s",s);
  free(s);
}
