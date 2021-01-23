#include "raytracer.h"

char *cylinder_tos(cylinder c)
{
  char buf[256];
  memset(buf,'\0',256);
  char* s1 = v3_tos(c.center);
  char* s2 = rgb_tos(c.shine);
  snprintf(buf,256,"cylinder(%s,%.2f,%.2f,<FN>,%s)",
	   s1,c.radius,c.height,s2);
  free(s1);
  free(s2);
  return strdup(buf);
}

char *object_tos(object obj)
{
  switch (obj.tag) {
  case SPHERE  : return sphere_tos(obj.o.s);
  case POSTER : return poster_tos(obj.o.p);
  case CYLINDER : return cylinder_tos(obj.o.cyl);
  }
  fprintf(stderr,"object_tos: bad tag %d",obj.tag);
  exit(1);
}

void object_show(FILE *f, object obj)
{
  char *s = object_tos(obj);
  fprintf(f,"%s",s);
  free(s);
}
