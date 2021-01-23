#include "raytracer.h"

char* hit_test_tos (hit_test t)
{
  char* string = malloc(sizeof(char*));
  if (0 > sprintf(string, "%d %lf %s %s %s %s",
		  t.miss,
		  t.t,
		  v3_tos(t.hit_point),
		  rgb_tos(t.surf),
		  rgb_tos(t.shine),
		  v3_tos(t.surf_norm)))
  {
    fprintf(stderr, "invalid inputs to hit_test_tos\n");
    free(string);
    return NULL;
  }
  else return string;
}

void hit_test_show (FILE* f, hit_test t)
{
  char* string = hit_test_tos(t);
  fprintf(f, "Hit Test: %s\n", hit_test_tos(t));
  free(string);
  return;
}


