#include "raytracer.h"

int main (int argc, char* argv[])
{
  stage stg = get_stage(3);
  render(stdout, stg);
  stage_show(stderr, stg);
  return 0;
}
