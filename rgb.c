#include "raytracer.h"

void check_color_component(char c, double x)
{
  if (x<0||x>1) {
    fprintf(stderr,"rgb_expr: %c out of range: %f\nhalting\n",c,x);
    exit(1);
  }
}

double clamp1(double x)
{
  return x>1?1:x;
}

rgb rgb_expr(double r, double g, double b)
{
  check_color_component('r',r);
  check_color_component('g',g);
  check_color_component('b',b);
  rgb color = {r,g,b};
  return color;
}

rgb rgb_modulate (rgb c1, rgb c2)
{
  return rgb_expr((c1.r * c2.r), (c1.g * c2.g), (c1.b * c2.b));
}

rgb rgb_scale (double s, rgb c)
{
  return rgb_expr(clamp1(s * c.r), clamp1(s * c.g), clamp1(s * c.b));
}

rgb rgb_add (rgb c1, rgb c2)
{
  return rgb_expr(clamp1((c1.r + c2.r)), clamp1((c1.g + c2.g)), clamp1((c1.b + c2.b)));
}

char *rgb_tos(rgb c)
{
  char buf[128];
  sprintf(buf,"rgb(%.2f,%.2f,%.2f)",c.r,c.g,c.b);
  return strdup(buf);
}

void rgb_show(FILE *f, rgb c)
{
  fprintf(f,"rgb(%.2f,%.2f,%.2f)",c.r,c.g,c.b);
}

byte discretize (double d)
{
  return (byte)(255.0*d);
}

void rgb_show_bytes (FILE* f, rgb c)
{
  fprintf(f, "%d %d %d\n", discretize(c.r), discretize(c.g), discretize(c.b));
}


  
