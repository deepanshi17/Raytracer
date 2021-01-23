#include <stdlib.h>
#include <stdio.h>
#include "raytracer.h"

double double_max (double a, double b)
{
  return (a > b)? a : b;
}

int shadow (v3 p, light l, obj_list* objs)
{
  if (objs == NULL)
    return 0;
  else {
  ray shadow_ray;
  shadow_ray.direction = l.direction;
  v3 nudge = v3_scale(0.0001, l.direction);
  v3 lifted = v3_add(p, nudge);
  shadow_ray.origin = lifted;
  hit_test result = intersect(shadow_ray, objs->obj);
  if (result.miss == 0)
    return 1;
  else return (shadow(p, l, objs->next));
  }
}

rgb lighting (scene s, ray r, hit_test h)
{
  rgb result;
  if (h.miss)
    return s.bg;
  if (shadow(h.hit_point, s.light, s.objects)) {
      result = rgb_modulate(h.surf, s.amb);
    }
  else {
    result = rgb_modulate(h.surf, rgb_add(s.amb,
		rgb_scale(double_max(0, v3_dot(h.surf_norm, s.light.direction)),
		          s.light.color)));
  }
  v3 n = h.surf_norm;
  v3 l = s.light.direction;
  rgb so = h.shine;
  v3 ro = v3_sub(v3_scale((2.0 * v3_dot(n, l)), n), l);
  v3 v = v3_neg(r.direction);
  rgb specular = {0,0,0};
  if (v3_dot(n, l) > 0)
    specular = rgb_scale(pow(double_max(v3_dot(ro, v), 0), 6.0), so);
  return rgb_add(result, specular);
}

v3 logical_loc (camera c, unsigned int x, unsigned int y)
{
  double x_offset = -1.0, y_offset = 1.0, pix;
  if (c.w < c.h) {
    x_offset = -(double)c.w/c.h;
    pix = (double)2/c.h;
  } else if (c.w > c.h) {
    y_offset = (double)c.h/c.w;
    pix = (double)2/c.w;
  } else
    pix = (double)2/c.w;
  return v3_expr(x_offset + (x+0.5)*pix, y_offset - (y+0.5)*pix, 0);
}
  
rgb trace_ray (scene s, ray r)
{
  int hit = 0;
  hit_test best, tmp;
  best.t = INFINITY;

  while(s.objects) {
    tmp = intersect(r, s.objects->obj);
    if (!tmp.miss && tmp.t < best.t) {
      best = tmp;
      hit = 1;
    } s.objects = s.objects->next;
  }
  if (hit)
    return lighting(s,r,best);
  return lighting(s,r,tmp);
}

void render (FILE* f, stage g)
{
  int i, j;
  printf("P3\n");
  printf("%d %d\n", g.c.w, g.c.h);
  printf("255\n");
  v3 dir;

  for (i=0; i < g.c.h; i++) {
    for (j=0; j < g.c.w; j++) {
      dir = v3_norm(v3_sub(logical_loc(g.c, j, i), g.c.loc));
      ray r = {g.c.loc, dir};
      rgb_show_bytes(f, trace_ray(g.s, r));
    }
  }
}
