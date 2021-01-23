#include "raytracer.h"

#define STAGE_FMT "*** begin stage\n\n%s\n\n%s\n*** end stage\n"

char *stage_tos(stage stg)
{
  char buf[1024];
  memset(buf,0,1024);
  char *c = camera_tos(stg.c);
  char *s = scene_tos(stg.s);
  snprintf(buf,1024,STAGE_FMT,c,s);
  free(c);
  free(s);
  return strdup(buf);
}

void stage_show(FILE *f, stage stg)
{
  char *s = stage_tos(stg);
  fprintf(f,"%s",s);
  free(s);
}


rgb const_color1 (v3 a, v3 b)
{
  return rgb_expr(0.5, 0.5, 0.8);
}

rgb const_color2 (v3 a, v3 b)
{
  return rgb_expr(1.0, 0, 0);
}

rgb k_silver (v3 cen, v3 loc)
{
  return rgb_expr(0.75, 0.75, 0.75);
}

cylinder cylinder_expr (v3 cn, double r, double h, rgb (*sc) (v3, v3), rgb sh)
{
  cylinder result;
  result.center = cn;
  result.radius = r;
  result.height = h;
  result.surface_color = sc;
  result.shine = sh;
  return result;
}


stage get_stage (unsigned int selection)
{
  stage stg;
  light l;
  poster p = poster_expr(v3_expr(-1, 1, 4), 2, 1, &const_color1, rgb_expr(0, 0, 0));
  sphere s = sphere_expr(v3_expr(0, 0, 3), 1.0, &const_color1, rgb_expr(0, 0, 0));
  object o, o2, o3;
  obj_list* objs;

  switch (selection) {
  case 1 :
    l.direction = v3_norm(v3_expr(-1, 1, -1));
    l.color = rgb_expr(1, 1, 1);
    o.tag = SPHERE;
    o.o.s = s;
    objs = ol_singleton(o);
    camera test_camera_1 = {v3_expr(0, 0, -5.0), 240, 240};
    scene test_scene_1 = {rgb_expr(0.3, 0.6, 1),
			  objs,
			  l,
			  rgb_expr(0.2, 0.2, 0.2)};
    stg.c = test_camera_1;
    stg.s = test_scene_1;
    break;

  case 2 :
    l.direction = v3_norm(v3_expr(1, -1, 1));
    l.color = rgb_expr(1, 1, 1);
    o.tag = POSTER;
    o.o.p = p;
    o2.tag = SPHERE;
    o2.o.s = s;
    objs = ol_cons(o, ol_cons(o2, NULL));
    camera test_camera_2 = {v3_expr(0, 0, -5.0), 240, 240};
    scene test_scene_2 = {rgb_expr(0.3, 0.6, 1.0),
			  objs,
			  l,
			  rgb_expr(0.2, 0.2, 0.2)};
    stg.c = test_camera_2;
    stg.s = test_scene_2;
    break;

  case 3 :
    l.direction = v3_expr(-0.58, 0.58, -0.58);
    l.color = rgb_expr(1, 1, 1);
    cylinder c1 = cylinder_expr(v3_expr(0, -3, 0), 0.05, 6, &const_color2, rgb_expr(0, 0, 0));
    cylinder c2 = cylinder_expr(v3_expr(0, 0.3, 0), 0.8, 0.05, &k_silver, rgb_expr(0,0,0));
    sphere s2 = sphere_expr(v3_expr(0, 0.3, 0), 0.4, &const_color1, rgb_expr(0, 0, 0));
    o.tag = CYLINDER;
    o.o.cyl = c1;
    o2.tag = CYLINDER;
    o2.o.cyl = c2;
    o3.tag = SPHERE;
    o3.o.s = s2;
    objs = ol_cons(o, ol_cons(o2, ol_cons(o3, NULL)));
    camera test_camera_3 = {v3_expr(0, 0, -5), 600, 700};
    scene test_scene_3 = {rgb_expr(0, 0, 0.25),
			  objs,
			  l,
			  rgb_expr(0.1, 0.1, 0.1)};
    stg.c = test_camera_3;
    stg.s = test_scene_3;
    break;

  case 4 :
    l.direction = v3_expr(-0.58, 0.58, -0.58);
    l.color = rgb_expr(1, 1, 1);
    cylinder c3 = cylinder_expr(v3_expr(0, -2, 0), 0.2, 4, &k_silver, rgb_expr(0, 0, 0));
    o.tag = CYLINDER;
    o.o.cyl = c3;
    objs = ol_singleton(o);
    camera test_camera_4 = {v3_expr(0, 0, -5), 120, 120};
    scene test_scene_4 = {rgb_expr(0.4, 0.2, 0.1),
			  objs,
			  l,
			  rgb_expr(0.3, 0.3, 0.3)};
    stg.c = test_camera_4;
    stg.s = test_scene_4;
    break;

  case 5 :
    l.direction = v3_expr(-0.58, 0.58, -0.58);
    l.color = rgb_expr(1, 1, 1);
    sphere s = sphere_expr(v3_expr(0, 0, 3), 1, &k_silver, rgb_expr(0.4, 0.4, 0.4));
    o.tag = SPHERE;
    o.o.s = s;
    objs = ol_singleton(o);
    camera test_camera_5 = {v3_expr(0, 0, -5), 300, 400};
    scene test_scene_5 = {rgb_expr(0, 0, 0),
			  objs,
			  l,
			  rgb_expr(0.1, 0.1, 0.1)};
    stg.c = test_camera_5;
    stg.s = test_scene_5;
    break;
  }
  return stg;
}

