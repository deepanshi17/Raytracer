#ifndef __RAYTRACER_H__
#define __RAYTRACER_H__

/* raytracer.h
 * CS152 Winter 2015
 * Project A
 * Feb 2015
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* ====================================== */ 
/* ====================================== */ 
/* === STRUCTURE AND TYPE DEFINITIONS === */ 
/* ====================================== */ 
/* ====================================== */ 

typedef unsigned char byte;

typedef struct {
  double x;
  double y;
  double z;
} v3;

typedef struct {
  /* each value must be on interval [0,1] */
  double r;
  double g;
  double b;
} rgb;

typedef struct {
  v3 origin;
  v3 direction;
} ray;

typedef struct {
  v3 center;
  double radius;
  rgb(*surface_color)(v3 center, v3 loc);
  rgb shine;
} sphere;

typedef struct {
  v3 upper_left;
  double w;
  double h;
  rgb(*surface_color)(v3 upper_left, v3 loc);
  rgb shine;
} poster;

typedef struct {
  v3 center;
  double radius;
  double height;
  rgb(*surface_color)(v3 center, v3 loc);
  rgb shine;
} cylinder;

union object_union {
  sphere s;
  poster p;
  cylinder cyl;
};

typedef struct {
  enum  { SPHERE, POSTER, CYLINDER } tag;
  union object_union o;
} object;

typedef struct obj_list obj_list;
/* convention: NULL is the empty object list */
struct obj_list {
  object obj;
  obj_list *next;
};

typedef struct {
  v3 direction; /* note: the vector points AT the light */
  rgb color;
} light;

typedef struct {
  rgb bg;
  obj_list *objects;
  light light;
  rgb amb;
} scene;

typedef struct {
  enum {HIT=0, MISS=1} miss;
  double t;
  v3 hit_point;
  rgb surf;
  rgb shine;
  v3 surf_norm;
} hit_test;

typedef struct {
  v3 loc;
  unsigned int h;
  unsigned int w;
} camera;

typedef struct {
  camera c;
  scene s;
} stage;

/* ================== */ 
/* ================== */ 
/* === OPERATIONS === */ 
/* ================== */ 
/* ================== */ 

/* Please note: the following is not necessarily every operation you
 * will need, nor will you necessarily use every operation in your
 * implementation.  Add and implement additional operations as you
 * require them.
 */

/* ================================== */ 
/* === vector operations => v3.c === */ 
/* ================================== */ 

v3 v3_expr(double x, double y, double z);

/* v3_add : add two vectors */
/* (x,y,z) + (X,Y,Z) = (x+X,y+Y,z+Z) */
v3 v3_add(v3 v1, v3 v2);

/* v3_sub : subtract the second vector from the first */
/* (x,y,z) - (X,Y,Z) = (x-X,y-Y,z-Z) */
v3 v3_sub(v3 v1, v3 v2);

/* v3_neg : negate a vector */
/* -(x,y,z) = (-x,-y,-z) */
v3 v3_neg(v3 v);

/* v3_scale : scale the components of a vector by scalar s */
/* s(x,y,z) = (sx,sy,sz) */
v3 v3_scale(double s, v3 v);

/* v3_dot : compute the dot product of two vectors */
/* (x,y,z) dot (X,Y,Z) = xX+yY+zZ */
double v3_dot(v3 v1, v3 v2);

/* v3_mag : compute the magnitude of a vector */
/* the square root of the sum of the squares of the components */
double v3_mag(v3 v);

/* v3_norm : normalize the vector */
/* i.e., a vector pointing in the same direction whose magnitude is 1 */
/* - normalize a vector by scaling by 1 over its magnitude */
/* - if a vector's magnitude is zero, it can't be normalized; return it */
v3 v3_norm(v3 v);

char *v3_tos(v3 v);

void v3_show(FILE *f, v3 v);

/* =============================== */ 
/* === rgb operations => rgb.c === */ 
/* =============================== */ 

/* reject if any component is outside [0.0,1.0] */
rgb rgb_expr(double r, double g, double b);

/* rgb_modulate : multiply the components or two rgbs */
/* (r,g,b) * (R,G,B) = (rR,gG,bB) */
rgb rgb_modulate(rgb c1, rgb c2);

/* rgb_scale : scale each component by given scalar, but */
/* don't allow any component in the result to fall outside [0,1] */
rgb rgb_scale(double s, rgb c);

/* add components of colors, but */
/* don't allow any component in the result to exceed 1 */
rgb rgb_add(rgb c1, rgb c2);

char *rgb_tos(rgb c);

void rgb_show(FILE *f, rgb c);

/* rgb_show_bytes */
/* print three integers on [0,255] with spaces in between and a newline */
/* suitable for use in a PPM */
void rgb_show_bytes(FILE *f, rgb c);

/* =============================== */ 
/* === ray operations => ray.c === */ 
/* =============================== */ 

/* see web page for ray_position code */
v3 ray_position(ray r, double t);

char *ray_tos(ray r);

void ray_show(FILE *f, ray r);

/* ===================================== */ 
/* === sphere operations => sphere.c === */ 
/* ===================================== */ 

sphere sphere_expr(v3 v, double r, rgb(*sc)(v3,v3), rgb sh);

char *sphere_tos(sphere s);

void sphere_show(FILE *f, sphere s);

/* ===================================== */ 
/* === slab operations => poster.c ===== */ 
/* ===================================== */ 

poster poster_expr(v3 ulc, double w, double h, rgb(*sc)(v3,v3), rgb sh);

char *poster_tos(poster s);

void poster_show(FILE *f, poster s);

/* ===================================== */ 
/* === object operations => object.c === */ 
/* ===================================== */ 

char *object_tos(object o);

void object_show(FILE *f, object o);

/* ======================================== */ 
/* === object list ops => obj_list.c === */ 
/* ======================================== */ 

obj_list *ol_cons(object o, obj_list *os);

obj_list *ol_singleton(object s);

char *ol_tos(obj_list *os);

void ol_show(FILE *f, obj_list *os);

void ol_free(obj_list *ss);

/* =================================== */ 
/* === light operations => light.c === */ 
/* =================================== */ 

char *light_tos(light l);

void light_show(FILE *f, light l);

/* =================================== */ 
/* === scene operations => scene.c === */ 
/* =================================== */ 

char *scene_tos(scene sc);

void scene_show(FILE *f, scene sc);

/* ========================================= */ 
/* === hit_test operations => hit_test.c === */ 
/* ========================================= */ 

char *hit_test_tos(hit_test t);

void hit_test_show(FILE *f, hit_test t);

/* ===================================== */ 
/* === camera operations => camera.c === */ 
/* ===================================== */ 

char *camera_tos(camera cm);

void camera_show(FILE *f, camera cm);

/* =================================== */ 
/* === stage operations => stage.c === */ 
/* =================================== */ 

char *stage_tos(stage stg);

void stage_show(FILE *f, stage stg);

/* The function get_stage just needs to return a valid stage. */
/* The argument may be used, or not, at your discretion. */
stage get_stage(unsigned int selection);

/* ================================================ */ 
/* ================================================ */ 
/* === RAY TRACER LOGIC and MAIN => raytracer.c === */ 
/* ================================================ */ 
/* ================================================ */ 

/* === project 1 operations === */

/* === intersect.c === */

hit_test intersect(ray r, object obj);

/* === raytracer.c === */

int shadow(v3 p, light l, obj_list *objs);

rgb lighting(scene s, ray r, hit_test h);

v3 logical_loc(camera c, unsigned int x, unsigned int y);

rgb trace_ray(scene s, ray r);

void render(FILE *f, stage g);

/* === rt.c === */

int main(int argc, char *argv[]);

#endif /* __RAYTRACER_H__ */
