#include "raytracer.h"

int within_poster(v3 point, poster s)
{
  double xmin = s.upper_left.x;
  double xmax = s.upper_left.x + s.w;
  double ymin = s.upper_left.y - s.h;
  double ymax = s.upper_left.y;
  return ((xmin <= point.x) && (point.x <= xmax) && (ymin <= point.y) && (point.y <= ymax));
}

int on_cylinder (v3 point, cylinder c)
{
  double ymin = c.center.y;
  double ymax = c.center.y + c.height;
  return ((ymin <= point.y) && (point.y <= ymax));
}

hit_test intersect_sphere (ray r, sphere s)
{
  hit_test result = {0};
  v3 sc = s.center;
  double sr = s.radius;
  v3 a = v3_sub(r.origin, sc);
  double b = v3_dot(a, r.direction);
  double c = v3_dot(a, a) - (sr * sr);
  double d = (b * b) - c;
  double t = -b - sqrt(d);
  result.miss = MISS;
  if (d > 0 && t > 0) {
    v3 intersection = ray_position(r, t);  
    result.miss = HIT;
    result.t = t;
    result.hit_point = intersection;
    result.surf = s.surface_color(sc, intersection);
    result.shine = s.shine;
    result.surf_norm = v3_norm(v3_sub(intersection, sc));
  }
  return result;
}

hit_test intersect_poster (ray r, poster p)
{
  hit_test result = {0};
  v3 ro = r.origin;
  v3 rd = r.direction;
  v3 n = v3_expr(0, 0, -1);
  double d = p.upper_left.z;
  double t = -(v3_dot(ro, n) + d) / v3_dot(rd, n);
  result.miss = MISS;
  v3 intersection = ray_position(r, t);
  if (t>0 && within_poster(intersection, p)) {
    result.miss = HIT;
    result.t = t;
    result.hit_point = intersection;
    result.surf = p.surface_color(p.upper_left, intersection);
    result.shine = p.shine;
    result.surf_norm = n;
  }
  return result;
}

hit_test intersect_cylinder (ray r, cylinder c) {
  hit_test result;
  v3 rp = v3_expr(r.direction.x,0,r.direction.z);
  double mp = v3_mag(rp);
  v3 np = v3_norm(rp);
  double xbar = r.origin.x - c.center.x;
  double zbar = r.origin.z - c.center.z;
  double a = pow(np.x,2) + pow(np.z,2);
  double b = 2 * ( (xbar*np.x) + (zbar*np.z) );
  double c1 = pow(xbar,2) + pow(zbar,2) - pow(c.radius,2);
  double d = pow(b,2) - (4*a*c1);
  result.miss = MISS;
  if(d >= 0) {
    double t_front = (-b - sqrt(d)) / (2*a);
    double t_back = (-b + sqrt(d)) / (2*a);
    v3 p_front = ray_position(r, t_front/mp);
    v3 p_back = ray_position(r, t_back/mp);
    if( (t_front < t_back) &&
	(t_front > 0) &&
	(p_front.y >= c.center.y) && (on_cylinder(p_front, c))) {
      result.miss = HIT;
      result.t = (t_front/mp);
      result.hit_point = p_front;
      result.surf = c.surface_color(c.center, p_front);
      result.shine = c.shine;
      v3 c2 = v3_expr(c.center.x, p_front.y, c.center.z);
      result.surf_norm = v3_norm(v3_sub(p_front,c2));
    } else if( (t_back > 0) &&
	       (p_back.y >= c.center.y) &&
	       (on_cylinder(p_back, c)) ) {
      result.miss = HIT;
      result.t = (t_back/mp);
      result.hit_point = p_back;
      result.surf = c.surface_color(c.center, p_back);
      result.shine = c.shine;
      v3 c3 = v3_expr(c.center.x, p_back.y, c.center.z);
      result.surf_norm = v3_norm(v3_sub(c3,p_back));
    }
  }
  return result;
}


hit_test intersect (ray r, object o)
{
  switch (o.tag) {
  case SPHERE :
    return intersect_sphere(r, o.o.s);
  case POSTER :
    return intersect_poster(r, o.o.p);
  case CYLINDER :
    return intersect_cylinder(r, o.o.cyl);
  }
  fprintf(stderr, "invalid object tag\n");
  exit(1);
}
