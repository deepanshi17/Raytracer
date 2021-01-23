#include <cmath>
#include <fstream>

struct Vec3 {
    double x, y, z;
    
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
    Vec3 operator + (const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
    Vec3 operator - (const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
    Vec3 operator * (double d) const { return Vec3(x*d, y*d, z*d); }
    Vec3 operator / (double d) const { return Vec3(x/d, y/d, z/d); }
    
    Vec3 normalize() const {
      double mg = sqrt(x*x + y*y + z*z);
      return Vec3(x/mg,y/mg,z/mg);
    }
  };
inline double dotp(const Vec3& p, const Vec3& q) {
    return(p.x*q.x + p.y*q.y + p.z*q.z);
}

struct Ray {
  Vec3 o,d;
  Ray(const Vec3& o, const Vec3& d) : o(o), d(d) {}
};

struct Sphere {
  Vec3 c;
  double r;
    
  Sphere(const Vec3& c, double r) : c(c), r(r) {}
    
  Vec3 getNormal(const Vec3& pi) const {
      return (pi - c)/r;
  }
    
  bool intersect(const Ray& ray, double &min) const {
    const Vec3 o = ray.o;
    const Vec3 d = ray.d;
    const Vec3 oc = o - c;
    const double b = 2 * dotp(oc, d);
    const double c = dotp(oc, oc) - r*r;
    double disc = b*b - 4 * c;
    if (disc < 1e-4) return false;
    disc = sqrt(disc);
    const double t0 = -b - disc;
    const double t1 = -b + disc;
    min = (t0 < t1) ? t0 : t1;
    return true;
  }
};

int main() {
    // Define height and width of image
    const int H = 500;
    const int W = 500;
    
    const Vec3 black(0, 0, 0);
    const Vec3 white(255, 255, 255);
    const Vec3 red(255, 0, 0);
    const Vec3 green(0, 255, 0);
    const Vec3 blue(0, 0, 255);
    
    const Sphere sphere(Vec3(W*0.5, H*0.5, 50), 50);
    const Sphere light(Vec3(0, 0, 50), 1);
    
    std::ofstream out("render.ppm");
    out << "P3\n" << W << ' ' << H << ' ' << "255\n";
    
    double min;
    Vec3 color(black);
    
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            color = black;
            
            const Ray ray(Vec3(x, y, 0), Vec3(0, 0, 1));
            if(sphere.intersect(ray, min)) {
                const Vec3 pi = ray.o + ray.d*min;
                const Vec3 L = light.c - pi;
                const Vec3 N = sphere.getNormal(pi);
                const double dt = dotp(L.normalize(), N. normalize());
                
                color = (blue + white*dt) * 0.5;
                
                // clamping
                color.x = (color.x > 255) ? 255 : (color.x < 0) ? 0 : color.x;
                color.y = (color.y > 255) ? 255 : (color.y < 0) ? 0 : color.y;
                color.z = (color.z > 255) ? 255 : (color.z < 0) ? 0 : color.z;
                
            }
            out << (int)color.x << ' '
            << (int)color.y << ' '
            << (int)color.z << '\n';
      }    }
    
}
