# Makefile (starter)

# CS152 Winter 2015
# Project A
# Feb 2015

.PHONY: clean

SRC = v3.c rgb.c ray.c sphere.c poster.c object.c obj_list.c light.c scene.c camera.c hit_test.c stage.c intersect.c raytracer.c rt.c

trace: raytracer.h $(SRC)
	gcc -O3 -Wall -o rt $(SRC) -lm

clean:
	rm -rf rt *.dSYM *.ppm
