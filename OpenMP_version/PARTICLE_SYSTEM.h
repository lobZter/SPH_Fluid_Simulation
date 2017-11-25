#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H
#define _USE_MATH_DEFINES
#include <cmath>

#include <omp.h>
#include "PARTICLE.h"
#include "WALL.h"
#include <vector>
//#include <tr1/tuple>
//#include <map>
#include "FIELD_3D.h"

#include "HW1.h"

#define h 0.0457 //0.02 //0.045

#define GAS_STIFFNESS 3.0 //20.0 // 461.5  // Nm/kg is gas constant of water vapor
#define REST_DENSITY 998.29 // kg/m^3 is rest density of water particle
#define PARTICLE_MASS 0.02 // kg
#define VISCOSITY 3.5 // 5.0 // 0.00089 // Ns/m^2 or Pa*s viscosity of water
#define SURFACE_TENSION 0.0728 // N/m
#define SURFACE_THRESHOLD 7.065
#define KERNEL_PARTICLES 20.0

#define GRAVITY_ACCELERATION -9.80665


#define WALL_K 10000.0 // wall spring constant
#define WALL_DAMPING -55.9 // wall damping constant

#define BOX_SIZE 0.8
#define MAX_PARTICLES 15000

#define INITIAL_SCENARIO SCENARIO_DAM

using namespace std;

class PARTICLE_SYSTEM {

public:
  PARTICLE_SYSTEM();
  ~PARTICLE_SYSTEM();

  void updateGrid();

  // draw to OGL
  void draw();

  void addParticle(const VEC3D& position);

  void addParticle(const VEC3D& position, const VEC3D& velocity);

  void stepVerlet(double dt);

  void stepVerletBrute(double dt);

  void calculateAcceleration();

  void calculateAccelerationBrute();

  void getNeighborParticles(vector<PARTICLE>& totalNeighborParticles, int x, int y, int z);

  void toggleGridVisble();

  void toggleSurfaceVisible();

  void toggleGravity();

  void toggleArrows();

  void toggleTumble();

  void setGravityVectorWithViewVector(VEC3D viewVector);

  void loadScenario(int scenario);

  string fpsf();

  string avgfpsf();

  string particlesf();

  string lowfpsf();

  string peakfpsf();

  FIELD_3D* grid;
  double surfaceThreshold;
  VEC3D gravityVector;

private:
  // list of particles, walls, and springs being simulated
  vector<PARTICLE> _particles;
  vector<WALL>     _walls;

  //unsigned int _particleCount;
  bool _isGridVisible;
  bool _tumble;

  double fps = 0.0;
  double low = 100000.0;
  double peak = 0.0;
  double avgfps = 0.0;
  int particles = 0;

  VEC3D boxSize;

};

#endif
