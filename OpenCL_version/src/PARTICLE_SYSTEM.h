#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <vector>
#include <CL/cl.h>
#include "PARTICLE.h"
#include "WALL.h"
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
#define WALL_DAMPING -0.9 // wall damping constant

#define BOX_SIZE 0.4
#define MAX_PARTICLES 3000

#define INITIAL_SCENARIO SCENARIO_DAM

#define PARTICLE_NUM 1458

typedef struct vec3 {
	double s0, s1, s2;
}VEC;

typedef struct particle {
	VEC gridIdx;
	VEC position;
	VEC velocity;
	VEC acceleration;
	double density;
	double pressure;
	int isSurface;
} PARTICLE_D;

typedef struct wall {
	VEC normal;
	VEC point;
} WALL_D;

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
  
  void collisionForce(PARTICLE& particle, VEC3D& f_collision);
  
  void getNeighborParticles(vector<PARTICLE>& totalNeighborParticles, int x, int y, int z);
  
  double Wpoly6(double radiusSquared);
  
  void Wpoly6Gradient(VEC3D& diffPosition, double radiusSquared, VEC3D& gradient);
  
  double Wpoly6Laplacian(double radiusSquared); 
  
  void WspikyGradient(VEC3D& diffPosition, double radiusSquared, VEC3D& gradient);
  
  double WviscosityLaplacian(double radiusSquared);
  
  void toggleGridVisble();
  
  void toggleSurfaceVisible();
  
  void toggleGravity();
  
  void toggleArrows();
  
  void toggleTumble();
  
  void generateFaucetParticleSet();
  
  void setGravityVectorWithViewVector(VEC3D viewVector);
  
  
  void loadScenario(int scenario);
  
  void initOpenCL();
  
  cl_int err;
  cl_platform_id platform;
  cl_device_id device;
  cl_context context;
  cl_command_queue commandqueue;
  cl_program program;
  cl_kernel kernel;
  cl_mem particles_d, walls_d;
  size_t local_work_size;
  size_t global_work_size;

  WALL_D walls[5];
  PARTICLE_D particles[PARTICLE_NUM];
  
  
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
  
  VEC3D boxSize;

};

#endif
