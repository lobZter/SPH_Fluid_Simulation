#include "PARTICLE_SYSTEM.h"
#include <omp.h>
#define M_PI 3.14159265

#include <time.h>
#include <sstream>
unsigned int iteration = 0;
int scenario;

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
PARTICLE_SYSTEM::PARTICLE_SYSTEM() :
_isGridVisible(false), surfaceThreshold(0.01), gravityVector(0.0,GRAVITY_ACCELERATION,0.0), grid(NULL)
{
  loadScenario(INITIAL_SCENARIO);
}

void PARTICLE_SYSTEM::loadScenario(int newScenario) {
	cout << "Simulating " << PARTICLE::count << " particles" << endl;
  scenario = newScenario;

  // remove all particles
  if (grid)
    delete grid;

  _walls.clear();

  // reset params
  PARTICLE::count = 0;
  particles = 0;
  iteration = 0;

  if (scenario == SCENARIO_DAM) {

    // create long grid
    boxSize.x = BOX_SIZE;
    boxSize.y = BOX_SIZE;
    boxSize.z = BOX_SIZE;

    int gridXRes = (int)ceil(boxSize.x/h);
    int gridYRes = (int)ceil(boxSize.y/h);
    int gridZRes = (int)ceil(boxSize.z/h);

    grid = new FIELD_3D(gridXRes, gridYRes, gridZRes);


    // add walls

    _walls.push_back(WALL(VEC3D(0,0,1), VEC3D(0,0,-boxSize.z/2.0))); // back
    _walls.push_back(WALL(VEC3D(0,0,-1), VEC3D(0,0,boxSize.z/2.0))); // front
    _walls.push_back(WALL(VEC3D(1,0,0), VEC3D(-boxSize.x/2.0,0,0)));     // left
    _walls.push_back(WALL(VEC3D(-1,0,0), VEC3D(boxSize.x/2.0,0,0)));     // right
    _walls.push_back(WALL(VEC3D(0,1,0), VEC3D(0,-boxSize.y/2.0,0))); // bottom

    vector<PARTICLE>& firstGridCell = (*grid)(0,0,0);

    // add particles

    for (double y = -boxSize.y/2.0; y < boxSize.y/2.0; y+= h/2.0) {
      for (double x = -boxSize.x/2.0; x < -boxSize.x/4.0; x += h/2.0) {
        for (double z = -boxSize.z/2.0; z < boxSize.z/2.0; z+= h/2.0) {
          firstGridCell.push_back(PARTICLE(VEC3D(x, y,z)));
        }
      }
    }

    cout << "Loaded dam scenario" << endl;
    cout << "Grid size is " << (*grid).xRes() << "x" << (*grid).yRes() << "x" << (*grid).zRes() << endl;
    cout << "Simulating " << PARTICLE::count << " particles" << endl;
    particles = PARTICLE::count;

  }
  else if (scenario == SCENARIO_CUBE) {

    // create cubed grid

    boxSize.x = BOX_SIZE*2.0;
    boxSize.y = BOX_SIZE;
    boxSize.z = BOX_SIZE/2.0;

    int gridXRes = (int)ceil(boxSize.x/h);
    int gridYRes = (int)ceil(boxSize.y/h);
    int gridZRes = (int)ceil(boxSize.z/h);

    grid = new FIELD_3D(gridXRes, gridYRes, gridZRes);

    // add walls

    _walls.push_back(WALL(VEC3D(0,0,1), VEC3D(0,0,-boxSize.z/2.0))); // back
    _walls.push_back(WALL(VEC3D(0,0,-1), VEC3D(0,0,boxSize.z/2.0))); // front
    _walls.push_back(WALL(VEC3D(1,0,0), VEC3D(-boxSize.x/2.0,0,0))); // left
    _walls.push_back(WALL(VEC3D(-1,0,0), VEC3D(boxSize.x/2.0,0,0))); // right
    _walls.push_back(WALL(VEC3D(0,1,0), VEC3D(0,-boxSize.y/2.0,0))); // bottom

    vector<PARTICLE>& firstGridCell = (*grid)(0,0,0);

    // add particles

    for (double y = 0; y < boxSize.y; y+= h/2.0) {
      for (double x = -boxSize.x/4.0; x < boxSize.x/4.0; x += h/2.0) {
        for (double z = -boxSize.z/4.0; z < boxSize.z/4.0; z+= h/2.0) {
          firstGridCell.push_back(PARTICLE(VEC3D(x,y,z)));
        }
      }
    }

    cout << "Loaded cube scenario" << endl;
    cout << "Grid size is " << (*grid).xRes() << "x" << (*grid).yRes() << "x" << (*grid).zRes() << endl;
    cout << "Simulating " << PARTICLE::count << " particles" << endl;

    particles = PARTICLE::count;
  }
  else if (scenario == SCENARIO_FAUCET) {

    // create cubed grid

    boxSize.x = BOX_SIZE;
    boxSize.y = BOX_SIZE;
    boxSize.z = BOX_SIZE;

    int gridXRes = (int)ceil(boxSize.x/h);
    int gridYRes = (int)ceil(boxSize.y/h);
    int gridZRes = (int)ceil(boxSize.z/h);

    grid = new FIELD_3D(gridXRes, gridYRes, gridZRes);

    // add walls

    _walls.push_back(WALL(VEC3D(0,0,1), VEC3D(0,0,-boxSize.z/2.0))); // back
    _walls.push_back(WALL(VEC3D(0,0,-1), VEC3D(0,0,boxSize.z/2.0))); // front
    _walls.push_back(WALL(VEC3D(1,0,0), VEC3D(-boxSize.x/2.0,0,0))); // left
    _walls.push_back(WALL(VEC3D(-1,0,0), VEC3D(boxSize.x/2.0,0,0))); // right
    _walls.push_back(WALL(VEC3D(0,1,0), VEC3D(0,-boxSize.y/2.0,0))); // bottom

    cout << "Loaded faucet scenario" << endl;
    cout << "Grid size is " << (*grid).xRes() << "x" << (*grid).yRes() << "x" << (*grid).zRes() << endl;
  }

  updateGrid();
}

int rainer = 0;
int rainerg = 2;
int rainerj = 0;

void PARTICLE_SYSTEM::addParticle(const VEC3D& position, const VEC3D& velocity) {
	(*grid)(0,0,0).push_back(PARTICLE(position, velocity));
	particles++;
}

void PARTICLE_SYSTEM::addParticle(const VEC3D& position) {
  addParticle(position, VEC3D());
}

PARTICLE_SYSTEM::~PARTICLE_SYSTEM()
{
  if (grid) delete grid;
}

void PARTICLE_SYSTEM::toggleGridVisble() {
  _isGridVisible = !_isGridVisible;
}

void PARTICLE_SYSTEM::toggleSurfaceVisible() {
  PARTICLE::isSurfaceVisible = !PARTICLE::isSurfaceVisible;
}

void PARTICLE_SYSTEM::toggleTumble() {
  _tumble = !_tumble;
}

void PARTICLE_SYSTEM::toggleGravity() {
  if (gravityVector.magnitude() > 0.0) {
    gravityVector = VEC3D(0,0,0);
  }
  else {
    gravityVector = VEC3D(0,GRAVITY_ACCELERATION,0);
  }
}

void PARTICLE_SYSTEM::toggleArrows() {
  PARTICLE::showArrows = !PARTICLE::showArrows;
}

void PARTICLE_SYSTEM::setGravityVectorWithViewVector(VEC3D viewVector) {

  if (_tumble)
    gravityVector = viewVector * GRAVITY_ACCELERATION;

}

// to update the grid cells particles are located in
// should be called right after particle positions are updated
void PARTICLE_SYSTEM::updateGrid() {
#pragma omp parallel for
  for (int x = 0; x < (*grid).xRes(); x++) {
    for (int y = 0; y < (*grid).yRes(); y++) {
      for (int z = 0; z < (*grid).zRes(); z++) {

        vector<PARTICLE>& particles = (*grid)(x,y,z);

        for (unsigned int p = 0; p < particles.size(); p++) {

          PARTICLE& particle = particles[p];

          int newGridCellX = (int)floor((particles[p].position().x+BOX_SIZE/2.0)/h);
          int newGridCellY = (int)floor((particles[p].position().y+BOX_SIZE/2.0)/h);
          int newGridCellZ = (int)floor((particles[p].position().z+BOX_SIZE/2.0)/h);

          if (newGridCellX < 0)
            newGridCellX = 0;
          else if (newGridCellX >= (*grid).xRes())
            newGridCellX = (*grid).xRes() - 1;
          if (newGridCellY < 0)
            newGridCellY = 0;
          else if (newGridCellY >= (*grid).yRes())
            newGridCellY = (*grid).yRes() - 1;
          if (newGridCellZ < 0)
            newGridCellZ = 0;
          else if (newGridCellZ >= (*grid).zRes())
            newGridCellZ = (*grid).zRes() - 1;

          // check if particle has moved

          if (x != newGridCellX || y != newGridCellY || z != newGridCellZ) {

            // move the particle to the new grid cell
            (*grid)(newGridCellX, newGridCellY, newGridCellZ).push_back(particles[p]);

            // remove it from it's previous grid cell

            particles[p] = particles.back();
            particles.pop_back();
            p--; // important! make sure to redo this index, since a new particle will (probably) be there
          }

        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// OGL drawing
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::draw()
{
  static VEC3F blackColor(0,0,0);
  static VEC3F blueColor(0,0,1);
  static VEC3F whiteColor(1,1,1);
  static VEC3F greyColor(0.2, 0.2, 0.2);
  static VEC3F lightGreyColor(0.8,0.8,0.8);
  //static VEC3F greenColor(34.0 / 255, 139.0 / 255, 34.0 / 255);
  static float shininess = 10.0;

  // draw the particles
  glEnable(GL_LIGHTING);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, blueColor);
  glMaterialfv(GL_FRONT, GL_SPECULAR, whiteColor);
  glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);

  for (int gridCellIndex = 0; gridCellIndex < (*grid).cellCount(); gridCellIndex++) {

    vector<PARTICLE>& particles = (*grid).data()[gridCellIndex];

    for (unsigned int p = 0; p < particles.size(); p++) {

      PARTICLE& particle = particles[p];

      particle.draw();

    }
  }

  glDisable(GL_LIGHTING);

  if (_isGridVisible) {
    glColor3fv(lightGreyColor);

    for (int x = 0; x < grid->xRes(); x++) {
      for (int y = 0; y < grid->yRes(); y++) {
        for (int z = 0; z < grid->zRes(); z++) {
          glPushMatrix();

          glTranslated(x*h-boxSize.x/2.0+h/2.0, y*h-boxSize.y/2.0+h/2.0, z*h-boxSize.z/2.0+h/2.0);
          glutWireCube(h);

          glPopMatrix();
        }
      }
    }
  }

  glColor3fv(greyColor);
  glPopMatrix();
  glScaled(boxSize.x, boxSize.y, boxSize.z);
  glutWireCube(1.0);
  glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
// Verlet integration
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::stepVerlet(double dt)
{

  calculateAcceleration();
#pragma omp parallel for
  for (int gridCellIndex = 0; gridCellIndex < (*grid).cellCount(); gridCellIndex++) {

    vector<PARTICLE>& particles = (*grid).data()[gridCellIndex];

    for (unsigned int p = 0; p < particles.size(); p++) {

      PARTICLE& particle = particles[p];

      VEC3D newPosition = particle.position() + particle.velocity()*dt + particle.acceleration()*dt*dt;
      VEC3D newVelocity = (newPosition - particle.position()) / dt;

      particle.position() = newPosition;
      particle.velocity() = newVelocity;
    }
  }

  if (scenario == SCENARIO_FAUCET && PARTICLE::count < MAX_PARTICLES) {
		VEC3D initialVelocity(-1.4,-1.8,0.4);

		/*addParticle(VEC3D(BOX_SIZE/2.0-h/2.0,BOX_SIZE+h*0.6, 0), initialVelocity);
		  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE, 0), initialVelocity);
		  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*-0.6, 0), initialVelocity);

		  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*0.3, h*0.6), initialVelocity);
		  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*0.3, h*-0.6), initialVelocity);

		  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*-0.3, h*0.6), initialVelocity);
		  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*-0.3, h*-0.6), initialVelocity);*/
		  if (rainer == 0) {
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.1,(BOX_SIZE+h*0.6)*1.5, 0.35), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.1,(BOX_SIZE+h*0.6)*1.5, -0.05), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.25,(BOX_SIZE+h*0.6)*1.5, 0.3), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.25,(BOX_SIZE+h*0.6)*1.5, -0.1), initialVelocity);


			}
			if (rainer == rainerg) {
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.4,(BOX_SIZE+h*0.6)*1.5, 0.25), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.4,(BOX_SIZE+h*0.6)*1.5, -0.15), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.55,(BOX_SIZE+h*0.6)*1.4, 0.2), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.55,(BOX_SIZE+h*0.6)*1.4, -0.15), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.7,(BOX_SIZE+h*0.6)*1.5, 0.15), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.7,(BOX_SIZE+h*0.6)*1.5, -0.2), initialVelocity);

			}
			if (rainer == rainerg*2) {
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.85,(BOX_SIZE+h*0.6)*1.5, 0.1), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.85,(BOX_SIZE+h*0.6)*1.5, -0.25), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*1,(BOX_SIZE+h*0.6)*1.5, 0.05), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-1,(BOX_SIZE+h*0.6)*1.5, -0.3), initialVelocity);

		}
			if (rainer == rainerg*3) {
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.1,(BOX_SIZE+h*0.6)*1.5, -0.35), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.1,(BOX_SIZE+h*0.6)*1.5, 0.05), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.25,(BOX_SIZE+h*0.6)*1.5, -0.3), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.25,(BOX_SIZE+h*0.6)*1.5, 0.1), initialVelocity);


		}
			if (rainer == rainerg*4) {
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.4,(BOX_SIZE+h*0.6)*1.5, -0.25), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.4,(BOX_SIZE+h*0.6)*1.5, 0.15), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.55,(BOX_SIZE+h*0.6)*1.4, -0.2), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.55,(BOX_SIZE+h*0.6)*1.4, 0.15), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.7,(BOX_SIZE+h*0.6)*1.5, -0.15), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.7,(BOX_SIZE+h*0.6)*1.5, 0.2), initialVelocity);

		}
			if (rainer == rainerg*5) {
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*0.85,(BOX_SIZE+h*0.6)*1.5, -0.1), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-0.85,(BOX_SIZE+h*0.6)*1.5, 0.25), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*1,(BOX_SIZE+h*0.6)*1.5, -0.05), initialVelocity);
		  addParticle(VEC3D((BOX_SIZE/2.0-h/2.0)*-1,(BOX_SIZE+h*0.6)*1.5, 0.3), initialVelocity);
			}
			rainer++;
			if (rainer > rainerg*6) rainer = 0;
  }

  updateGrid();

  iteration++;
}

/*
 Calculate the acceleration of each particle using a grid optimized approach.
 For each particle, only particles in the same grid cell and the (26) neighboring grid cells must be considered,
 since any particle beyond a grid cell distance away contributes no force.
*/

void PARTICLE_SYSTEM::calculateAcceleration() {
clock_t t1,t2;
    t1=clock();
	static double hSquared = h*h;
  ///////////////////
  // STEP 1: UPDATE DENSITY & PRESSURE OF EACH PARTICLE
#pragma omp parallel for
  for (int x = 0; x < (*grid).xRes(); x++) {
    for (int y = 0; y < (*grid).yRes(); y++) {
      for (int z = 0; z < (*grid).zRes(); z++) {

        vector<PARTICLE>& particles = (*grid)(x,y,z);

        for (unsigned int p = 0; p < particles.size(); p++) {

          PARTICLE& particle = particles[p];

          particle.density() = 0.0;

          // now iteratate through neighbors

          for (int offsetX = -1; offsetX <= 1 && x+offsetX < (*grid).xRes(); offsetX++) {
            if (x+offsetX < 0) continue;

            for (int offsetY = -1; offsetY <= 1 && y+offsetY < (*grid).yRes(); offsetY++) {
              if (y+offsetY < 0) continue;

              for (int offsetZ = -1; offsetZ <= 1 && z+offsetZ < (*grid).zRes(); offsetZ++) {
                if (z+offsetZ < 0) continue;

                vector<PARTICLE>& neighborGridCellParticles = (*grid)(x+offsetX, y+offsetY, z+offsetZ);

                for (int i = 0; i < neighborGridCellParticles.size(); i++) {

                  VEC3D diffPosition = particle.position() - neighborGridCellParticles[i].position();

                  double radiusSquared = diffPosition.dot(diffPosition);

                  if (radiusSquared <= hSquared) {
					  static double coefficient = 315.0/(64.0*M_PI*pow(h,9));

					  particle.density() += coefficient * pow(hSquared-radiusSquared, 3);
				  }

                }
              }
            }
          }

          particle.density() *= PARTICLE_MASS;
          particle.pressure() = GAS_STIFFNESS * (particle.density() - REST_DENSITY);
        }
      }
    }
  }

  ///////////////////
  // STEP 2: COMPUTE FORCES FOR ALL PARTICLES
#pragma omp parallel for
  for (int x = 0; x < (*grid).xRes(); x++) {
    for (int y = 0; y < (*grid).yRes(); y++) {
      for (int z = 0; z < (*grid).zRes(); z++) {

        vector<PARTICLE>& particles = (*grid)(x,y,z);

        for (unsigned int p = 0; p < particles.size(); p++) {

          PARTICLE& particle = particles[p];

          //cout << "particle id: " << particle.id() << endl;

          VEC3D f_pressure,
          f_viscosity,
          f_surface,
          f_gravity = gravityVector * particle.density(),
          colorFieldNormal;

          double colorFieldLaplacian = 0.0;

          // now iteratate through neighbors

          for (int offsetX = -1; offsetX <= 1 && x+offsetX < (*grid).xRes(); offsetX++) {
            if (x+offsetX < 0) continue;

            for (int offsetY = -1; offsetY <= 1 && y+offsetY < (*grid).yRes(); offsetY++) {
              if (y+offsetY < 0) continue;

              for (int offsetZ = -1; offsetZ <= 1 && z+offsetZ < (*grid).zRes(); offsetZ++) {
                if (z+offsetZ < 0) continue;

                vector<PARTICLE>& neighborGridCellParticles = (*grid)(x+offsetX, y+offsetY, z+offsetZ);

                for (unsigned int i = 0; i < neighborGridCellParticles.size(); i++) {

                  PARTICLE& neighbor = neighborGridCellParticles[i];

                  //if (particle.id() == neighbor.id()) continue; // SKIPPING COMPARISON OF THE SAME PARTICLE

                  VEC3D diffPosition = particle.position() - neighbor.position();
                  double radiusSquared = diffPosition.dot(diffPosition);

                  if (radiusSquared <= hSquared) {
					static double coefficient2 = -945.0/(32.0*M_PI*pow(h,9));
					static double coefficient = -45.0/(M_PI*pow(h,6));

					double radius = sqrt(radiusSquared);

                    if (particle.id() != neighbor.id()) {
                      f_pressure += (particle.pressure()/pow(particle.density(),2)+neighbor.pressure()/pow(neighbor.density(),2))*(coefficient * pow(h-radius, 2) * diffPosition / radius);
                      f_viscosity += (neighbor.velocity() - particle.velocity()) * -1 * coefficient * (h - radius) / neighbor.density();

                    }

                    colorFieldNormal += coefficient2 * pow(hSquared-radiusSquared, 2) * diffPosition / neighbor.density();
                    colorFieldLaplacian += coefficient2 * (hSquared-radiusSquared) * (3.0*hSquared - 7.0*radiusSquared) / neighbor.density();
                  }
                }
              }
            }
          } // end of neighbor grid cell iteration

          f_pressure *= -PARTICLE_MASS * particle.density();
          f_viscosity *= VISCOSITY * PARTICLE_MASS;
          colorFieldNormal *= PARTICLE_MASS;
          particle.normal = -1.0 * colorFieldNormal;
          colorFieldLaplacian *= PARTICLE_MASS;

          // surface tension force

          double colorFieldNormalMagnitude = colorFieldNormal.magnitude();

          if (colorFieldNormalMagnitude > SURFACE_THRESHOLD) {

            particle.flag() = true;
            f_surface = -SURFACE_TENSION * colorFieldNormal / colorFieldNormalMagnitude * colorFieldLaplacian;

          }
          else {
            particle.flag() = false;
          }

          // ADD IN SPH FORCES
          particle.acceleration() = (f_pressure + f_viscosity + f_surface + f_gravity) / particle.density();

          // EXTERNAL FORCES HERE (USER INTERACTION, SWIRL)
          for (unsigned int i = 0; i < _walls.size(); i++) {

			double d = (_walls[i].point() - particle.position()).dot(_walls[i].normal()) + 0.01; // particle radius

			if (d > 0.0) {
			  particle.acceleration() += WALL_K * _walls[i].normal() * d;
			  particle.acceleration() += WALL_DAMPING * particle.velocity().dot(_walls[i].normal()) * _walls[i].normal();
			}
		  }
        }
      }
    }
  }

t2=clock();
    float diff = 1.0/(((float)t2-(float)t1)/CLOCKS_PER_SEC);

    if (avgfps != 0.0)
        if (diff > 500)
            avgfps += 500;
        else
            avgfps += diff;
    else
        if (diff > 500)
            avgfps = 1000;
        else
            avgfps = diff + diff;

    fps = diff;
    if (fps > peak)
        if (fps > 500)
            peak = 500;
        else
            peak = fps;
    if (fps < low)
        low = fps;
    avgfps = avgfps / 2.0;
}

string PARTICLE_SYSTEM::particlesf() {
    std::ostringstream sstream;
    sstream << particles;
    return sstream.str();
}

string PARTICLE_SYSTEM::fpsf() {
    std::ostringstream sstream;
    sstream << fps;
    return sstream.str();
}

string PARTICLE_SYSTEM::avgfpsf() {
    std::ostringstream sstream;
    sstream << avgfps;
    return sstream.str();
}

string PARTICLE_SYSTEM::lowfpsf() {
    std::ostringstream sstream;
    double t = low;
    if (low == 100000.0)
        t = 0.0;
    sstream << t;
    return sstream.str();
}

string PARTICLE_SYSTEM::peakfpsf() {
    std::ostringstream sstream;
    sstream << peak;
    return sstream.str();
}

