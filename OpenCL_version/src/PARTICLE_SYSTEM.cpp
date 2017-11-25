#include "PARTICLE_SYSTEM.h"
using namespace std;

int read_kernel_from_file(char *filename, char **source, size_t *len);
char filename[] = "kernel.c";
char *source[1];
size_t src_len[1];

PARTICLE_SYSTEM::PARTICLE_SYSTEM() : 
_isGridVisible(false), surfaceThreshold(0.01), gravityVector(0.0,GRAVITY_ACCELERATION,0.0), grid(NULL)
{
	initOpenCL();
	loadScenario(INITIAL_SCENARIO);
}

void PARTICLE_SYSTEM::initOpenCL() {

	printf("Initialize OpenCL object and context\n");
	// get first platfrom
	err = clGetPlatformIDs(1, &platform, NULL);
	if (err != CL_SUCCESS) {
		std::cerr << "Unable to get platforms\n";
	}
	// get first GPU device on platform
	err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
	if (err != CL_SUCCESS) {
		std::cerr << "Unable to get devices\n";
	}
	// create context for selected device
	context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
	if (err != CL_SUCCESS) {
		std::cerr << "Can't create OpenCL context\n";
	}
	// create commadqueue
	commandqueue = clCreateCommandQueue(context, device, NULL, &err);
	if (err != CL_SUCCESS) {
		std::cerr << "Can't create command queue\n";
		clReleaseContext(context);
	}
	// create program
	err = read_kernel_from_file(filename, &source[0], &src_len[0]);
	if (err != CL_SUCCESS) {
		printf("read_kernel_from_file() failed: %d\n", err);
		clReleaseContext(context);
	}

	program = clCreateProgramWithSource(context, 1, (const char **)source, (size_t *)src_len, &err);
	if (!program || err) {
		printf("clCreateProgramWithSource() failed: %d\n", err);
	}
	free(source[0]);
	//std::cout << source << std::endl;
	/*program = clCreateProgramWithSource(context, 1, (const char **)&kernel_src, NULL, &err);
	if (err != CL_SUCCESS) {
	std::cerr << "Can't create program\n";
	clReleaseCommandQueue(commandqueue);
	clReleaseContext(context);
	}*/
	// build program for all devices associated with program
	err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
	if (err != CL_SUCCESS) {
		std::cerr << "Can't build program\n" << err << std::endl;
		if (err == CL_BUILD_PROGRAM_FAILURE) {
			// Determine the size of the log
			size_t log_size;
			clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
			// Allocate memory for the log
			char *log = (char *)malloc(log_size);
			// Get the log
			clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
			// Print the log
			printf("%s\n", log);
		}
		clReleaseCommandQueue(commandqueue);
		clReleaseContext(context);
	}
	// use func "histogram" as the kernel
	kernel = clCreateKernel(program, "stepVerlet", &err);
	if (err != CL_SUCCESS) {
		std::cerr << "Can't load kernel\n";
		clReleaseProgram(program);
		clReleaseCommandQueue(commandqueue);
		clReleaseContext(context);
	}
}

void PARTICLE_SYSTEM::loadScenario(int newScenario) {
    
    boxSize.x = BOX_SIZE*2.0;
    boxSize.y = BOX_SIZE;
    boxSize.z = BOX_SIZE/2.0;
    
	walls[0].normal.s0 = 0.0; walls[0].normal.s1 = 0.0; walls[0].normal.s2 = 1.0;
	walls[1].normal.s0 = 0.0; walls[1].normal.s1 = 0.0; walls[1].normal.s2 = -1.0;
	walls[2].normal.s0 = 1.0; walls[2].normal.s1 = 0.0; walls[2].normal.s2 = 0.0;
	walls[3].normal.s0 = -1.0; walls[3].normal.s1 = 0.0; walls[3].normal.s2 = 0.0;
	walls[4].normal.s0 = 0.0; walls[4].normal.s1 = 1.0; walls[4].normal.s2 = 0.0;

	walls[0].point.s0 = 0.0; walls[0].point.s1 = 0.0; walls[0].point.s2 = -boxSize.z/2;
	walls[1].point.s0 = 0.0; walls[1].point.s1 = 0.0; walls[1].point.s2 = boxSize.z/2;
	walls[2].point.s0 = -boxSize.x/2; walls[2].point.s1 = 0.0; walls[2].point.s2 = 0.0;
	walls[3].point.s0 = boxSize.x/2; walls[3].point.s1 = 0.0; walls[3].point.s2 = 0.0;
	walls[4].point.s0 = 0.0; walls[4].point.s1 = -boxSize.y/2; walls[4].point.s2 = 0.0;

	printf("%f %f %f	%f %f %f \n", walls[0].normal.s0, walls[0].normal.s1, walls[0].normal.s2, walls[0].point.s0, walls[0].point.s1, walls[0].point.s2);
	printf("%f %f %f	%f %f %f \n", walls[1].normal.s0, walls[1].normal.s1, walls[1].normal.s2, walls[1].point.s0, walls[1].point.s1, walls[1].point.s2);
	printf("%f %f %f	%f %f %f \n", walls[2].normal.s0, walls[2].normal.s1, walls[2].normal.s2, walls[2].point.s0, walls[2].point.s1, walls[2].point.s2);
	printf("%f %f %f	%f %f %f \n", walls[3].normal.s0, walls[3].normal.s1, walls[3].normal.s2, walls[3].point.s0, walls[3].point.s1, walls[3].point.s2);
	printf("%f %f %f	%f %f %f \n", walls[4].normal.s0, walls[4].normal.s1, walls[4].normal.s2, walls[4].point.s0, walls[4].point.s1, walls[4].point.s2);

	memset(particles, 0x00, sizeof(particles));
	int particleIdx = 0;
    for (double y = -boxSize.y/2.0; y < boxSize.y/2.0; y+= h/2.0) {
      for (double x = -boxSize.x/2.0; x < -boxSize.x/4.0; x += h/2.0) {
        for (double z = -boxSize.z/2.0; z < boxSize.z/2.0; z+= h/2.0) {
			particles[particleIdx].position.s0 = x;
			particles[particleIdx].position.s1 = y;
			particles[particleIdx].position.s2 = z;
			particleIdx++;
        }
      }
    }
  
	// create buffer in gpu ram
	particles_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(particles), NULL, &err);
	if (err != CL_SUCCESS) {
		std::cerr << "Can't create OpenCL buffer\n";
		clReleaseMemObject(particles_d);
		clReleaseCommandQueue(commandqueue);
		clReleaseContext(context);
	}
	walls_d = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(walls), NULL, &err);
	if (err != CL_SUCCESS) {
		std::cerr << "Can't create OpenCL buffer\n";
		clReleaseMemObject(particles_d);
		clReleaseMemObject(walls_d);
		clReleaseCommandQueue(commandqueue);
		clReleaseContext(context);
	}

	// init walls
	err = clEnqueueWriteBuffer(commandqueue, walls_d, CL_TRUE, 0, sizeof(walls), walls, NULL, NULL, NULL);
	if (err != CL_SUCCESS) {
		std::cerr << "Write OpenCL buffer error\n";
		clReleaseMemObject(particles_d);
		clReleaseMemObject(walls_d);
		clReleaseCommandQueue(commandqueue);
		clReleaseContext(context);
	}
	err = clEnqueueWriteBuffer(commandqueue, particles_d, CL_TRUE, 0, sizeof(particles), particles, NULL, NULL, NULL);
}

void PARTICLE_SYSTEM::draw() 
{ 
  static VEC3F blackColor(0,0,0); 
  static VEC3F blueColor(0,0,1); 
  static VEC3F whiteColor(1,1,1);
  static VEC3F greyColor(0.2, 0.2, 0.2);
  static VEC3F lightGreyColor(0.8,0.8,0.8);
  static VEC3F greenColor(34.0 / 255, 139.0 / 255, 34.0 / 255);
  static VEC3F lightBlueColor(0.01, 0.25, 1.0);
  static float shininess = 10.0;

  // draw the particles
  glEnable(GL_LIGHTING);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, lightBlueColor);
  glMaterialfv(GL_FRONT, GL_SPECULAR, whiteColor);
  glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);
  
  
  for (int i = 0; i < PARTICLE_NUM; i++) {
	  glPushMatrix();
	  glTranslated(particles[i].position.s0, particles[i].position.s1, particles[i].position.s2);
	  glutSolidSphere(0.015, 10, 10);
	  glPopMatrix();
  }
    
  glMaterialfv(GL_FRONT, GL_DIFFUSE, blueColor);
  glColor3fv(greyColor);
  
  glPopMatrix();
  glScaled(boxSize.x, boxSize.y, boxSize.z);
  glutWireCube(1.0);
  glPopMatrix();
}

void PARTICLE_SYSTEM::stepVerlet(double dt)
{
	clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&particles_d);
	clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&walls_d);

	

	global_work_size = PARTICLE_NUM;
	err = clEnqueueNDRangeKernel(commandqueue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	if (err != CL_SUCCESS) {
		printf("clEnqueueNDRangeKernel() failed: %d\n", err);
	}

	err = clEnqueueReadBuffer(commandqueue, particles_d, CL_TRUE, 0, sizeof(particles), particles, NULL, NULL, NULL);
	if (err != CL_SUCCESS) {
		std::cerr << "Can't read back data\n";
	}
}

int read_kernel_from_file(char *filename, char **source, size_t *len) {
	//need to put the kernel function into a file named xxx.cl
	struct stat statbuf;
	FILE *fp;
	size_t file_len;

	fp = fopen(filename, "r");
	if (fp == 0) return -1; //failed reading file

	stat(filename, &statbuf);
	file_len = (size_t)statbuf.st_size;
	*len = file_len;
	*source = (char *)malloc(file_len + 1);
	fread(*source, file_len, 1, fp);
	(*source)[file_len] = '\0';

	fclose(fp);
	return 0;
}

PARTICLE_SYSTEM::~PARTICLE_SYSTEM()
{
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
		gravityVector = VEC3D(0, 0, 0);
	}
	else {
		gravityVector = VEC3D(0, GRAVITY_ACCELERATION, 0);
	}
}

void PARTICLE_SYSTEM::toggleArrows() {
	PARTICLE::showArrows = !PARTICLE::showArrows;
}

void PARTICLE_SYSTEM::setGravityVectorWithViewVector(VEC3D viewVector) {

	if (_tumble)
		gravityVector = viewVector * GRAVITY_ACCELERATION;

}