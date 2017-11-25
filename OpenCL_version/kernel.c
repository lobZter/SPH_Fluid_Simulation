#define PARTICLE_NUM 1458
#define WALL_NUM 5

#define REST_DENSITY 998.29 // kg/m^3
#define PARTICLE_MASS 0.02 // kg
#define VISCOSITY 3.5 // Ns/m^2 or Pa*s
#define SURFACE_TENSION 0.0728 // N/m 
#define SURFACE_THRESHOLD 7.065
#define GAS_STIFFNESS 3.0 // Nm/kg
#define KERNEL_PARTICLES 20.0
#define h 0.0457 // m
#define GRAVITY_ACCELERATION -9.80665

#define WALL_K 10000.0 // wall spring constant
#define WALL_DAMPING -0.9 // wall damping constant

#define dt 0.003

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
} PARTICLE;

typedef struct wall {
	VEC normal;
	VEC point;
} WALL;

double Wpoly6(double radiusSquared)
{
	double coefficient = 315.0 / (64.0*M_PI*pow(h, 9));
	double hSquared = h*h;

	return coefficient * pow(hSquared - radiusSquared, 3);
}

VEC Wpoly6Gradient(VEC diffPosition, double radiusSquared)
{
	double coefficient = -945.0 / (32.0*M_PI*pow(h, 9));
	double hSquared = h*h;

	VEC result;
	result.s0 = diffPosition.s0 * coefficient * pow(hSquared - radiusSquared, 2);
	result.s1 = diffPosition.s1 * coefficient * pow(hSquared - radiusSquared, 2);
	result.s2 = diffPosition.s2 * coefficient * pow(hSquared - radiusSquared, 2);

	return result;
}

double Wpoly6Laplacian(double radiusSquared)
{
	double coefficient = -945.0 / (32.0*M_PI*pow(h, 9));
	double hSquared = h*h;

	return coefficient * (hSquared - radiusSquared) * (3.0*hSquared - 7.0*radiusSquared);
}

VEC WspikyGradient(VEC diffPosition, double radiusSquared)
{
	double coefficient = -45.0 / (M_PI*pow(h, 6));
	double radius = sqrt(radiusSquared);

	VEC result;
	result.s0 = coefficient * pow(h - radius, 2) * diffPosition.s0 / radius;
	result.s1 = coefficient * pow(h - radius, 2) * diffPosition.s1 / radius;
	result.s2 = coefficient * pow(h - radius, 2) * diffPosition.s2 / radius;

	return result;
}

double WviscosityLaplacian(double radiusSquared)
{
	double coefficient = 45.0 / (M_PI*pow(h, 6));
	double radius = sqrt(radiusSquared);

	return coefficient * (h - radius);
}

double dotProduct(VEC a, VEC b) {
	return a.s0 * b.s0 + a.s1 * b.s1 + a.s2 * b.s2;
}


__kernel void stepVerlet(
	__global PARTICLE particles[PARTICLE_NUM],
	__global WALL walls[WALL_NUM]) 
{
	int i, neighborIdx, particleIdx = get_global_id(0);

	//printf("%f %f %f	%f %f %f \n", walls[0].normal.s0, walls[0].normal.s1, walls[0].normal.s2, walls[0].point.s0, walls[0].point.s1, walls[0].point.s2);
	//printf("%f %f %f	%f %f %f \n", walls[1].normal.s0, walls[1].normal.s1, walls[1].normal.s2, walls[1].point.s0, walls[1].point.s1, walls[1].point.s2);
	//printf("%f %f %f	%f %f %f \n", walls[2].normal.s0, walls[2].normal.s1, walls[2].normal.s2, walls[2].point.s0, walls[2].point.s1, walls[2].point.s2);
	//printf("%f %f %f	%f %f %f \n", walls[3].normal.s0, walls[3].normal.s1, walls[3].normal.s2, walls[3].point.s0, walls[3].point.s1, walls[3].point.s2);
	//printf("%f %f %f	%f %f %f \n", walls[4].normal.s0, walls[4].normal.s1, walls[4].normal.s2, walls[4].point.s0, walls[4].point.s1, walls[4].point.s2);

	// STEP 1: UPDATE DENSITY & PRESSURE OF EACH PARTICLE

	particles[particleIdx].density = 0.0;

	for (neighborIdx = 0; neighborIdx < PARTICLE_NUM; neighborIdx++)
	{
		if (neighborIdx == particleIdx) continue;

		VEC diffPosition;
		diffPosition.s0 = particles[particleIdx].position.s0 - particles[neighborIdx].position.s0;
		diffPosition.s1 = particles[particleIdx].position.s1 - particles[neighborIdx].position.s1;
		diffPosition.s2 = particles[particleIdx].position.s2 - particles[neighborIdx].position.s2;

		double radiusSquared = dotProduct(diffPosition, diffPosition);

		if (radiusSquared <= h*h)
			particles[particleIdx].density += Wpoly6(radiusSquared);
	
	}

	particles[particleIdx].density *= PARTICLE_MASS;
	particles[particleIdx].pressure = GAS_STIFFNESS * (particles[particleIdx].density - REST_DENSITY);

	barrier(CLK_GLOBAL_MEM_FENCE);
	// STEP 2: COMPUTE FORCES FOR ALL PARTICLES

	VEC f_pressure; f_pressure.s0 = 0; f_pressure.s1 = 0; f_pressure.s2 = 0;
	VEC	f_viscosity; f_viscosity.s0 = 0; f_viscosity.s1 = 0; f_viscosity.s2 = 0;
	VEC	f_surface; f_surface.s0 = 0; f_surface.s1 = 0; f_surface.s2 = 0;
	VEC	f_gravity;  f_gravity.s0 = 0; f_gravity.s1 = particles[particleIdx].density * -9.80665; f_gravity.s2 = 0;
	VEC	colorFieldNormal; colorFieldNormal.s0 = 0; colorFieldNormal.s1 = 0; colorFieldNormal.s2 = 0;
	double colorFieldLaplacian = 0;
	double colorFieldNormalMagnitude;
	
	for (neighborIdx = 0; neighborIdx < PARTICLE_NUM; neighborIdx++)
	{
		if (neighborIdx == particleIdx) continue;

		VEC diffPosition;
		diffPosition.s0 = particles[particleIdx].position.s0 - particles[neighborIdx].position.s0;
		diffPosition.s1 = particles[particleIdx].position.s1 - particles[neighborIdx].position.s1;
		diffPosition.s2 = particles[particleIdx].position.s2 - particles[neighborIdx].position.s2;
		double radiusSquared = dotProduct(diffPosition, diffPosition);

		if (radiusSquared <= h*h) {

			VEC poly6Gradient, spikyGradient;

			poly6Gradient = Wpoly6Gradient(diffPosition, radiusSquared);
			spikyGradient = WspikyGradient(diffPosition, radiusSquared);

			f_pressure.s0 += (particles[particleIdx].pressure / pow(particles[particleIdx].density, 2) + particles[neighborIdx].pressure / pow(particles[neighborIdx].density, 2)) * spikyGradient.s0;
			f_pressure.s1 += (particles[particleIdx].pressure / pow(particles[particleIdx].density, 2) + particles[neighborIdx].pressure / pow(particles[neighborIdx].density, 2)) * spikyGradient.s1;
			f_pressure.s2 += (particles[particleIdx].pressure / pow(particles[particleIdx].density, 2) + particles[neighborIdx].pressure / pow(particles[neighborIdx].density, 2)) * spikyGradient.s2;
			
			f_viscosity.s0 += (particles[neighborIdx].velocity.s0 - particles[particleIdx].velocity.s0) *WviscosityLaplacian(radiusSquared) / particles[neighborIdx].density;
			f_viscosity.s1 += (particles[neighborIdx].velocity.s1 - particles[particleIdx].velocity.s1) *WviscosityLaplacian(radiusSquared) / particles[neighborIdx].density;
			f_viscosity.s2 += (particles[neighborIdx].velocity.s2 - particles[particleIdx].velocity.s2) *WviscosityLaplacian(radiusSquared) / particles[neighborIdx].density;
		
			colorFieldNormal.s0 += poly6Gradient.s0 / particles[neighborIdx].density;
			colorFieldNormal.s1 += poly6Gradient.s1 / particles[neighborIdx].density;
			colorFieldNormal.s2 += poly6Gradient.s2 / particles[neighborIdx].density;
			colorFieldLaplacian += Wpoly6Laplacian(radiusSquared) / particles[neighborIdx].density;
		}
			

	}

	f_pressure.s0 *= -PARTICLE_MASS * particles[particleIdx].density;
	f_pressure.s1 *= -PARTICLE_MASS * particles[particleIdx].density;
	f_pressure.s2 *= -PARTICLE_MASS * particles[particleIdx].density;
	f_viscosity.s0 *= VISCOSITY * PARTICLE_MASS;
	f_viscosity.s1 *= VISCOSITY * PARTICLE_MASS;
	f_viscosity.s2 *= VISCOSITY * PARTICLE_MASS;
	colorFieldNormal.s0 *= PARTICLE_MASS;
	colorFieldNormal.s1 *= PARTICLE_MASS;
	colorFieldNormal.s2 *= PARTICLE_MASS;
	colorFieldLaplacian *= PARTICLE_MASS;
	
	colorFieldNormalMagnitude = sqrt(colorFieldNormal.s0 * colorFieldNormal.s0 + colorFieldNormal.s1 * colorFieldNormal.s1 + colorFieldNormal.s2 * colorFieldNormal.s2);
	if (colorFieldNormalMagnitude > SURFACE_THRESHOLD) {
		particles[particleIdx].isSurface = 1;
		f_surface.s0 = -SURFACE_TENSION * colorFieldNormal.s0 / colorFieldNormalMagnitude * colorFieldLaplacian;
		f_surface.s1 = -SURFACE_TENSION * colorFieldNormal.s1 / colorFieldNormalMagnitude * colorFieldLaplacian;
		f_surface.s2 = -SURFACE_TENSION * colorFieldNormal.s2 / colorFieldNormalMagnitude * colorFieldLaplacian;
	} else {
		particles[particleIdx].isSurface = 0;
	}

	particles[particleIdx].acceleration.s0 = (f_pressure.s0 + f_viscosity.s0 + f_surface.s0 + f_gravity.s0) / particles[particleIdx].density;
	particles[particleIdx].acceleration.s1 = (f_pressure.s1 + f_viscosity.s1 + f_surface.s1 + f_gravity.s1) / particles[particleIdx].density;
	particles[particleIdx].acceleration.s2 = (f_pressure.s2 + f_viscosity.s2 + f_surface.s2 + f_gravity.s2) / particles[particleIdx].density;

	for (i = 0; i < WALL_NUM; i++) {
		VEC diff;
		diff.s0 = walls[i].point.s0 - particles[particleIdx].position.s0;
		diff.s1 = walls[i].point.s1 - particles[particleIdx].position.s1;
		diff.s2 = walls[i].point.s2 - particles[particleIdx].position.s2;

		double d = dotProduct(diff, walls[i].normal) + 0.01;
		if (d > 0.0) {
			particles[particleIdx].acceleration.s0 += WALL_K * walls[i].normal.s0 * d;
			particles[particleIdx].acceleration.s1 += WALL_K * walls[i].normal.s1 * d;
			particles[particleIdx].acceleration.s2 += WALL_K * walls[i].normal.s2 * d;
			particles[particleIdx].acceleration.s0 += WALL_DAMPING * walls[i].normal.s0 * dotProduct(particles[particleIdx].velocity, walls[i].normal);
			particles[particleIdx].acceleration.s1 += WALL_DAMPING * walls[i].normal.s1 * dotProduct(particles[particleIdx].velocity, walls[i].normal);
			particles[particleIdx].acceleration.s2 += WALL_DAMPING * walls[i].normal.s2 * dotProduct(particles[particleIdx].velocity, walls[i].normal);
		}
	}

	barrier(CLK_GLOBAL_MEM_FENCE);
	// STEP 3: INTEGRAL

	VEC newPosition;
	newPosition.s0 = particles[particleIdx].position.s0 + particles[particleIdx].velocity.s0*dt + particles[particleIdx].acceleration.s0*dt*dt;
	newPosition.s1 = particles[particleIdx].position.s1 + particles[particleIdx].velocity.s1*dt + particles[particleIdx].acceleration.s1*dt*dt;
	newPosition.s2 = particles[particleIdx].position.s2 + particles[particleIdx].velocity.s2*dt + particles[particleIdx].acceleration.s2*dt*dt;

	VEC newVelocity;
	newVelocity.s0 = (newPosition.s0 - particles[particleIdx].position.s0) / dt;
	newVelocity.s1 = (newPosition.s1 - particles[particleIdx].position.s1) / dt;
	newVelocity.s2 = (newPosition.s2 - particles[particleIdx].position.s2) / dt;

	particles[particleIdx].position = newPosition;
	particles[particleIdx].velocity = newVelocity;

}