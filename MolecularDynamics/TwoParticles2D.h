#pragma once
#include "Particle.h"
#include <vector>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>


using namespace std;

class TwoParticles2D
{
public:
	Particle p1;
	Particle p2;
	double tmax;
	double dt;
	double coeff;
	bool init;
	bool ran;
	bool periodic_bound;
	double x_max;
	double x_min;
	double y_max;
	double y_min;
	double epsilon;
	double sigma;

	void run();
	void csv_output(string filename);
	void xyz_output(string filename);
	void console_output();
	TwoParticles2D(Particle p1, Particle p2, double tmax = 1, 
		double coeff = 0.005, bool periodic_bound = false, double epsilon = 1, double sigma = 1);
	~TwoParticles2D();
private:
	vector<double> _time;
	vector<double> _kinetic;
	vector<double> _potential;
	Particle* _virtual_particle;
	double _t;
	double _r6;
	double _r12;
	double _r2;
	void create_virtual_particles();
	void destroy_virtual_particles();
	double force_x(Particle& p1, Particle& p2);
	double force_y(Particle& p1, Particle& p2);
	double potential();
	double square_distance(Particle& p1, Particle& p2);
	void update_r2_r6_r12();
	void measure(double U);
	void initiate_verlet(double fx, double fy);
	void verlet_update_pos(double fx, double fy);
	void keep_in_bound(Particle& p);
	void update_velocity(Particle& p);
	string pos_string(Particle& p, size_t i);
	void get_U_fx_fy(double& U, double& fx, double& fy);
};


