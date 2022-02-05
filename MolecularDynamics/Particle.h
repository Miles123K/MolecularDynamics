#pragma once

#include <vector>

using namespace std;

class Particle
{
public:
	double m;
	double x;
	double x_minus_dt;
	double x_plus_dt;
	double y;
	double y_minus_dt;
	double y_plus_dt;
	double vx;
	double vy;
	bool init;
	vector<double> x_path;
	vector<double> vx_path;
	vector<double> y_path;
	vector<double> vy_path;
	
	Particle(double m = 1, double x = 0, double y = 0, double vx = 0, double vy = 0);
	double kinetic_energy();
	void record_path();
	void record_minus_dt();
	void update_velocity(double dt);
};

