#include "TwoParticles2D.h"


void TwoParticles2D::run()
{
	if (init == false) {
		throw std::runtime_error("initialization required");
	}

	double U = 0.0;
	double fx = 0.0;
	double fy = 0.0;

	get_U_fx_fy(U, fx, fy);
	// must not change within an iteration
	
	p1.x_minus_dt = p1.x;
	p1.y_minus_dt = p1.y;
	p2.x_minus_dt = p2.x;
	p2.y_minus_dt = p2.y;
	p1.record_path();
	p2.record_path();
	measure(U);
	initiate_verlet(fx, fy);
	keep_in_bound(p1);
	keep_in_bound(p2);
	
	
	get_U_fx_fy(U, fx, fy);


	while (_t < tmax)
	{
		p1.record_path();
		p2.record_path();
		measure(U);
		verlet_update_pos(fx, fy);
		keep_in_bound(p1);
		keep_in_bound(p2);

		get_U_fx_fy(U, fx, fy);
		_t += dt;
	}
	ran = true;
}

void TwoParticles2D::console_output()
{
	cout << "Total Time Steps: " << _time.size() << endl;
	cout << setw(10) << "time: " << setw(15) << "p1 pos" << setw(15) << "p1 vel"
		<< setw(15) << "p2 pos" << setw(15) << "p2 vel" << setw(10) << "energy" << endl;
	for (size_t i = 0; i != _time.size(); i++) {
		string p1_pos = pos_string(p1, i);
		string p2_pos = pos_string(p2, i);
		cout << setw(10) << _time[i] << setw(15) << p1_pos << setw(15) << " "
			<< setw(15) << p2_pos << setw(15) << " " << setw(10) << _potential[i] + _kinetic[i] << endl;
	}

}

TwoParticles2D::TwoParticles2D(Particle p1, Particle p2, double tmax, double coeff, bool periodic_bound, double epsilon, double sigma)
{
	this->p1 = p1;
	this->p2 = p2;
	this->tmax = tmax;
	this->coeff = coeff;
	this->periodic_bound = periodic_bound;
	this->ran = false;
	this->epsilon = 1;
	this->sigma = 1;
	this->init = true;
	this->x_max = 10 * sigma;
	this->y_max = 10 * sigma;
	this->x_min = 0;
	this->y_min = 0;
	this->dt = tmax * coeff;
	this->_virtual_particle = nullptr;
	this->_r12 = 0;
	this->_r6 = 0;
	this->_r2 = square_distance(p1, p2);
	this->_t = 0;
}

TwoParticles2D::~TwoParticles2D()
{
	if (_virtual_particle != nullptr) {
		delete _virtual_particle;
	}
}

void TwoParticles2D::create_virtual_particles()
{
	if (periodic_bound == false) {
		_virtual_particle = nullptr;
		return;
	}
	double Lx = x_max - x_min;
	double Ly = y_max - y_min;
	if ((abs(p1.x - p2.x) > Lx / 2) || (abs(p1.y - p2.y) > Ly / 2)) {
		_virtual_particle = new Particle(p2);
		if (abs(p1.x - p2.x) > Lx / 2) {
			if (p1.x - p2.x < 0) {
				_virtual_particle->x = p2.x - Lx;
			}
			else {
				_virtual_particle->x = p2.x + Lx;
			}
		}
		if (abs(p1.y - p2.y) > Ly / 2) {
			if (p1.y - p2.y < 0) {
				_virtual_particle->y = p2.y - Ly;
			}
			else {
				_virtual_particle->y = p2.y + Ly;
			}
		}
	}
}

void TwoParticles2D::destroy_virtual_particles()
{
	if (_virtual_particle != nullptr) {
		delete _virtual_particle;
		_virtual_particle = nullptr;
	}
}

double TwoParticles2D::force_x(Particle& p1, Particle& p2)
{
	
	return (p2.x - p1.x) * (24 * pow(sigma, 6) * epsilon * (_r6 - 2 * pow(sigma, 6))) / (_r12 * _r2);
}

double TwoParticles2D::force_y(Particle& p1, Particle& p2)
{
	
	return (p2.y - p1.y) * (24 * pow(sigma, 6) * epsilon * (_r6 - 2 * pow(sigma, 6))) / (_r12 * _r2);
}

double TwoParticles2D::potential()
{
	return 4 * epsilon * (pow(sigma, 6) / _r6) * ((pow(sigma, 6) / _r6) - 1);
}

double TwoParticles2D::square_distance(Particle& p1, Particle& p2)
{
	return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

void TwoParticles2D::update_r2_r6_r12()
{
	_r2 = square_distance(p1, p2);
	_r6 = pow(_r2, 3);
	_r12 = _r6 * _r6;
}

void TwoParticles2D::measure(double U)
{
	_potential.push_back(U);
	_kinetic.push_back(p1.kinetic_energy() + p2.kinetic_energy());
	_time.push_back(_t);
}

void TwoParticles2D::initiate_verlet(double fx, double fy)
{
	p1.x += p1.vx * dt;
	p1.vx += fx * dt;
	p1.y += p1.vy * dt;
	p1.vy += fy * dt;
	p2.x += p2.vx * dt;
	p2.vx += -fx * dt;
	p2.y += p2.vy * dt;
	p2.vy += -fy * dt;
}

void TwoParticles2D::verlet_update_pos(double fx, double fy)
{
	double x1_last = p1.x;
	double x2_last = p2.x;
	double y1_last = p1.y;
	double y2_last = p2.y;

	p1.x = 2.0 * p1.x - p1.x_minus_dt + fx * dt * dt;
	p1.y = 2.0 * p1.y - p1.y_minus_dt + fy * dt * dt;
	p2.x = 2.0 * p2.x - p2.x_minus_dt - fx * dt * dt;
	p2.y = 2.0 * p2.y - p2.y_minus_dt - fy * dt * dt;

	p1.x_minus_dt = x1_last;
	p2.x_minus_dt = x2_last;
	p1.y_minus_dt = y1_last;
	p2.y_minus_dt = y2_last;
	// telescope one dt forward in time
	p1.x_plus_dt = 2.0 * p1.x - p1.x_minus_dt + fx * dt * dt;
	p1.y_plus_dt = 2.0 * p1.y - p1.y_minus_dt + fy * dt * dt;
	p2.x_plus_dt = 2.0 * p2.x - p2.x_minus_dt + fx * dt * dt;
	p2.y_plus_dt = 2.0 * p2.y - p2.y_minus_dt + fy * dt * dt;
	p1.update_velocity(dt);
	p2.update_velocity(dt);
}

void TwoParticles2D::keep_in_bound(Particle& p)
{
	double Lx = x_max - x_min;
	double Ly = y_max - y_min;
	if (p.x > x_max) {
		p.x = p.x - Lx;
	}
	if (p.x < x_min) {
		p.x = p.x + Lx;
	}
	if (p.y > y_max) {
		p.y = p.y - Ly;
	}
	if (p.y < y_min) {
		p.y = p.y + Ly;
	}
}

string TwoParticles2D::pos_string(Particle& p, size_t i)
{
	ostringstream oss1;
	oss1 << "(" << p.x_path[i] << "," << p.y_path[i] << ")";
	return oss1.str();
}

void TwoParticles2D::get_U_fx_fy(double& U, double& fx, double& fy)
{
	create_virtual_particles();
	if (_virtual_particle == nullptr) {
		update_r2_r6_r12();
		U = potential();
		fx = force_x(p1, p2);
		fy = force_y(p1, p2);
	}
	else {
		_r2 = square_distance(p1, *_virtual_particle);
		_r6 = pow(_r2, 3);
		_r12 = _r6 * _r6;
		U = potential();
		fx = force_x(p1, *_virtual_particle);
		fy = force_y(p1, *_virtual_particle);
	}
	destroy_virtual_particles();
}



