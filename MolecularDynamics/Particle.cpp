#include "Particle.h"

Particle::Particle(double m, double x, double y, double vx, double vy)
{
    this->m = m;
    this->x = x;
    this->y = y;
    this->vx = vx;
    this->vy = vy;
    this->init = true;
    this->x_minus_dt = 0.0;
    this->y_minus_dt = 0.0;
    this->x_plus_dt = 0.0;
    this->y_plus_dt = 0.0;
}

double Particle::kinetic_energy()
{
    return 0.5 * m * (vx * vx + vy * vy);
}

void Particle::record_path()
{
    x_path.push_back(x);
    y_path.push_back(y);
    vx_path.push_back(vx);
    vy_path.push_back(vy);
}

void Particle::record_minus_dt()
{
    x_minus_dt = x;
    y_minus_dt = y;
}

void Particle::update_velocity(double dt)
{
    vx = (x_plus_dt - x_minus_dt) / (2 * dt);
    vy = (y_plus_dt - y_minus_dt) / (2 * dt);
}


