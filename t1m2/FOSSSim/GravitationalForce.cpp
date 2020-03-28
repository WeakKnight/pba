#include "GravitationalForce.h"

void GravitationalForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );

  // Your code goes here!
  
  VectorXs a(2);
  a.setZero();
  
  VectorXs b(2);
  b.setZero();
  
  a[0] = x[m_particles.first * 2];
  a[1] = x[m_particles.first * 2 + 1];
  
  b[0] = x[m_particles.second * 2]; 
  b[1] = x[m_particles.second * 2 + 1];
  
  scalar a_m = m[m_particles.first * 2];
  scalar b_m = m[m_particles.second * 2];
  
  scalar len = (a - b).norm();
  
  E += (m_G * a_m * b_m / len);
}

void GravitationalForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );

  VectorXs a(2);
  a.setZero();
  
  VectorXs b(2);
  b.setZero();
  
  // Your code goes here!
  a[0] = x[m_particles.first * 2];
  a[1] = x[m_particles.first * 2 + 1];
  
  b[0] = x[m_particles.second * 2]; 
  b[1] = x[m_particles.second * 2 + 1];
  
  scalar a_m = m[m_particles.first * 2];
  scalar b_m = m[m_particles.second * 2];
  
  scalar lenSquare = (a - b).squaredNorm();
  
  VectorXs unit(2);
  unit.setZero();
  unit = (a - b).normalized();
  
  VectorXs a_gradE(2);
  a_gradE.setZero();
  
  VectorXs b_gradE(2);
  b_gradE.setZero();
  
  a_gradE = (m_G * a_m * b_m / lenSquare) * unit;
  b_gradE = (-1.0 * m_G * a_m * b_m / lenSquare) * unit;
  
  gradE[m_particles.first * 2] += a_gradE[0];
  gradE[m_particles.first * 2 + 1] += a_gradE[1];
  
  gradE[m_particles.second * 2] += b_gradE[0];
  gradE[m_particles.second * 2 + 1] += b_gradE[1];
}
