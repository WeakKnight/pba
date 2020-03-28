#include "SpringForce.h"

void SpringForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Your code goes here!
  VectorXs a(2);
  a.setZero();

  VectorXs b(2);
  b.setZero();

  a[0] = x[m_endpoints.first * 2];
  a[1] = x[m_endpoints.first * 2 + 1];

  b[0] = x[m_endpoints.second * 2]; 
  b[1] = x[m_endpoints.second * 2 + 1];

  scalar singleEnergy = 0.5 * m_k * ((a - b).norm() - m_l0) * ((a - b).norm() - m_l0);
  E += singleEnergy;
}

void SpringForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Your code goes here!
    
  VectorXs a(2);
  a.setZero();

  VectorXs b(2);
  b.setZero();

  a[0] = x[m_endpoints.first * 2];
  a[1] = x[m_endpoints.first * 2 + 1];

  b[0] = x[m_endpoints.second * 2]; 
  b[1] = x[m_endpoints.second * 2 + 1];

  // from b to a
  VectorXs unit(2);
  unit.setZero();
  unit = (a - b).normalized();

  scalar len = (a - b).norm();

  VectorXs a_gradE(2);
  a_gradE.setZero();

  VectorXs b_gradE(2);
  b_gradE.setZero();

  a_gradE = m_k * (len - m_l0) * unit;
  b_gradE = -1.0 * m_k * (len - m_l0) * unit;

  gradE[m_endpoints.first * 2] += a_gradE[0];
  gradE[m_endpoints.first * 2 + 1] += a_gradE[1];

  gradE[m_endpoints.second * 2] += b_gradE[0];
  gradE[m_endpoints.second * 2 + 1] += b_gradE[1];

  // damping
  VectorXs a_dampE(2);
  a_dampE.setZero();

  VectorXs b_dampE(2);
  b_dampE.setZero();
  
  VectorXs a_v(2);
  a_v.setZero();

  VectorXs b_v(2);
  b_v.setZero();

  a_v[0] = v[m_endpoints.first * 2];
  a_v[1] = v[m_endpoints.first * 2 + 1];

  b_v[0] = v[m_endpoints.second * 2]; 
  b_v[1] = v[m_endpoints.second * 2 + 1];

  VectorXs aMinusB = a_v - b_v;
  scalar dotProduct = unit[0] * aMinusB[0] + unit[1] * aMinusB[1];
  
  a_dampE = m_b * dotProduct * unit;
  b_dampE = -1.0 * m_b * dotProduct * unit;

  gradE[m_endpoints.first * 2] += a_dampE[0];
  gradE[m_endpoints.first * 2 + 1] += a_dampE[1];

  gradE[m_endpoints.second * 2] += b_dampE[0];
  gradE[m_endpoints.second * 2 + 1] += b_dampE[1];
}
