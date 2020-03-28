#include "SymplecticEuler.h"
#include <stdio.h>

bool SymplecticEuler::stepScene( TwoDScene& scene, scalar dt )
{
  // Your code goes here!
  VectorXs& X = scene.getX();
  VectorXs& V = scene.getV();
  const VectorXs& M = scene.getM();
  
  VectorXs F(X.size());
  F.setZero();

  scene.accumulateGradU(F);
  F = -F;
  
  for(int i = 0; i < (X.size() / 2); i++)
  {
      if(scene.isFixed(i))
      {
          continue;
      }
    
      VectorXs vThis(2);
      vThis.setZero();
    
      vThis[0] = V[i * 2]; 
      vThis[1] = V[i * 2 + 1];
      
      VectorXs xThis(2);
      xThis.setZero();
    
      xThis[0] = X[i * 2];
      xThis[1] = X[i * 2 + 1];
    
      VectorXs mInverse(2);
      mInverse.setZero();
      mInverse[0] = 1.0 / M[i * 2];
      mInverse[1] = 1.0 / M[i * 2 + 1];
      
      VectorXs fThis(2);
      fThis.setZero();
      fThis[0] = F[i * 2]; 
      fThis[1] = F[i * 2 + 1];
      
      VectorXs newV(2);
      newV.setZero();
    
      VectorXs newX(2);
      newX.setZero();
    
      newV[0] = vThis[0] + dt * mInverse[0] * fThis[0];
      newV[1] = vThis[1] + dt * mInverse[1] * fThis[1];
      newX[0] = xThis[0] + dt * newV[0];
      newX[1] = xThis[1] + dt * newV[1];
    
      X[i * 2] = newX[0];
      X[i * 2 + 1] = newX[1];
    
      V[i * 2] = newV[0];
      V[i * 2 + 1] = newV[1];
  }
  
  return true;
}






