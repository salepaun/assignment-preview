 #include "StableFluidsSim.h"
#include <Eigen/LU>

void StableFluidsSim::diffuseD(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
  assert((*x0 == *x0).all());
  
  scalar a = diff * dt * N * N;
  *x = *x0;
  
  for (int k = 0; k < 20; k++)
  {
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // Your code goes here!
        // Do diffuse for ([1, N], [1, N])
      }
    }
  }
}

void StableFluidsSim::diffuseU(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
  assert((*x0 == *x0).all());
  
  scalar a = diff * dt * N * N;
  *x = *x0;
  
  for (int k = 0; k < 20; k++)
  {
    for (int i = 1; i <= N; i++)
    {
      for (int j = 0; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // Your code goes here! 
        // Do diffuse for ([1, N], [0, N]), note the case when (j == 0) or (j == N) need special treatment
        if (j == 0)
          ;
        else if (j == N)
          ;
        else
          ;
      }
    }
  }
}

void StableFluidsSim::diffuseV(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
  assert((*x0 == *x0).all());
  
  scalar a = diff * dt * N * N;
  *x = *x0;
  
  for (int k = 0; k < 20; k++)
  {
    for (int i = 0; i <= N; i++)
    {
      for (int j = 1; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // Your code goes here!
        if (i == 0)
          ; 
        else if (i == N)
          ; 
        else
          ;
      }
    }
  }
}

void StableFluidsSim::advectD(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt)
{
  assert((*x0 == *x0).all());
  assert((*u == *u).all());
  assert((*v == *v).all());
  
  // Your code goes here!
  // Advect for ([1, N], [1, N])
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
    }
  }
}

scalar StableFluidsSim::interpolateD(ArrayXs * d, scalar i, scalar j)
{
    // Your code goes here!
    // Note the indices should be CLAMP-ed to [0, m_N], since we have to use (i + 1) and (j + 1)
}

scalar StableFluidsSim::interpolateU(ArrayXs * u, scalar i, scalar j)
{
    // Your code goes here! 
    // Note the i index should be CLAMP-ed to [0, m_N], while j index should be CLAMP-ed to [0, m_N-1], since we have to use (i + 1) and (j + 1)
}

scalar StableFluidsSim::interpolateV(ArrayXs * v, scalar i, scalar j)
{
    // Your code goes here!
}

void StableFluidsSim::advectU(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt)
{
  assert((*x0 == *x0).all());
  assert((*u == *u).all());
  assert((*v == *v).all());
  
  for (int i = 1; i <= N; i++)
  {
    for (int j = 0; j <= N; j++)
    {
      // Your code goes here!
      
      // Add the origin of U grid to the coordinate before sampling, for example, sample at (i + 0, j + 0.5) when you need backtracing the old velocity at (i, j)
      
      // Now you have the backward-traced velocity, minus it from the current position (i + 0, j + 0.5), then sample the velocity again.
    }
  }
}

void StableFluidsSim::advectV(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt)
{
    assert((*x0 == *x0).all());
    assert((*u == *u).all());
    assert((*v == *v).all());
  
  for (int i = 0; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      // Your code goes here!
    }
  }
}

void StableFluidsSim::project(int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0)
{
  if (VERBOSE) std::cout << "u0: " << std::endl << *u0 << std::endl << std::endl;
  if (VERBOSE) std::cout << "v0: " << std::endl << *v0 << std::endl << std::endl;

  ArrayXs div(N + 2, N + 2);
  ArrayXs p(N + 2, N + 2);
  div.setZero();
  p.setZero();
  scalar h = 1.0 / N;
  
  // Your code goes here!
  
  // set solid boundary conditions, 0 the most top and bottom row / left and right column of u0, v0
  
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      // compute divergence of the velocity field, note the divergence field is available from ([1, N], [1, N])
    }
  }
  
  for (int k = 0; k < 20; k++)
  {
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // solve for pressure inside the region ([1, N], [1, N])
      }
    }
  }
  
  (*u) = (*u0);
  (*v) = (*v0);
  
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j < N; j++)
    {
      // apply pressure to correct velocities ([1, N], [1, N)) for u, ([1, N), [1, N]) for v
    }
  }
}

