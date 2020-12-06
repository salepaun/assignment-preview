 #include "StableFluidsSim.h"
#include <Eigen/LU>


#define LERP(a,b,x) (1-x)*a + x*b;




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
        (*x)(i,j) = (*x0)(i,j) + a * ((*x)(i-1,j)+(*x)(i+1,j)+(*x)(i,j-1)+(*x)(i,j+1)) / (1+4*a);
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
          (*x)(i,j) = (*x0)(i,j) + a * ((*x)(i-1,j)+(*x)(i+1,j)+(*x)(i,j+1)) / (1+3*a);
        else if (j == N)
          (*x)(i,j) = (*x0)(i,j) + a * ((*x)(i-1,j)+(*x)(i+1,j)+(*x)(i,j-1)) / (1+3*a);
        else
          (*x)(i,j) = (*x0)(i,j) + a * ((*x)(i-1,j)+(*x)(i+1,j)+(*x)(i,j-1)+(*x)(i,j+1)) / (1+4*a);
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
          (*x)(i,j) = (*x0)(i,j) + a * ((*x)(i+1,j)+(*x)(i,j-1)+(*x)(i,j+1)) / (1+3*a);
        else if (i == N)
          (*x)(i,j) = (*x0)(i,j) + a * ((*x)(i-1,j)+(*x)(i,j-1)+(*x)(i,j+1)) / (1+3*a);
        else
          (*x)(i,j) = (*x0)(i,j) + a * ((*x)(i-1,j)+(*x)(i+1,j)+(*x)(i,j-1)+(*x)(i,j+1)) / (1+4*a);
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
  double dt0, px, py;

  dt0 = dt*N;

  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      px = i - dt0 * (*u)(i,j);
      py = j - dt0 * (*v)(i,j);
      (*x)(i,j) = interpolateD(x0, px, py);
    }
  }
}

scalar StableFluidsSim::interpolateD(ArrayXs * d, scalar i, scalar j)
{
  // Your code goes here!
  // Note the indices should be CLAMP-ed to [0, m_N], since we have to use (i + 1) and (j + 1)
  int i0, j0, i1, j1;

  /*
  double s0, s1, t0, t1;
  i0 = CLAMP(int(i), 0, m_N);
  i1 = int(i) + 1;
  j0 = CLAMP(int(j), 0, m_N);
  j1 = int(j) + 1;

  s1 = i-i0; s0 = 1-s1; t1 = j-j0; t0=1-t1;
  return s0 * (t0*(*d)(i0,j0) + t1*(*d)(i0,j1)) + s1 * (t0*(*d)(i1,j0) + t1*(*d)(i1,j1));
  */

  double s, t, d0, d1;

  i0 = CLAMP(int(i), 0, m_N);
  i1 = i0 + 1;
  j0 = CLAMP(int(j), 0, m_N);
  j1 = j0 + 1;
  s = CLAMP(i-i0, 0, 1);
  t = CLAMP(j-j0, 0, 1);
  d0 = LERP((*d)(i0,j0), (*d)(i1,j0), s);
  d1 = LERP((*d)(i0,j1), (*d)(i1,j1), s);
      
  return LERP(d0, d1, t);
}

scalar StableFluidsSim::interpolateU(ArrayXs * u, scalar i, scalar j)
{
  // Your code goes here! 
  // Note the i index should be CLAMP-ed to [0, m_N], while j index should be CLAMP-ed to [0, m_N-1], since we have to use (i + 1) and (j + 1)
  int i0, j0, i1, j1;
  double s, t, u0, u1;

  i0 = CLAMP(int(i), 0, m_N);
  i1 = i0 + 1;
  j0 = CLAMP(int(j-0.5), 0, m_N-1);
  j1 = j0 + 1;
  s = CLAMP(i-i0, 0, 1);
  t = CLAMP(j-j0-0.5, 0, 1);
  u0 = LERP((*u)(i0,j0), (*u)(i1,j0), s);
  u1 = LERP((*u)(i0,j1), (*u)(i1,j1), s);

  return LERP(u0, u1, t);
}

scalar StableFluidsSim::interpolateV(ArrayXs * v, scalar i, scalar j)
{
  // Your code goes here!
  int i0, j0, i1, j1;
  double s, t, v0, v1;

  i0 = CLAMP(int(i-0.5), 0, m_N-1);
  i1 = i0 + 1;
  j0 = CLAMP(int(j), 0, m_N);
  j1 = j0 + 1;
  s = CLAMP(i-i0-0.5, 0, 1);
  t = CLAMP(j-j0, 0, 1);
  v0 = LERP((*v)(i0,j0), (*v)(i1,j0), s);
  v1 = LERP((*v)(i0,j1), (*v)(i1,j1), s);

  return LERP(v0, v1, t);
}

void StableFluidsSim::advectU(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt)
{
  assert((*x0 == *x0).all());
  assert((*u == *u).all());
  assert((*v == *v).all());
  
  double ux, uy, dt0;
  dt0 = dt*N;

  for (int i = 1; i <= N; i++)
  {
    for (int j = 0; j <= N; j++)
    {
      // Your code goes here!
      
      // Add the origin of U grid to the coordinate before sampling, for example, sample at (i + 0, j + 0.5) when you need backtracing the old velocity at (i, j)
      ux = i - dt0 * (*u)(i,j);
      uy = j - dt0 * (*v)(i,j);
      (*x)(i,j) = interpolateU(x0, ux, uy);
      
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

