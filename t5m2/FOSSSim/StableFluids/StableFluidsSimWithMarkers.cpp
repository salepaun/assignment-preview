#include "StableFluidsSimWithMarkers.h"
#include <Eigen/LU>

//#define VERBOSE (1)
#define SOLVER  (0)     
#define VOLUME_COMPENSATION_ETA (0.01)

StableFluidsSimWithMarkers::StableFluidsSimWithMarkers( const int& rows, const int& cols, const ArrayXb & has_solid, const scalar& diff, const scalar& visc)
: StableFluidsSim(rows, cols, diff, visc)
, m_has_fluid(m_N + 2, m_N + 2)
, m_has_solid(has_solid)
{
    assert(rows==cols);
    
    clear();
}

StableFluidsSimWithMarkers::~StableFluidsSimWithMarkers()
{

}

//////////////////////////////////////////////////////////////////////////////////////

bool StableFluidsSimWithMarkers::in_solid(const Vector2s & pos)
{
    int i = std::max(1, std::min(m_N, (int)(pos.x() * m_N + 0.5)));
    int j = std::max(1, std::min(m_N, (int)(pos.y() * m_N + 0.5)));
    return m_has_solid(i, j);
}

void StableFluidsSimWithMarkers::detect_fluid(int N, ArrayXb & fluids, const std::vector<Vector2s> & particles)
{
    fluids.setZero();
    for (size_t p = 0; p < particles.size(); p++)
    {
        int i = std::max(1, std::min(N, (int)(particles[p].x() * N + 0.5)));
        int j = std::max(1, std::min(N, (int)(particles[p].y() * N + 0.5)));
        fluids(i, j) = true;
    }
}

void StableFluidsSimWithMarkers::solveByJacobi(ArrayXs & p, ArrayXs & div)
{
  // TODO: Implement the Jacobi Solver with Ping-pong buffers
  int N = m_N;
  scalar h = 1.0 / N;
  
  ArrayXs p_copy = p;
  
  ArrayXs* ps[2] = {&p, &p_copy};
  
  for (int k = 0; k < 500; k++)
  {
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= N; j++)
      {
        // you may use (*ps[!(k & 1)]) to refer to the next state and (*ps[(k & 1)]) for current state
      }
    }
  }
}

void StableFluidsSimWithMarkers::solveByGaussSeidel(ArrayXs & p, ArrayXs & div)
{
    int N = m_N;
    scalar h = 1.0 / N;

    for (int k = 0; k < 500; k++)
    {
        for (int i = 1; i <= N; i++)
        {
            for (int j = 1; j <= N; j++)
            {
                if (m_has_solid(i, j) || !m_has_fluid(i, j))
                {
                    p(i, j) = 0;
                    continue;
                }
                
                scalar numer = div(i, j) * (h * h);
                int denom = -4;
                
                if (m_has_solid(i - 1, j))      denom++;
                else if (m_has_fluid(i - 1, j)) numer -= p(i - 1, j);
                if (m_has_solid(i + 1, j))      denom++;
                else if (m_has_fluid(i + 1, j)) numer -= p(i + 1, j);
                if (m_has_solid(i, j - 1))      denom++;
                else if (m_has_fluid(i, j - 1)) numer -= p(i, j - 1);
                if (m_has_solid(i, j + 1))      denom++;
                else if (m_has_fluid(i, j + 1)) numer -= p(i, j + 1);
                
                assert(denom != 0);
                p(i, j) = numer / denom;
            }
        }
    }
}

void StableFluidsSimWithMarkers::solveByEigenLU(ArrayXs & p, ArrayXs & div)
{
    int N = m_N;
    scalar h = 1.0 / N;
    
    MatrixXs A = MatrixXs::Zero(N * N, N * N);
    VectorXs b = VectorXs::Zero(N * N);
    
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (m_has_solid(i + 1, j + 1))
            {
                // solid cell
                A(i * N + j, i * N + j) = 1;
                b(i * N + j) = 0;
                
            } else if (!m_has_fluid(i + 1, j + 1))
            {
                // air cell
                A(i * N + j, i * N + j) = 1;
                b(i * N + j) = 0;
                
            } else
            {
                // fluid cell
                b(i * N + j) = div(i + 1, j + 1) / (N * N);
                A(i * N + j, i * N + j) = -4;
                
                if      (m_has_solid(i + 2, j + 1)) A(i * N + j, i * N + j)++;
                else if (m_has_fluid(i + 2, j + 1)) A(i * N + j, (i + 1) * N + j) = 1;
                if      (m_has_solid(i + 0, j + 1)) A(i * N + j, i * N + j)++;
                else if (m_has_fluid(i + 0, j + 1)) A(i * N + j, (i - 1) * N + j) = 1;
                if      (m_has_solid(i + 1, j + 2)) A(i * N + j, i * N + j)++;
                else if (m_has_fluid(i + 1, j + 2)) A(i * N + j, i * N + (j + 1)) = 1;
                if      (m_has_solid(i + 1, j + 0)) A(i * N + j, i * N + j)++;
                else if (m_has_fluid(i + 1, j + 0)) A(i * N + j, i * N + (j - 1)) = 1;
                
            }
        }
    }
    
    VectorXs pv = A.fullPivLu().solve(b);
    
    for (int i = 1; i <= N; i++)
        for (int j = 1; j <= N; j++)
            p(i, j) = pv((i - 1) * N + (j - 1));
}

void StableFluidsSimWithMarkers::solveByCG(ArrayXs & p, ArrayXs & div)
{
    
}

void StableFluidsSimWithMarkers::solveByPCG(ArrayXs & p, ArrayXs & div)
{
    
}

void StableFluidsSimWithMarkers::project(int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0)
{
  //////////////////////////////////
  // Your code goes here!
  
  ArrayXs div(N + 2, N + 2);
  ArrayXs p(N + 2, N + 2);
  div.setZero();
  p.setZero();
  scalar h = 1.0 / N;
  
  // set velocities on solid boundaries to zero
  for (int i = 0; i <= N + 1; i++)
  {
    for (int j = 0; j <= N + 1; j++)
    {
      
    }
  }
    // compute the discrete divergence of the velocity field, with modifications to compensate for volume change
    // Hint: use fluid_cell_count() as current and m_original_fluid_cell_count as the initial number of liquid cells
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
    }
  }
  
    // solve for pressure
  solveByJacobi(p, div);
  
    // apply the pressure gradient to the velocity field
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j < N; j++)
    {
    }
  }
  
    // set solid boundary conditions again to overwrite the wrong results coming from applying pressures above
  for (int i = 0; i <= N + 1; i++)
  {
    for (int j = 0; j <= N + 1; j++)
    {
    }
  }
    // propagate liquid velocities to nearby air cells, up to 5 layers
    // step 1: test if both adjacent cells are neither solids nor fluids; if so, set velocity to a sentinel value (SCALAR_INFINITY)
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j < N; j++)
    {
    }
  }
  
    // step 2: if the velocity of current cell is SCALAR_INFINITY, propagate all the neighbor cells that don't have a SCALAR_INFINITY to current cell and average them by dividing the number of propagated cells.
  for (int k = 0; k < 5; k++)
  {
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j < N; j++)
      {
        int ip = CLAMP(i + 1, 1, N);
        int im = CLAMP(i - 1, 1, N);
        int jp = CLAMP(j + 1, 1, N - 1);
        int jm = CLAMP(j - 1, 1, N - 1);
        
        if ((*u)(i, j) == SCALAR_INFINITY)
        {
        }
        if ((*v)(j, i) == SCALAR_INFINITY)
        {
        }
      }
    }
  }
  
    // step 3: replace all the remaining velocities = SCALAR_INFINITY by 0.
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j < N; j++)
    {
    }
  }
}

void StableFluidsSimWithMarkers::dens_step(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar diff, scalar dt)
{
    ArrayXs * outu = u;
    ArrayXs * outv = v;
    
    add_source(N, x, x0, dt); 
    
    SWAP(x0, x); 
    diffuseD(N, x, x0, diff, dt); 
    
    SWAP(x0, x); 
    advectD(N, x, x0, u, v, dt);
    
    if (outu != u)
        *outu = *u;    
    if (outv != v)
        *outv = *v;
    
}

void StableFluidsSimWithMarkers::vel_step(int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0, scalar visc, scalar dt)
{
    ArrayXs * outu = u;
    ArrayXs * outv = v;
    
    add_source(N, u, u0, dt); 
    add_source(N, v, v0, dt); 
    
    // gravity
    SWAP(u0, u); 
    SWAP(v0, v); 
    scalar g = 10.0;
    *u = *u0;
    *v = *v0 + g * dt;
    
    SWAP(u0, u);
    SWAP(v0, v);
    project(N, u, v, u0, v0);
  
    m_uAfterDiffusion.setZero();
    m_vAfterDiffusion.setZero();
    m_uAfterDiffusion = *u;
    m_vAfterDiffusion = *v;
    // diffusion
//    SWAP(u0, u); 
//    SWAP(v0, v); 
//    diffuseU(N, u, u0, visc, dt); 
//    diffuseV(N, v, v0, visc, dt); 
    
//    SWAP(u0, u);
//    SWAP(v0, v);
//    project(N, u, v, u0, v0);
    
    // advection
    SWAP(u0, u); 
    SWAP(v0, v); 
    advectU(N, u, u0, u0, v0, dt); 
    advectV(N, v, v0, u0, v0, dt); 
  
    m_uAfterAdvect.setZero();
    m_vAfterAdvect.setZero();
    m_uAfterAdvect = *u;
    m_vAfterAdvect = *v;
  
    SWAP(u0, u);
    SWAP(v0, v);
    project(N, u, v, u0, v0);
    
    if (outu != u)
        *outu = *u;    
    if (outv != v)
        *outv = *v;
    
}

//////////////////////////////////////////////////////////////////////////////////////

void StableFluidsSimWithMarkers::advect_particles(int N, ArrayXs & u, ArrayXs & v, scalar dt)
{
    std::vector<Vector2s> & particles = m_particles;
  
  //////////////////////////////////
  // Your code goes here!
  // Hint: Use interpolateV and interpolateU to get the new position of a particle
  // then CLAMP it to the range [0.5 / N, (N + 0.4999) / N]
    const double c_1 = 0.5;
    const double c_2 = 0.4999;
    
    
    detect_fluid(N, m_has_fluid, particles);
}

//////////////////////////////////////////////////////////////////////////////////////

void StableFluidsSimWithMarkers::stepSystem( const scalar& dt)
{
    if (VERBOSE) std::cout << "step" << std::endl;
    if (VERBOSE) std::cout << "fluid cell count = " << fluid_cell_count() << std::endl;

    ArrayXs new_d(m_N + 2, m_N + 2);
    ArrayXs new_u(m_N + 2, m_N + 1);
    ArrayXs new_v(m_N + 1, m_N + 2);
    
    new_d.setZero();
    new_u.setZero();
    new_v.setZero();
    
    vel_step(m_N, &new_u, &new_v, &m_u, &m_v, m_visc, dt);
    advect_particles(m_N, new_u, new_v, dt);
    
    dens_step(m_N, &new_d, &m_d, &new_u, &new_v, m_diff, dt);
    
    m_d = new_d;
    m_u = new_u;
    m_v = new_v;
    
}

void StableFluidsSimWithMarkers::clear()
{
    StableFluidsSim::clear();
    
    m_particles.clear();
    for (int i = m_N / 2; i <= m_N; i++)
    {
        for (int j = 1; j <= m_N; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                if (!m_has_solid(i, j))
                {
                    m_particles.push_back(Vector2s(i + (double)rand() / RAND_MAX - 0.5, j + (double)rand() / RAND_MAX - 0.5) / m_N);
                }
            }
        }
    }
    
    detect_fluid(m_N, m_has_fluid, m_particles);
    m_original_fluid_cell_count = fluid_cell_count();
}


