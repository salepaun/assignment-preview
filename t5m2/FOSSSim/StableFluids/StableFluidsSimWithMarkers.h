#ifndef __STABLE_FLUIDS_SIM_WITH_MARKERS_H__
#define __STABLE_FLUIDS_SIM_WITH_MARKERS_H__

#include "StableFluidsSim.h"

class StableFluidsSimWithMarkers : public StableFluidsSim
{
public:
    
    // Initialize the fluid simulation. Assumes rows == cols.
    //   rows: Number of rows in the Eulerian grid.
    //   cols: Number of cols in the Eulerian grid.
    //   diff: Diffusion coefficient for the passive markers.
    //   visc: Viscosity of the fluid. 
    StableFluidsSimWithMarkers( const int& rows, const int& cols, const ArrayXb & has_solid, const scalar& diff = 0.0001, const scalar& visc = 0.00000001 );
    
    // Deallocates memory used during the simulation.
    ~StableFluidsSimWithMarkers();
    
    // Integrates the system forward in time by dt. 
    void stepSystem( const scalar& dt );
    
    // Clears the density and velocity fields, and reseed particles
    void clear();
    
    // Free surface tracking
    const std::vector<Vector2s> & particles() const { return m_particles; }
    const ArrayXb & has_fluid() const { return m_has_fluid; }
    const ArrayXb & has_solid() const { return m_has_solid; }
    int fluid_cell_count() const { int count = 0; for (int i = 0; i < m_has_fluid.rows(); i++) for (int j = 0; j < m_has_fluid.cols(); j++) if (m_has_fluid(i, j)) count++; return count; }
    
protected:
    // Time stepping utilities
    virtual void dens_step(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar diff, scalar dt);
    virtual void vel_step(int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0, scalar visc, scalar dt);
    
    virtual void project(int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0);

    virtual void advect_particles(int N, ArrayXs & u, ArrayXs & v, scalar dt);
    virtual void detect_fluid(int N, ArrayXb & fluids, const std::vector<Vector2s> & particles);
    
    bool in_solid(const Vector2s & pos);    // note that pos.x is vertical coordinate, and pos.y is horizontal coordinate.
    
    // Implemenation for different solvers for the pressure solve
    void solveByJacobi      (ArrayXs & p, ArrayXs & div);
    void solveByGaussSeidel (ArrayXs & p, ArrayXs & div);
    void solveByEigenLU     (ArrayXs & p, ArrayXs & div);
    void solveByCG          (ArrayXs & p, ArrayXs & div);
    void solveByPCG         (ArrayXs & p, ArrayXs & div);
    
private:    
    // marker particles
    std::vector<Vector2s> m_particles;
    ArrayXb m_has_fluid;
    const ArrayXb m_has_solid;
    
    int m_original_fluid_cell_count;
    
};

#endif
