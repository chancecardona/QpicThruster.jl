# Potential solver for a particle-in-cell example program
# Based on the Gauss-Seidel method
#
# For more, visit http://www.particleincell.com/2010/es-pic-method/
# and http://www.particleincell.com/2011/particle-in-cell-example/
###################################################################

using LinearAlgebra

# ϕ is the potential
# ϕ₀ is the reference potential
# ϕₚ is the wall potential
# n_d is the number density of the particles
# n_0 is the density of air (kg/m^3)
# T_e is the temperature of the electrons (eV)
# plate_dims is the dimensions
# A is the Finite Difference Stencil and should be passed by const ref.
function eval_2d_potential_GS(ϕ, ϕ_0, ϕ_p, n_d, n_0, T_e, plate_dims)
    global A, ϵ₀, q_e
    
    tol = 0.1      # solver tolerance
    
    # get nx from size of density
    nx = size(n_d,1)
    ny = size(n_d,2)
    nn = length(n_d)
    
    # convert density and potential into column vectors
    # TODO: better way of doing this?
    b_0 = reshape(n_d, length(n_d),1)
    x = reshape(ϕ, length(ϕ),1)
  
    # Residue of error
    R = 0.0

    # solve
    for t = 1:2000
        
        # recalculate rhs
        b = b_0 .- n_0 * exp.( (x .- ϕ_0) / T_e )     # add boltzmann term for electrons
        b = -b * q_e / ϵ₀
     
        # set boundaries
        b[1:nx] .= 0                 # zero electric field on y=0
        b[nn-nx+1:nn] .= 0           # zero electric field on y=L
        b[nx:nx:nn] .= 0             # zero electric field on x=L
        b[1:nx:nn] .= ϕ_0           # fixed potential on x=0
    
        # set potential on fixed nodes
        for j = plate_dims[2,1]:plate_dims[2,2]
            b[[plate_dims[1,1] : plate_dims[1,2]] .+ (j-1) * nx] = ones(plate_dims[1,2] - plate_dims[1,1] + 1, 1) .* ϕ_p      # wall potential
        end
    
        # update nodes
    	for i = 1:nn
            x[i] = (b[i] - (A[i,1:i-1]' * x[1:i-1]) - (A[i,i+1:nn]' * x[i+1:nn])) / A[i,i]
        end
        
        # compute residue to check for convergence, do only every 10 iterations
        if mod(t, 10) == 0
            R = norm(b - A*x)       # set the residue
            if (R <= tol)           # converged
                print("  GS converged in $t iterations with norm $R\n")
                break
            end
        end
    end
    
    # check if the solver converged to the specified tolerance
    if (R > tol)
        print("  GS failed to converge!!\n")
    end
    
    # return solution as a nx*ny array
    x = reshape(x, nx, ny)
    return x
end
