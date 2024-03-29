#!/bin/env julia

# Code originally adapted from https://www.particleincell.com/2011/particle-in-cell-example/ (credit Lubos), by Chance Cardona in 2024

using Plots
using Printf

include("eval_2d_potential_GS.jl") # Solver for Maxwell's Equations (or Poisson's)
include("geometry.jl") # Geometry of Object / Environment
include("constants.jl") # Physical Constants


function electrostatic_PIC()
    """ Treats electrons are a fluid (Boltzmann Relation). 
        Doesnt account for magnetic or relativistic effects. """
    # Setup 
    n_0 = 1e12      # Equilibrium density of electrons at a potential (kg/m^3)
    ϕ_0 = 0         # Reference Potential (gauge field variant for the potential)
    T_e = 1         # Electron temp (eV)
    T_i = 0.1       # Ion temp (eV)
    v_drift = 7000  # Ion injection velocity (7km/s)
    ϕ_p = -5        # Wall Potential

    # Calculate Plasma Parameters
    l_D = sqrt(ϵ₀ * T_e / (n_0 * q_e))    	# Debye length
    v_th = sqrt(2 * q_e * T_i / m_ion)		# Thermal velocity with T_i in eV
    ## Geometry Parts
    # Set simulation domain
    nx = 16               # number of nodes in x direction
    ny = 10               # number of nodes in y direction
    nt = 200              # number of time steps
    dh = l_D              # cell size (cells are square)
    np_insert = (ny-1) * 15 # insert 15 particles per cell

    nn = nx * ny            # total number of nodes
    dt = 0.1*dh / v_drift   # time step, at vdrift move 0.10dx
    L_x = (nx-1)*dh        # domain length in x direction
    L_y = (ny-1)*dh        # domain length in y direction

    # Specify plate dimensions (2x2 matrix)
    plate_dims = Array{Int}(undef, 2, 2)
    plate_dims[1, :] = [floor(nx / 3), floor(nx / 3)+2] # x range
    plate_dims[2, :] = [1, floor(ny / 2)]              # y range

    # create an object domain for visualization
    object = rectangle_from_coords(plate_dims[1,1], plate_dims[1,2], plate_dims[2,1], plate_dims[2,2])

    ## End Geometry

    # calculate specific weight
    flux = n_0*v_drift*L_y     # flux of entering particles
    npt = flux*dt              # number of real particles created per timestep
    sp_wt = npt/np_insert       # specific weight, real particles per macroparticle
    q_mp = 1                   # macroparticle charge
    max_part = 20000           # buffer size
    
    # allocate particle array
    part_x = zeros(max_part, 2) # particle positions
    part_v = zeros(max_part, 2) # particle velocities

    # set up multiplication matrix for potential solver
    # (here we are setting up the Finite Difference stencil)
    global A = zeros(nn, nn)              # allocate empty nn * nn matrix

    # set regular stencil on internal nodes
    for j = 2:ny-1                   # only internal nodes
        for i = 2:nx-1
            u = (j-1)*nx + i         # unknown (row index)
            # ny = 10
            # nx = 16
            # j = range(2,9)
            # i = range(2,15)
            # u=18 -> problem. 16 + i(2) = 18.
            A[u,u]    = -4/(dh*dh)   # phi(i,j)
            A[u,u-1]  = 1/(dh*dh)    # phi(i-1,j)
            A[u,u+1]  = 1/(dh*dh)    # phi(i+1,j)
            A[u,u-nx] = 1/(dh*dh)    # phi(i,j-1)
            A[u,u+nx] = 1/(dh*dh)    # phi(i,j+1)
        end  
    end
    
    # neumann boundary on y=0
    for i = 1:nx
        u = i
        A[u,u] = -1/dh              # phi(i,j)
        A[u,u+nx] = 1/dh            # phi(i,j+1)
    end
    
    # neumann boundary on y=L_y
    for i = 1:nx
        u = (ny-1)*nx + i
        A[u,u-nx] = 1/dh            # phi(i,j-1)
        A[u,u] = -1/dh              # phi(i,j)
    end
    
    # neumann boundary on x=L_x
    for j = 1:ny
        u = (j-1)*nx + nx
        A[u,:] = zeros(1,nn)       # clear row
        A[u,u-1] = 1/dh            # phi(i-1,j)
        A[u,u] = -1/dh             # phi(i,j)
    end
    
    # dirichlet boundary on x=0
    for j = 1:ny
        u = (j-1)*nx + 1
        A[u,:] = zeros(1,nn)       # clear row
        A[u,u] = 1                 # phi(i,j)
    end
    
    # dirichlet boundary on nodes corresponding to the plate
    for j = plate_dims[2,1]:plate_dims[2,2]
        for i = plate_dims[1,1]:plate_dims[1,2]
            u = (j-1) * nx + i
            A[u,:] = zeros(1,nn)     # clear row
            A[u,u] = 1               # phi(i,j)
        end
    end
    
    # initialize
    ϕ = ones(nx,ny) * ϕ_0         # set initial potential to phi0
    np = 0                        # clear number of particles
    
    print("Solving potential for the first time. Please be patient, this could take a while.")
    
    # Loop through all times
    for t = 1:nt
        # reset field quantities
        n_d = zeros(nx, ny) # Number density
        E_x = zeros(nx, ny) # E field x component
        E_y = zeros(nx, ny) # E field y component
        q_dist = zeros(nx, ny) # Charge distribution

        # 1. Compute Charge Density   	
        # Charge density is a scalar field that varies across space.
        # Compute it for each infinitesimal mesh cell volume by summing the particles
        # in the volume, then distributing it to the 4 (if 2D) (8 if 3D) actual mesh-grid nodes proportional
        # to how near the particle is to the corresponding node.
	    for p = 1:np                         # loop over particles
            x_ind = 1 + (part_x[p,1] / dh)      # x (fractional included) index of particle's cell
            i::Int = floor(x_ind)
	    	hx = x_ind - i                    # the remainder
            
            y_ind = 1 + (part_x[p,2] / dh)      # y (fractional included) index of particle's cell
            j::Int = floor(y_ind)
            hy = y_ind - j                    # the remainder

            # interpolate charge to nodes via the Scatter Operation
	    	q_dist[i, j] = q_dist[i, j] + (1-hx) * (1-hy)
	    	q_dist[i+1, j] = q_dist[i+1, j] + hx * (1-hy)
            q_dist[i, j+1] = q_dist[i, j+1] + (1-hx) * hy
            q_dist[i+1, j+1] = q_dist[i+1, j+1] + hx * hy
	    end 

	    # calculate density
	    n_d = sp_wt * q_mp * q_dist / (dh^2)
        
        # apply boundaries
        n_d[1,:] = 2 * n_d[1,:]      # double density since only half volume contributing
        n_d[nx,:] = 2 * n_d[nx,:]
        n_d[:,1] = 2 * n_d[:,1]
        n_d[:,ny] = 2 * n_d[:,ny]
        
        # add density floor for plotting and to help the solver
        n_d .+= 1e4

        # 2. Compute E Potential
        ϕ = eval_2d_potential_GS(ϕ, ϕ_0, ϕ_p, n_d, n_0, T_e, plate_dims)

        # 3. Compute E Field
        E_x[2:nx-1,:] = ϕ[1:nx-2,:] - ϕ[3:nx,:] # Central Difference on internal nodes X
        E_y[:,2:ny-1] = ϕ[:,1:ny-2] - ϕ[:,3:ny] # Central Difference on internal nodes Y
        E_x[1,:] = 2*(ϕ[1,:] - ϕ[2,:]) # Forward Difference on X=0
        E_y[:,1] = 2*(ϕ[:,1] - ϕ[:,2]) # Forward Difference on Y=0
        E_x[nx,:] = 2*(ϕ[nx-1,:] - ϕ[nx,:]) # Backward Difference on X=L_x
        E_y[:,ny] = 2*(ϕ[:,ny-1] - ϕ[:,ny]) # Backward Difference on Y=L_y
        E_x = E_x / (2*dh) # Divide by nominator
        E_y = E_y / (2*dh)

        # 4. Generate Particles
        # Check array limits
        if (np + np_insert >= max_part)
            np_insert = max_part - np
        end

        # Assumes Birdsall approximation of maxwellian sampling distribution
        part_x[np+1 : np+np_insert, 1] = rand(np_insert,1)*dh # x position (within leftmost cell column)
        part_x[np+1 : np+np_insert, 2] = rand(np_insert,1)*L_y # y position (anywhere in mesh y)
        M = 1
        part_v[np+1 : np+np_insert, 1] = v_drift .+ (-1.5 .+ rand(np_insert, 1) + rand(np_insert, 1) + rand(np_insert, 1)) * v_th
        part_v[np+1 : np+np_insert, 2] = sqrt(M/12) .* ( sum( rand(np_insert, 1) for _ in 1:M ) .- M/2 ) * v_th

        np += np_insert
 
        # 5. Move Particles
        p = 1
        while (p <= np)
            f_i = 1 + part_x[p, 1] / dh
            i::Int = floor(f_i) # i index of particle's cell
            hx = f_i - i # fractional x position in cell
            f_j = 1 + part_x[p, 2] / dh
            j::Int = floor(f_j) # j index of particle's cell
            hy = f_j - j # fractional y position in cell

            # Gather Exlectric Field
            E = [E_x[i, j], E_y[i, j]] * (1-hx) * (1-hy) # (i,j) contribution
            E += [E_x[i+1, j], E_y[i+1, j]] * hx * (1-hy) # (i+1,j) contribution
            E += [E_x[i, j+1], E_y[i, j+1]] * hx * (1-hy) # (i, j+1) contribution
            E += [E_x[i+1, j+1], E_y[i+1, j+1]] * hx * (1-hy) # (i+1, j+1) contribution

            # Calculate Position and Velocity from Force
            F = q_e * E # Lorentz Force F = qE
            a = F / m_ion
            part_v[p,:] = part_v[p, :] + a * dt
            part_x[p,:] = part_x[p, :] + part_v[p, :] * dt

            # Process Boundaries (where is particle at)
            #   Relective Boundary (on bottom)
            if (part_x[p, 2] < 0) # y < 0
                part_x[p, 2] = -part_x[p, 2] # Flip position back into domain
                part_v[p, 2] = -part_v[p, 2] # Reverse y velocity
            end
            #   Inside Plate (our obstacle object)
            in_box = false
            if ((i >= plate_dims[1,1] && i < plate_dims[1, 2]) &&
                (j >= plate_dims[2,1] && j < plate_dims[2, 2]))
                in_box = true
            end  
            #   Absorbing Boundary (left, right, top, or in object)
            if (part_x[p, 1] < 0 || 
                part_x[p, 1] >= L_x ||
                part_x[p, 2] >= L_y ||
                in_box)
                part_x[p, :] = part_x[np, :] # Kill particle by replacing it with the last particle in particle array
                part_v[p, :] = part_v[np, :] # Kill particle v by replacing it with the last particle v in particle array
                np -= 1 # Reduce particle count
                p -= 1 # Reduce particle index so this entry is processed again
            end
            p += 1 # Move to next particle
        end

        # 6. Output (Plot)
        if (mod(t,25) == 0 || t == 1 || t == nt)      # plot only every 25 time steps, or the first/last.
            # Density Plot
            density_ticks = collect(1e11:1e11:1.1e12)
            p_density = contour(n_d', 
                                fill=true, 
                                colorbar=true, clims=(1e11, 1.1e12), 
                                # TODO: add colorbar formatting for GR
                                #colorbar_ticks=(density_ticks, [ @sprintf("%.1e", tick) for tick in density_ticks ]),
                                title=@sprintf("Density %i", t))
            # Add the geometry-object outline
            plot!(object[:,2], object[:,1], linecolor=:black, linewidth=2, legend=false)
            # Potential Plot
            p_potential = contour(ϕ', fill=true, colorbar=true, title=@sprintf("Potential %i", t))
            # Add the geometry-object outline
            plot!(object[:,2], object[:,1], linecolor=:black, linewidth=2, legend=false)

            # Combine and show
            combined_plot = plot(p_density, p_potential, layout=(1,2))
            # Draw
            display(combined_plot)
            sleep(1.5)
        end
        
        print("Time Step $t (t=$(t*dt)s), Particles $np:")
    end # finish iterating through times
    print("Complete!\n")
    sleep(15)
end
