
# Constants
const ϵ₀ = 0.854e-12        # Permittivity of Free Space
const qₑ = 1.602e-19        # Elementary electron charge
const κ = 1.381e-23         # Boltzmann Constant
const AMU = 1.661e-27       # Atomic Mass Unit
const m_oxygen = 32*AMU     # mass of molecular oxygen
const m_ion = m_oxygen      # mass of an ion


function electrostatic_PIC()
    # Setup 
    n_0 = 1e12      # Density (1/m^3)
    ϕ_0 = 0         # Reference Potential (gauge field variant)
    T_e = 1         # Electron temp (eV)
    T_i = 0.1       # Ion Velocity in eV
    v_drift = 7000  # Ion injection velocity (7km/s)
    ϕ_p = -5        # Wall Potential

    # Calculate Plasma Parameters
    l_D = sqrt(ϵ₀ * T_e / (n_0 * qₑ))    	# Debye length
    v_th = sqrt(2 * qₑ * T_i / m_ion)		# Thermal velocity with T_i in eV
    ## Geometry Parts
    # Set simulation domain
    nx = 16               # number of nodes in x direction
    ny = 10               # number of nodes in y direction
    nt = 200              # number of time steps
    dh = l_D              # cell size
    np_insert = (ny-1) * 15 # insert 15 particles per cell

    nn = nx * ny            # total number of nodes
    dt = 0.1*dh / v_drift   # time step, at vdrift move 0.10dx
    L_x = (nx-1)*dh        # domain length in x direction
    L_y = (ny-1)*dh        # domain length in y direction

    # Specify plate dimensions (2x2 matrix)
    box = Array{Int}(undef, 2, 2)
    box[1, :] = [floor(nx / 3), floor(nx / 3)+2] # x range
    box[2, :] = [1, floor(ny / 2)];              # y range

    # create an object domain for visualization
    object = zeros(nx,ny);
    for j = box[2,1]:box[2,2]
        # For object values from x_min to x_max, and y_min to y_max, set to 1
        object[box[1,1]:box[1,2],j] = ones(box[1,2] - box[1,1] + 1, 1);
    end

    ## End Geometry

    # calculate specific weight
    flux = n_0*v_drift*L_y;     # flux of entering particles
    npt = flux*dt;              # number of real particles created per timestep
    spwt = npt/np_insert;       # specific weight, real particles per macroparticle
    mp_q = 1;                   # macroparticle charge
    max_part = 20000;           # buffer size
    
    # allocate particle array
    part_x = zeros(max_part, 2); # particle positions
    part_v = zeros(max_part, 2); # particle velocities

    # set up multiplication matrix for potential solver
    # here we are setting up the Finite Difference stencil

    A = zeros(nn, nn);              # allocate empty nn * nn matrix

    # set regular stencil on internal nodes
    for j = 2:ny-1                   # only internal nodes
        for i = 2:nx-1
            u = (j-1)*nx + i;         # unknown (row index)
            # ny = 10
            # nx = 16
            # j = range(2,9)
            # i = range(2,15)
            # u=18 -> problem. 16 + i(2) = 18.
            A[u,u]    = -4/(dh*dh);   # phi(i,j)
            A[u,u-1]  = 1/(dh*dh);    # phi(i-1,j)
            A[u,u+1]  = 1/(dh*dh);    # phi(i+1,j)
            A[u,u-nx] = 1/(dh*dh);    # phi(i,j-1)
            A[u,u+nx] = 1/(dh*dh);    # phi(i,j+1)
        end  
    end
    
    # neumann boundary on y=0
    for i = 1:nx
        u = i;
        A[u,u] = -1/dh;              # phi(i,j)
        A[u,u+nx] = 1/dh;            # phi(i,j+1)
    end
    
    # neumann boundary on y=L_y
    for i = 1:nx
        u = (ny-1)*nx + i;
        A[u,u-nx] = 1/dh;            # phi(i,j-1)
        A[u,u] = -1/dh;              # phi(i,j)
    end
    
    # neumann boundary on x=L_x
    for j = 1:ny
        u = (j-1)*nx + nx;
        A[u,:] = zeros(1,nn);       # clear row
        A[u,u-1] = 1/dh;            # phi(i-1,j)
        A[u,u] = -1/dh;             # phi(i,j)
    end
    
    # dirichlet boundary on x=0
    for j = 1:ny
        u = (j-1)*nx + 1;
        A[u,:] = zeros(1,nn);       # clear row
        A[u,u] = 1;                 # phi(i,j)
    end
    
    # dirichlet boundary on nodes corresponding to the plate
    for j = box[2,1]:box[2,2]
        for i = box[1,1]:box[1,2]
            u = (j-1)*nx + i;
            A[u,:] = zeros(1,nn);     # clear row
            A[u,u] = 1;               # phi(i,j)
        end
    end
    
    # initialize
    ϕ = ones(nx,ny) * ϕ_0;         # set initial potential to phi0
    np = 0;                        # clear number of particles
    
    print("Solving potential for the first time. Please be patient, this could take a while.")
    
    # Loop
        # Compute Charge Density
        # Compute E Potential
        # Compute E Field
        # Move Particles
        # Generate Particles
        # Output
        # Break
end
