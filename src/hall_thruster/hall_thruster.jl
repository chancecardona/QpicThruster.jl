# This uses https://docs.juliahub.com/General/HallThruster/stable/physics/ for now. 
# This just documents the use of this external code base.

using Revise, HallThruster
using Plots

function run_landmark(duration = 1e-3; ncells = 200, nsave = 2, dt = 0.7e-8, CFL = 0.8, case = 1)

    domain = (0.0, 0.05)

    #Landmark cases loss frequencies
    αϵ_in, αϵ_out = if case == 1
        (1.0, 1.0)
    elseif case == 2
        (0.5, 1.0)
    elseif case == 3
        (0.4, 1.0)
    end

    scheme = HallThruster.HyperbolicScheme(
        flux_function = HallThruster.rusanov,
        limiter = HallThruster.minmod,
        reconstruct = true
    )

    ϵ_anode = 3.0
    ϵ_cathode = 3.0

    # Define the length over which the anomalous transport, wall collision frequency, and wall loss coefficients change
    transition_function = HallThruster.LinearTransition(1e-3, 0.0)

    config = HallThruster.Config(;
        ncharge = 1,
        anode_Te = 2/3 * ϵ_anode,
        cathode_Te = 2/3 * ϵ_cathode,
        discharge_voltage = 300.0,
        ionization_model = HallThruster.LandmarkIonizationLookup(),
        excitation_model = HallThruster.LandmarkExcitationLookup(),
        electron_neutral_model = HallThruster.LandmarkElectronNeutral(),
        electron_ion_collisions = false,
        wall_loss_model = HallThruster.ConstantSheathPotential(20, αϵ_in, αϵ_out),
        LANDMARK = true,
        neutral_velocity = 150.0,
        neutral_temperature = 0.0,
        ion_temperature = 0.0,
        thruster = HallThruster.SPT_100,
        anode_mass_flow_rate = 5e-6,
        scheme,
        domain,
        transition_function,
        electron_pressure_coupled = true,
        ion_wall_losses = false,
        anom_model = HallThruster.TwoZoneBohm(1/160, 1/16),
        anode_boundary_condition = :dirichlet,
    )

    @time sol = HallThruster.run_simulation(
        config; duration, grid = EvenGrid(ncells), nsave,
        dt, dtmin = dt / 100, dtmax = dt * 100, adaptive = true, CFL
    )
    return sol
end

sol_case1 = run_landmark(1e-3; ncells=150, nsave=1000, case = 1, CFL = 0.9)

p = plot(sol_case1, label = "HallThruster.jl", legend = :outertop);
display(p)

landmark_1, landmark_2, landmark_hybrid = HallThruster.load_landmark_data(1);


# Compare to landmark

function compare_to_landmark(sol, case)
    p = plot(sol, label = "HallThruster.jl", legend = :outertop, lc = :black)

    landmark_1, landmark_2, landmark_hybrid = HallThruster.load_landmark_data(case)
    plot!(p, landmark_1, ls = :dash, lc = RGB(0.0, 1.0, 0.0), label = "LANDMARK, fluid, δ = 0.5 mm")
    plot!(p, landmark_2, ls = :dash, lc = RGB(1.0, 0.0, 0.0), label = "LANDMARK, fluid, δ = 1.0 mm")
    plot!(p, landmark_hybrid, ls = :dash, lc = RGB(0.0, 0.0, 1.0), label = "LANDMARK, Hybrid-PIC")
    return p
end

p = compare_to_landmark(sol_case1, 1); display(p)

avg_case1 = time_average(sol_case1)
p = compare_to_landmark(avg_case1, 1); display(p)

sol_case2 = run_landmark(1e-3; ncells=200, nsave=1000, case = 2, dt = 1.0e-8, CFL = 0.9)
avg_case2 = time_average(sol_case2)
p = compare_to_landmark(avg_case2, 2); display(p)

sol_case3 = run_landmark(1e-3; ncells=200, nsave=1000, case = 3, dt = 1.0e-8, CFL = 0.9)
avg_case3 = time_average(sol_case3)
p = compare_to_landmark(avg_case3, 3); display(p)

# Postprocess to get Thrust and other relevant values
function plot_current(sol)
    (;A_ch, index, config) = sol.params
    e = HallThruster.e
    mi = config.propellant.m
    t_ms = sol.t * 1000 # convert from seconds to milliseconds
    N = length(t_ms)

    current = HallThruster.discharge_current(sol)
    ion_current = HallThruster.ion_current(sol)
    electron_current = HallThruster.electron_current(sol)

    p = plot(t_ms, current, label = "Discharge current", xlabel = "Time (ms)", ylabel = "Current (A)", linewidth = 1.5, ylims = (-10, 25))
    plot!(p, t_ms, ion_current, label = "Ion current (cathode)")
    plot!(p, t_ms, electron_current, label = "Electron current (cathode)")
    return p
end

p = plot_current(sol_case1); display(p)

p = plot_current(sol_case2); display(p)

p = plot_current(sol_case3); display(p)
