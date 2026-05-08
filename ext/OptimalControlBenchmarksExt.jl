using OptimalControlBenchmarks
using BlockSQP2
from BlockSQP2.NLPlayouts import get_layout

function solve_with_corleone(benchmark, optimizer::BlockSQP2.Optimizer, grids, plot_res = true)

    data = benchmark(grids)

    oc_problem = data.system
    scaled_grids = data.grids
    num_states, num_controls = data.dims

    # Extract control variable
    controls = inputs(oc_problem)
    control_grid = scaled_grids.control_grid
    shooting_grid = scaled_grids.shooting_grid

    control_map = [
        c => control_grid for c in controls
    ]

    dynopt = CorleoneDynamicOptProblem(
        oc_problem,
        [],
        control_map...;
        algorithm = Tsit5(),
        shooting = shooting_grid
    )
    
    layout = get_layout(dynopt.layer)
    
    optprob = OptimizationProblem(
        dynopt,
        AutoForwardDiff(),
        Val(:ComponentArrays)
    )

    sol = solve(
        optprob,
        optimizer,
        max_iter = 100,
        opt_tol = 5.0e-6,
        layout = layout
    )

    if plot_res
        # plotting
        u_opt = ComponentVector(sol.u, optprob.f.f.ax)
        opt_traj, _ = dynopt.layer(
            nothing, u_opt,
            LuxCore.initialstates(Random.default_rng(), dynopt.layer)
        )

        f = plot_oc_problem(opt_traj, num_states, num_controls)
        display(f)
    end

    return sol

end