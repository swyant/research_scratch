using Trapz
function calibrate(ecalib_pred, ecalib_ref, calib_uncertainty, alpha=0.05)
    num_calib = length(ecalib_pred)
    calib_scores = abs.(ecalib_pred .- ecalib_ref) ./ calib_uncertainty
    q_hat = quantile(calib_scores, clamp(ceil((num_calib+1)*(1-alpha))/num_calib, 0.0, 1.0))
    q_hat
end

function compute_miscalibration_area(expected_ps, observed_ps)
    area = 0.0
    #for i in 2:length(expected_ps)-1
    #    trap = abs(trapz(expected_ps[i-1:i+1], observed_ps[i-1:i+1]) -
    #             trapz(expected_ps[i-1:i+1], expected_ps[i-1:i+1]))
    for i in 2:length(expected_ps)
        trap = abs(trapz(expected_ps[i-1:i], observed_ps[i-1:i]) -
                 trapz(expected_ps[i-1:i], expected_ps[i-1:i]))
        area += trap
    end
    area
end

# converted from Medford jupyter notebook via Claude
function make_calibration_plot(expected_ps, observed_ps; width=600)
    # Convert to percentages
    expected_ps = expected_ps .* 100
    observed_ps = observed_ps .* 100

    fig = Figure(resolution=(width, width))
    ax = Axis(fig[1, 1],
        aspect=DataAspect(),
        xlabel="Expected conf. level",
        ylabel="Observed conf. level",
        limits=(0, 100, 0, 100)
    )

    # Main line
    lines!(ax, 1.0 .- expected_ps, observed_ps)

    # Diagonal reference line
    lines!(ax, 1.0 .-expected_ps, 1.0 .-expected_ps, linestyle=:dash, alpha=0.4)

    # Filled area between curves
    band!(ax, expected_ps, expected_ps, observed_ps, color=(:blue, 0.2))

    # Configure ticks - approximately 4 ticks on each axis
    ax.xticks = 0:10:100
    ax.yticks = 0:10:100

    # Add percentage signs to ticks
    ax.xtickformat = xs -> ["$(Int(x))%" for x in xs]
    ax.ytickformat = xs -> ["$(Int(x))%" for x in xs]

    ## Add text for miscalibration area
    #text!(ax, "miscalc. area = $(round(area, digits=3))",
    #    position=(8, 2),
    #    align=(:left, :bottom)
    #)

    return fig
end

# Claude
function parity_plot(etest_ref, etest_pred, qhat_scored;
                     title="Parity Plot",
                     xlabel="Reference Values",
                     ylabel="Predicted Values",
                     figsize=(600, 600))
    # Create figure and axis
    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1],
    title=title,
    xlabel=xlabel,
    ylabel=ylabel,
    limits = (-5.0,-4.0,-5.0,-4.0))

    # Calculate min and max for setting plot limits
    min_val = min(minimum(etest_pred), minimum(etest_ref))
    max_val = max(maximum(etest_pred), maximum(etest_ref))

    # Add diagonal reference line
    lines!(ax, [min_val, max_val], [min_val, max_val],
    color=:red,
    linestyle=:dash,
    label="Perfect Prediction")

    # Plot scatter with error bars
    errorbars!(ax, etest_ref, etest_pred, qhat_scored,
    whiskerwidth=1,  # Width of error bar caps
    color=:cyan3)

    # Scatter plot of points
    scatter!(ax, etest_ref, etest_pred,
    color=:teal,
    markersize=10)

    # Set equal aspect ratio
    #ax.aspect = DataAspect()

    # Add legend
    axislegend(ax)

    return fig
end

# Claude
function parity_plot(etest_ref, etest_pred, low_bounds, high_bounds;
                     title="Parity Plot",
                     xlabel="Reference Values",
                     ylabel="Predicted Values",
                     figsize=(600, 600))
    # Create figure and axis
    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1],
    title=title,
    xlabel=xlabel,
    ylabel=ylabel,
    limits = (-100.0,100.0,-100.0,100.0))

    # Calculate min and max for setting plot limits
    min_val = min(minimum(etest_pred), minimum(etest_ref))
    max_val = max(maximum(etest_pred), maximum(etest_ref))

    # Add diagonal reference line
    lines!(ax, [min_val, max_val], [min_val, max_val],
    color=:red,
    linestyle=:dash,
    label="Perfect Prediction")

    # Plot scatter with error bars
    errorbars!(ax, etest_ref, etest_pred, low_bounds, high_bounds,
    whiskerwidth=1,  # Width of error bar caps
    color=:cyan3)

    # Scatter plot of points
    scatter!(ax, etest_ref, etest_pred,
    color=:teal,
    markersize=10)

    # Set equal aspect ratio
    #ax.aspect = DataAspect()

    # Add legend
    axislegend(ax)

    return fig
end


function uncertainty_vs_residuals(uncertainty, residuals;
                                 title  ="Uncertainty vs. residuals",
                                 xlabel = "Distance",
                                 ylabel = "Residuals",
                                 figsize = (600,600),
                                 limits = (0,0.2,-0.0005,0.01))
    fig = Figure(size=figsize)
    ax  = Axis(fig[1,1],
               title=title,
               xlabel=xlabel,
               ylabel=ylabel,
               limits=limits)

    hlines!(ax, 0.0, color=:red, linestyle=:dash)

    scatter!(ax, uncertainty, residuals, markersize=3)

    #ax.aspect=DataAspect()
    fig
end

function generate_predicted_alphas(calib_scores, test_uq, test_abs_residuals)
    num_calib = length(calib_scores)
    num_test = length(test_uq)
    alpha_compls = collect(range(0.01,0.99,step=0.01))
    alphas = 1 .- alpha_compls # i.e. iterate 0.99..0.01, but will then plot as 1-0.99...1-0.01
    
    predicted_alphas = Float64[]
    for alpha in alphas
        qh = quantile(calib_scores, clamp(ceil((num_calib+1)*(1-alpha))/num_calib, 0.0, 1.0))
    
        qh_scores = qh*test_uq
        predicted_alpha = sum(test_abs_residuals .> qh_scores) / num_test
        push!(predicted_alphas, predicted_alpha)
    end
    predicted_alphas
end
#qhat_scores = q_hat*test_feature_distances[test_idxs_wrt_test]
#coverage = sum(test_abs_residuals .> qhat_scores) / num_test


# from subsampling_dpp.jl in PL.jl examples
function concat_dataset(confs::Vector{DataSet})
    N = length(confs)
    confs_vec = [[confs[i][j] for j = 1:length(confs[i])] for i = 1:N]
    confs_all = reduce(vcat, confs_vec)
    return DataSet(confs_all)
end
