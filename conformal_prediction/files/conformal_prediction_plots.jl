using ColorSchemes

#### Parity Plot 
function custom_parity_plot(etest_ref, etest_pred, qhat_scored;
                            title="Parity Plot Subset",
                            xlabel="Reference Energies (eV)",
                            ylabel="Predicted Energies(eV)",
                            width=600,
                            colormap=:viridis,
                            color_value=0.6,
                            marker_size=10,
                            line_width=3.0,
                            text_size=18,
                            label_size=22,
                            grid_visible=false,
                            grid_color=(:gray, 0.3),
                            grid_linewidth=0.5,
                            errorbar_color=nothing,
                            marker_color=nothing,
                            diagonal_color=:red,
                            diagonal_alpha=0.6,
                            diagonal_style=:dash,
                            equal_aspect=false)

    # Get color from colormap if no specific colors provided
    diagonal_color = base_color = get(ColorSchemes.colorschemes[colormap], color_value)
    marker_color = isnothing(marker_color) ?
                    get(ColorSchemes.colorschemes[colormap], color_value+0.15) : marker_color
    errorbar_color = isnothing(marker_color) ?
                    get(ColorSchemes.colorschemes[colormap], color_value+0.25) : marker_color
    
    
    # Create figure and axis with better formatting
    fig = Figure(resolution=(1.75*width, width), fontsize=text_size, figure_padding=30)
    
    # Calculate min and max for setting plot limits
    min_val = min(minimum(etest_pred), minimum(etest_ref))
    max_val = max(maximum(etest_pred), maximum(etest_ref))
    
    min_val = -3257.35
    max_val = -3257.15
    
    ## Add a small buffer to the limits to avoid cutting off points or labels
    #buffer = (max_val - min_val) * 0.05
    #plot_min = min_val - buffer
    #plot_max = max_val + buffer
    
    ax = Axis(fig[1, 1],
              #title=title,
              xlabel=xlabel,
              ylabel=ylabel,
              #limits=(plot_min, plot_max, plot_min, plot_max),
              limits=(min_val, max_val, min_val-0.5, max_val+0.5),
              titlesize=label_size,
              xlabelsize=label_size,
              ylabelsize=label_size,
              xticklabelsize=text_size,
              yticklabelsize=text_size,
              spinewidth=1.5,
              xgridvisible=grid_visible,
              ygridvisible=grid_visible,
              xgridcolor=grid_color,
              ygridcolor=grid_color,
              xgridwidth=grid_linewidth,
              ygridwidth=grid_linewidth
              )

     
    # Add diagonal reference line
    #lines!(ax, [plot_min, plot_max], [plot_min, plot_max],
    lines!(ax, [min_val,max_val], [min_val,max_val], 
           color=diagonal_color,
           linestyle=diagonal_style,
           linewidth=line_width - 1,  # Slightly thinner than main points
           alpha=diagonal_alpha,
           label="Perfect Prediction")
    
    # Plot scatter with error bars
    errorbars!(ax, etest_ref, etest_pred, qhat_scored,
               whiskerwidth=6,  # Width of error bar caps
               color=errorbar_color)
    
    # Scatter plot of points
    scatter!(ax, etest_ref, etest_pred,
             color=marker_color,
             markersize=marker_size)
    
    # Set equal aspect ratio (usually important for parity plots)
    if equal_aspect
        ax.aspect = DataAspect()
    end
    
    # Add legend with better formatting
    #axislegend(ax, position=:lt, framevisible=true, framecolor=(:black, 0.2),
    #padding=(10, 10, 10, 10), labelsize=text_size-2)
    
    return fig
end

##### Calibration Plots 


# Initial Plots
function make_custom_calibration_plot(expected_ps, observed_ps;
                                      width=600,
                                      colormap=:viridis,
                                      color_value=0.6,  # Value between 0-1 in the colormap
                                      main_line_width=3.0,
                                      band_alpha=0.2,
                                      axis_color=:black,
                                      text_size=18,
                                      label_size=22,
                                      grid_visible=true,
                                      grid_color=(:gray, 0.3),
                                      grid_width=0.5)
    # Convert to percentages
    #expected_ps = expected_ps .* 100
    #observed_ps = observed_ps .* 100

    expected_ps = (1.0 .- expected_ps).* 100
    observed_ps = (1.0 .- observed_ps).* 100
    # Get color from colormap
    color = get(ColorSchemes.colorschemes[colormap], color_value)
    band_color = (color, band_alpha)

    fig = Figure(resolution=(width, width), fontsize=text_size, figure_padding=30)
    ax = Axis(fig[1, 1],
        aspect=DataAspect(),
        xlabel="Expected Confidence Level",
        ylabel="Observed Confidence Level",
        limits=(0, 100, 0, 100),
        xlabelsize=label_size,
        ylabelsize=label_size,
        xticklabelsize=text_size,
        yticklabelsize=text_size,
        spinewidth=1.5,
        xgridvisible=grid_visible,
        ygridvisible=grid_visible,
        xgridcolor=grid_color,
        ygridcolor=grid_color,
        xgridwidth=grid_width,
        ygridwidth=grid_width
    )

    #Set spine and tick colors
    ax.bottomspinecolor = axis_color
    ax.leftspinecolor = axis_color
    ax.rightspinecolor = axis_color
    ax.topspinecolor = axis_color

    ax.xticklabelcolor = axis_color
    ax.yticklabelcolor = axis_color
    ax.xlabelcolor = axis_color
    ax.ylabelcolor = axis_color

    # Main line - made bolder
    lines!(ax, expected_ps, observed_ps, color=:black, linewidth=main_line_width)

    # Diagonal reference line
    lines!(ax, expected_ps, expected_ps, linestyle=:dash, color=:black, alpha=0.6, linewidth=1.5)

    # Filled area between curves
    band!(ax, expected_ps, expected_ps, observed_ps, color=band_color)
    #band!(ax, expected_ps, expected_ps, observed_ps, color=(:blue, 0.2))

    # Configure ticks
    ax.xticks = 0:20:100
    ax.yticks = 0:20:100

    # Add percentage signs to ticks
    ax.xtickformat = xs -> ["$(Int(x))%" for x in xs]
    ax.ytickformat = xs -> ["$(Int(x))%" for x in xs]

    return fig
end

# from ediff plots
function make_custom_calibration_plot1(expected_ps, observed_ps;
                                      width=600,
                                      colormap=:viridis,
                                      color_value=0.6,  # Value between 0-1 in the colormap
                                      main_line_width=3.0,
                                      band_alpha=0.2,
                                      axis_color=:black,
                                      text_size=18,
                                      label_size=22,
                                      grid_visible=true,
                                      grid_color=(:gray, 0.3),
                                      grid_width=0.5)
    # Convert to percentages
    #expected_ps = expected_ps .* 100
    #observed_ps = observed_ps .* 100

    expected_ps = (1.0 .- expected_ps).* 100
    observed_ps = (1.0 .- observed_ps).* 100

    # Get color from colormap
    colormap = :managua
    #axis_color = get(ColorSchemes.colorschemes[colormap], 0.4)
    #grid_color = (axis_color, 0.3)
    base_band_color = get(ColorSchemes.colorschemes[colormap], 0.5)
    band_color = (base_band_color, band_alpha)

    line_color = get(ColorSchemes.colorschemes[colormap], 0.5)

    # Get color from colormap
    #color = get(ColorSchemes.colorschemes[colormap], color_value)
    #band_color = (color, band_alpha)
    #line_color=:black

    fig = Figure(resolution=(width, width), fontsize=text_size, figure_padding=30)
    ax = Axis(fig[1, 1],
        aspect=DataAspect(),
        xlabel="Expected Confidence Level",
        ylabel="Observed Confidence Level",
        limits=(0, 100, 0, 100),
        xlabelsize=label_size,
        ylabelsize=label_size,
        xticklabelsize=text_size,
        yticklabelsize=text_size,
        spinewidth=1.5,
        xgridvisible=grid_visible,
        ygridvisible=grid_visible,
        xgridcolor=grid_color,
        ygridcolor=grid_color,
        xgridwidth=grid_width,
        ygridwidth=grid_width
    )

    #Set spine and tick colors
    ax.bottomspinecolor = axis_color
    ax.leftspinecolor = axis_color
    ax.rightspinecolor = axis_color
    ax.topspinecolor = axis_color

    ax.xticklabelcolor = axis_color
    ax.yticklabelcolor = axis_color
    ax.xlabelcolor = axis_color
    ax.ylabelcolor = axis_color

    # Main line - made bolder
    lines!(ax, expected_ps, observed_ps, color=line_color, linewidth=main_line_width)

    # Diagonal reference line
    lines!(ax, expected_ps, expected_ps, linestyle=:dash, color=line_color, alpha=0.6, linewidth=1.5)

    # Filled area between curves
    band!(ax, expected_ps, expected_ps, observed_ps, color=band_color)
    #band!(ax, expected_ps, expected_ps, observed_ps, color=(:blue, 0.2))

    # Configure ticks
    ax.xticks = 0:20:100
    ax.yticks = 0:20:100

    # Add percentage signs to ticks
    ax.xtickformat = xs -> ["$(Int(x))%" for x in xs]
    ax.ytickformat = xs -> ["$(Int(x))%" for x in xs]

    return fig
end

using ColorSchemes

function make_custom_calibration_plot2(expected_ps, observed_ps;
                                      width=600,
                                      colormap=:viridis,
                                      color_value=0.6,  # Value between 0-1 in the colormap
                                      main_line_width=3.0,
                                      band_alpha=0.2,
                                      axis_color=:black,
                                      text_size=18,
                                      label_size=22,
                                      grid_visible=true,
                                      grid_color=(:gray, 0.3),
                                      grid_width=0.5)
    # Convert to percentages
    #expected_ps = expected_ps .* 100
    #observed_ps = observed_ps .* 100

    expected_ps = (1.0 .- expected_ps).* 100
    observed_ps = (1.0 .- observed_ps).* 100

    # Get color from colormap
    colormap = :managua
    #axis_color = get(ColorSchemes.colorschemes[colormap], 0.4)
    #grid_color = (axis_color, 0.3)
    base_band_color = get(ColorSchemes.colorschemes[colormap], 0.35)
    band_color = (base_band_color, band_alpha)

    line_color = get(ColorSchemes.colorschemes[colormap], 0.35)

    # Get color from colormap
    #color = get(ColorSchemes.colorschemes[colormap], color_value)
    #band_color = (color, band_alpha)
    #line_color=:black

    fig = Figure(resolution=(width, width), fontsize=text_size, figure_padding=40)
    ax = Axis(fig[1, 1],
        aspect=DataAspect(),
        #xlabel="Expected Confidence Level",
        #ylabel="Observed Confidence Level",
        limits=(0, 100, 0, 100),
        xlabelsize=label_size,
        ylabelsize=label_size,
        xticklabelsize=text_size,
        yticklabelsize=text_size,
        spinewidth=1.5,
        xgridvisible=grid_visible,
        ygridvisible=grid_visible,
        xgridcolor=grid_color,
        ygridcolor=grid_color,
        xgridwidth=grid_width,
        ygridwidth=grid_width
    )

    #Set spine and tick colors
    ax.bottomspinecolor = axis_color
    ax.leftspinecolor = axis_color
    ax.rightspinecolor = axis_color
    ax.topspinecolor = axis_color

    ax.xticklabelcolor = axis_color
    ax.yticklabelcolor = axis_color
    ax.xlabelcolor = axis_color
    ax.ylabelcolor = axis_color

    # Main line - made bolder
    lines!(ax, expected_ps, observed_ps, color=line_color, linewidth=main_line_width)

    # Diagonal reference line
    lines!(ax, expected_ps, expected_ps, linestyle=:dash, color=line_color, alpha=0.6, linewidth=1.5)

    # Filled area between curves
    band!(ax, expected_ps, expected_ps, observed_ps, color=band_color)
    #band!(ax, expected_ps, expected_ps, observed_ps, color=(:blue, 0.2))

    # Configure ticks
    ax.xticks = 0:20:100
    ax.yticks = 0:20:100

    # Add percentage signs to ticks
    ax.xtickformat = xs -> ["$(Int(x))%" for x in xs]
    ax.ytickformat = xs -> ["$(Int(x))%" for x in xs]

    return fig
end


#### Histograms
function custom_histogram(data;
    width=600,
    bins=500,
    colormap=:viridis,
    color_value=0.6,
    title="Histogram",
    xlabel="Value",
    ylabel="Frequency",
    fill_alpha=0.8,
    edge_linewidth=1.0,
    axis_color=:black,
    text_size=18,
    label_size=22,
    grid_visible=false,
    grid_color=(:gray, 0.3),
    grid_linewidth=0.5,
    bar_color=nothing,
    edge_color=nothing,
    normalize=false,
    kde=false,
    kde_linewidth=3.0,
    kde_color=:black)

# Get color from colormap if no specific colors provided
base_color = get(ColorSchemes.colorschemes[colormap], color_value)
bar_color = isnothing(bar_color) ? (base_color, fill_alpha) : bar_color
#edge_color = isnothing(edge_color) ? darker(base_color, 0.2) : edge_color

# Create figure and axis with better formatting
fig = Figure(resolution=(width, width), fontsize=text_size)

# Calculate sensible limits with buffer
data_min = minimum(data)
data_max = maximum(data)
#buffer = (data_max - data_min) * 0.05
#x_min = data_min - buffer
#x_max = data_max + buffer
#x_min = -0.01
#x_max = 1.0
x_min= data_min
x_max = data_max

# Create axis with formatting
ax = Axis(fig[1, 1],
#title=title,
xlabel=xlabel,
ylabel=ylabel,
xlabelsize=label_size,
ylabelsize=label_size,
titlesize=label_size,
xticklabelsize=text_size,
yticklabelsize=text_size,
spinewidth=1.5,
xgridvisible=grid_visible,
ygridvisible=grid_visible,
xgridcolor=grid_color,
ygridcolor=grid_color,
xgridwidth=grid_linewidth,
ygridwidth=grid_linewidth
)

# Set spine and tick colors

ax.bottomspinecolor = axis_color
ax.leftspinecolor = axis_color
ax.rightspinecolor = axis_color
ax.topspinecolor = axis_color

#ax.xticks = 0:0.2:1.0
#ax.yticks = 0:10:50

ax.xticklabelcolor = axis_color
ax.yticklabelcolor = axis_color
ax.xlabelcolor = axis_color
ax.ylabelcolor = axis_color
ax.titlecolor = axis_color

# Add extra padding to avoid cutting off labels
#fig.margin = 20

# Create the histogram
hist = hist!(ax, data,
bins=bins,
color=bar_color,
#strokecolor=edge_color,
strokecolor=bar_color,
strokewidth=edge_linewidth,
normalization=normalize ? :pdf : :none)

# Optionally add KDE curve
if kde
density = kde!(ax, data,
color=kde_color,
linewidth=kde_linewidth,
label="KDE")

# Add legend if KDE is used
axislegend(ax, position=:rt, framevisible=true,
framecolor=(:black, 0.2),
padding=(10, 10, 10, 10),
labelsize=text_size-2)
end

# Adjust x limits
ax.limits = (x_min, x_max, nothing, nothing)

return fig
end


function custom_histogram2(data1, data2;
                           width=600,
                           bins=10,
                           colormap=:viridis,
                           color_value=0.0,
                           title="Histogram",
                           xlabel="Score Values",
                           ylabel="Frequency",
                           fill_alpha=0.5,
                           edge_linewidth=1.0,
                           axis_color=:black,
                           text_size=18,
                           label_size=22,
                           grid_visible=false,
                           grid_color=(:gray, 0.3),
                           grid_linewidth=0.5,
                           bar_color=nothing,
                           edge_color=nothing,
                           normalize=false,
                           kde=false,
                           kde_linewidth=3.0,
                           kde_color=:black)

    # Get color from colormap if no specific colors provided
    base_color1 = get(ColorSchemes.colorschemes[colormap], color_value)
    bar_color1 = isnothing(bar_color) ? (base_color1, fill_alpha) : bar_color
    
    base_color2 = get(ColorSchemes.colorschemes[colormap], color_value+0.3)
    bar_color2 = isnothing(bar_color) ? (base_color2, 0.8) : bar_color
    
    #base_color1 = :black
    #bar_color1 = (:black, fill_alpha)
    #
    #base_color2 = :black
    #bar_color2 = (:black, fill_alpha)
    
    #edge_color = isnothing(edge_color) ? darker(base_color, 0.2) : edge_color
    
    # Create figure and axis with better formatting
    fig = Figure(resolution=(width, width), fontsize=text_size)
    
    # Calculate sensible limits with buffer
    #data1 = [datum for datum in data1 if datum < 2]
    #data2 = [datum for datum in data2 if datum < 2]
    data_min = minimum([data1;data2])
    data_max = maximum([data1;data2])
    #buffer = (data_max - data_min) * 0.05
    #x_min = data_min - buffer
    #x_max = data_max + buffer
    x_min = data_min
    x_max = data_max
    #x_min = -0.01
    #x_max = 1.0
    
    # Create axis with formatting
    ax = Axis(fig[1, 1],
    #title=title,
    xlabel=xlabel,
    ylabel=ylabel,
    xlabelsize=label_size,
    ylabelsize=label_size,
    titlesize=label_size,
    xticklabelsize=text_size,
    yticklabelsize=text_size,
    spinewidth=1.5,
    xgridvisible=grid_visible,
    ygridvisible=grid_visible,
    xgridcolor=grid_color,
    ygridcolor=grid_color,
    xgridwidth=grid_linewidth,
    ygridwidth=grid_linewidth
    )
    
    # Set spine and tick colors
    
    ax.bottomspinecolor = axis_color
    ax.leftspinecolor = axis_color
    ax.rightspinecolor = axis_color
    ax.topspinecolor = axis_color
    
    #ax.xticks = 0:0.2:1.0
    #ax.yticks = 0:10:50
    
    ax.xticklabelcolor = axis_color
    ax.yticklabelcolor = axis_color
    ax.xlabelcolor = axis_color
    ax.ylabelcolor = axis_color
    ax.titlecolor = axis_color
    
    # Add extra padding to avoid cutting off labels
    #fig.margin = 20
    
    # Create the histogram
    hist = hist!(ax, data1,
    bins=bins,
    color=bar_color1,
    #strokecolor=edge_color,
    strokecolor=bar_color1,
    strokewidth=edge_linewidth,
    normalization=normalize ? :pdf : :none)
    
    hist = hist!(ax, data2,
    bins=bins,
    color=bar_color2,
    #strokecolor=edge_color,
    strokecolor=bar_color2,
    strokewidth=edge_linewidth,
    normalization=normalize ? :pdf : :none)
    
    # Optionally add KDE curve
    #if kde
    #density = kde!(ax, data,
    #color=kde_color,
    #linewidth=kde_linewidth,
    #label="KDE")
    #
    ## Add legend if KDE is used
    #axislegend(ax, position=:rt, framevisible=true,
    #framecolor=(:black, 0.2),
    #padding=(10, 10, 10, 10),
    #labelsize=text_size-2)
    #end
    
    # Adjust x limits
    ax.limits = (x_min, x_max, nothing, nothing)
    
    return fig
end


##### Another Parity Plot 
function alt_parity_plot(qhat_uq, res;
                        title="Parity Plot Subset",
                        xlabel="Heuristic Uncertainty (eV)",
                        ylabel="Residuals(eV)",
                        width=600,
                        colormap=:viridis,
                        color_value=0.6,
                        marker_size=10,
                        line_width=3.0,
                        axis_color=:black,
                        text_size=18,
                        label_size=22,
                        grid_visible=false,
                        grid_color=(:gray, 0.3),
                        grid_linewidth=0.5,
                        errorbar_color=nothing,
                        marker_color=nothing,
                        diagonal_color=:red,
                        diagonal_alpha=0.6,
                        diagonal_style=:dash)



    # Create figure and axis with better formatting
    fig = Figure(resolution=(width, width), fontsize=text_size, figure_padding=30)
    
    # Calculate min and max for setting plot limits
    min_val = min(minimum(qhat_uq), minimum(res))
    max_val = max(maximum(qhat_uq), maximum(res))
    #
    #min_val = -3257.35
    #max_val = -3257.15
    
    ## Add a small buffer to the limits to avoid cutting off points or labels
    #buffer = (max_val - min_val) * 0.05
    #plot_min = min_val - buffer
    #plot_max = max_val + buffer
    
    ax = Axis(fig[1, 1],
    #title=title,
    xlabel=xlabel,
    ylabel=ylabel,
    #limits=(plot_min, plot_max, plot_min, plot_max),
    #limits=(min_val, max_val, min_val-0.5, max_val+0.5),
    #limits=(0.0,1.0,0.0,1.0),
    titlesize=label_size,
    xlabelsize=label_size,
    ylabelsize=label_size,
    xticklabelsize=text_size,
    yticklabelsize=text_size,
    spinewidth=1.5,
    xgridvisible=grid_visible,
    ygridvisible=grid_visible,
    xgridcolor=grid_color,
    ygridcolor=grid_color,
    xgridwidth=grid_linewidth,
    ygridwidth=grid_linewidth
    )
    
    # Set spine and tick colors
    ax.bottomspinecolor = axis_color
    ax.leftspinecolor = axis_color
    ax.rightspinecolor = axis_color
    ax.topspinecolor = axis_color
    
    ax.xticklabelcolor = axis_color
    ax.yticklabelcolor = axis_color
    ax.xlabelcolor = axis_color
    ax.ylabelcolor = axis_color
    ax.titlecolor = axis_color
    
    # Add diagonal reference line
    lines!(ax, [min_val,max_val], [min_val,max_val],
    color=diagonal_color,
    linestyle=diagonal_style,
    linewidth=line_width - 1,  # Slightly thinner than main points
    alpha=diagonal_alpha,
    label="Perfect Prediction")
    
    
    # Scatter plot of points
    covered_qhatuq = Float64[]
    covered_res = Float64[]
    uncovered_qhatuq = Float64[]
    uncovered_res = Float64[]
    for i in eachindex(qhat_uq)
        if res[i] > qhat_uq[i]
            push!(uncovered_qhatuq, qhat_uq[i])
            push!(uncovered_res, res[i])
        else
            push!(covered_qhatuq, qhat_uq[i])
            push!(covered_res, res[i])
        end
    end
    
    #marker_color=:black
    #scatter!(ax,qhat_uq, res,
    #color=marker_color,
    #markersize=marker_size)
    
    @show uncovered_qhatuq
    @show uncovered_res
    scatter!(ax,covered_qhatuq, covered_res,
    color=:black,
    markersize=marker_size)
    
    
    scatter!(ax,uncovered_qhatuq, uncovered_res,
    color=:red,
    markersize=marker_size)
    # Set equal aspect ratio (usually important for parity plots)
    #ax.aspect = DataAspect()
    
    # Add legend with better formatting
    #axislegend(ax, position=:lt, framevisible=true, framecolor=(:black, 0.2),
    #padding=(10, 10, 10, 10), labelsize=text_size-2)
    
    return fig
end

##### Another Parity Plot 
function basic_parity_plot(ref, pred;
                        title="Parity Plot Subset",
                        xlabel="Reference Heuristic Uncertainty (eV)",
                        ylabel="Predicted Heuristic Uncertainty",
                        width=600,
                        colormap=:viridis,
                        color_value=0.6,
                        marker_size=10,
                        line_width=3.0,
                        text_size=18,
                        label_size=22,
                        grid_visible=false,
                        grid_color=(:gray, 0.3),
                        grid_linewidth=0.5,
                        errorbar_color=nothing,
                        marker_color=nothing,
                        diagonal_color=:red,
                        diagonal_alpha=0.6,
                        diagonal_style=:dash,
                        min_val = nothing,
                        max_val =nothing)


    # Create figure and axis with better formatting
    fig = Figure(resolution=(width, width), fontsize=text_size, figure_padding=30)
    
    # Calculate min and max for setting plot limits
    if isnothing(min_val) || isnothing(max_val)
        min_val = min(minimum(ref), minimum(pred))
        max_val = max(maximum(ref), maximum(pred))
    end
    #
    #min_val = -3257.35
    #max_val = -3257.15
    
    ## Add a small buffer to the limits to avoid cutting off points or labels
    buffer = (max_val - min_val) * 0.05
    plot_min = min_val - buffer
    plot_max = max_val + buffer
    
    ax = Axis(fig[1, 1],
             #title=title,
             xlabel=xlabel,
             ylabel=ylabel,
             limits=(plot_min, plot_max, plot_min, plot_max),
             #limits=(min_val, max_val, min_val-0.5, max_val+0.5),
             #limits=(0.0,1.0,0.0,1.0),
             titlesize=label_size,
             xlabelsize=label_size,
             ylabelsize=label_size,
             xticklabelsize=text_size,
             yticklabelsize=text_size,
             spinewidth=1.5,
             xgridvisible=grid_visible,
             ygridvisible=grid_visible,
             xgridcolor=grid_color,
             ygridcolor=grid_color,
             xgridwidth=grid_linewidth,
             ygridwidth=grid_linewidth
             )
    
    
    # Add diagonal reference line
    lines!(ax, [min_val,max_val], [min_val,max_val],
           color=diagonal_color,
           linestyle=diagonal_style,
           linewidth=line_width - 1,  # Slightly thinner than main points
           alpha=diagonal_alpha,
           label="Perfect Prediction")
    
    
    scatter!(ref,pred,
            color=:black,
            markersize=marker_size)
    
    
    # Set equal aspect ratio (usually important for parity plots)
    ax.aspect = DataAspect()
    
   
    return fig
end
