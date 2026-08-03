module Variograms

using LinearAlgebra, LsqFit, StatsBase, DataFrames, Distances, Dates

# Function to calculate pairwise distances
function pairwise_distances(lon, lat, y)
    n = length(y)
    h = zeros(eltype(lon), n, n)
    γ = zeros(eltype(y), n, n)

    for i in 1:n
        # Pre-extract tuples to avoid heap allocations in the inner loop
        p1 = (lon[i], lat[i])
        y_i = y[i]
        
        # Only compute the upper triangle (j > i)
        for j in (i+1):n
            p2 = (lon[j], lat[j])
            
            # If using Distances.haversine(p1, p2, r), it defaults to meters. 
            # If your custom haversine already returns meters, keep the /1000.
            d = haversine(p1, p2) / 1000 
            g = 0.5 * (y_i - y[j])^2
            
            # Symmetric assignment
            h[i, j] = d
            h[j, i] = d
            
            γ[i, j] = g
            γ[j, i] = g
        end
    end

    return h, γ
end

function empirical_spatial_variogram(lon, lat, y; n_bins=50)
    n = length(y)

    # -----------------------------------------------------------
    # First pass: determine maximum pairwise distance
    # -----------------------------------------------------------
    h_max = zero(eltype(lon))

    for i in 1:n-1
        p1 = (lon[i], lat[i])

        for j in i+1:n
            d = haversine(p1, (lon[j], lat[j])) / 1000

            if d > h_max
                h_max = d
            end
        end
    end

    bin_width = h_max / n_bins

    bin_sums = zeros(Float64, n_bins)
    bin_counts = zeros(Int, n_bins)

    # -----------------------------------------------------------
    # Second pass: compute semivariances and bin immediately
    # -----------------------------------------------------------
    for i in 1:n-1
        p1 = (lon[i], lat[i])
        yi = y[i]

        for j in i+1:n
            d = haversine(p1, (lon[j], lat[j])) / 1000

            γ = 0.5 * (yi - y[j])^2

            bin = d == h_max ? n_bins : floor(Int, d / bin_width) + 1

            @inbounds begin
                bin_sums[bin] += γ
                bin_counts[bin] += 1
            end
        end
    end

    # -----------------------------------------------------------
    # Construct DataFrame
    # -----------------------------------------------------------
    h_mid = Float64[]
    variogram = Float64[]
    count = Int[]

    for b in 1:n_bins
        if bin_counts[b] > 0
            push!(h_mid, (b - 0.5) * bin_width)
            push!(variogram, bin_sums[b] / bin_counts[b])
            push!(count, bin_counts[b])
        end
    end

    return DataFrame(
        h_mid = h_mid,
        variogram = variogram,
        count = count,
    )
end

function empirical_temporal_variogram(times, y; n_bins=50)
    n = length(y)

    # -----------------------------------------------------------
    # First pass: determine maximum time lag (hours)
    # -----------------------------------------------------------
    h_max = 0.0

    for i in 1:n-1
        for j in i+1:n
            d = abs(Dates.value(times[j] - times[i])) / (1000 * 60 * 60)  # hours

            if d > h_max
                h_max = d
            end
        end
    end

    bin_width = h_max / n_bins

    bin_sums = zeros(Float64, n_bins)
    bin_counts = zeros(Int, n_bins)

    # -----------------------------------------------------------
    # Second pass: compute semivariances and bin immediately
    # -----------------------------------------------------------
    for i in 1:n-1
        yi = y[i]

        for j in i+1:n
            d = abs(Dates.value(times[j] - times[i])) / (1000 * 60 * 60)  # hours

            γ = 0.5 * (yi - y[j])^2

            bin = d == h_max ? n_bins : floor(Int, d / bin_width) + 1

            @inbounds begin
                bin_sums[bin] += γ
                bin_counts[bin] += 1
            end
        end
    end

    # -----------------------------------------------------------
    # Construct DataFrame
    # -----------------------------------------------------------
    h_mid = Float64[]
    variogram = Float64[]
    count = Int[]

    for b in 1:n_bins
        if bin_counts[b] > 0
            push!(h_mid, (b - 0.5) * bin_width)
            push!(variogram, bin_sums[b] / bin_counts[b])
            push!(count, bin_counts[b])
        end
    end

    return DataFrame(
        h_mid = h_mid,          # time lag (hours)
        variogram = variogram,
        count = count,
    )
end

function covariance(x1, x2, λₓ, σ)
    return (σ^2) .* exp(-norm.(x1.-x2)/λₓ);
end

function spatio_variogram(h,p)
    return (p[1]^2) .* (1 .- exp.(-h./ p[2]))
end

function spatio_RBF_variogram(h,p)
    return (p[1]^2) .* (1 .- exp.(-(h.^2)./ (2 .*(p[2]).^2)))
end

function temporal_variogram(t,p)
    return p[2].*t .- p[3].*(cos.(2 .*pi.*t./12.5) .- 1)
end

function spatiotemporal_variogram(lambda,beta_0,beta_1,beta_2,h,u,l)
    return ((lambda^2) .* exp(-h./l)) .* (beta_0 - beta_1.*u + beta_2.*(cos.(pi.*u./12.5) - 1))
end

"Fit hyperparameters to spatial measurements"
function hp_fit(lon, lat, y)
    emp_vario = empirical_spatial_variogram(lon, lat, y);
    param_fit = curve_fit(spatio_variogram, emp_vario.h_mid, emp_vario.variogram, [1.0,1.0]);
    return param_fit.param[1], param_fit.param[2]
end

"Fit hyperparameters to temporal measurements"
function hp_fit(times, y)
    emp_vario = empirical_temporal_variogram(times, y);
    param_fit = curve_fit(temporal_variogram, emp_vario.h_mid, emp_vario.variogram, [1.0,1.0,1.0]);
    return param_fit.param[1], param_fit.param[2], param_fit.param[3]
end

"Fit hyperparameters to spatiotemporal measurements"
function hp_fit(measurements)
    h, u, γ = pairwise_distances_st(measurements);
    emp_vario = empirical_variogram_st(h, u, γ);
    lags = [emp_vario.h_mid emp_vario.u_mid]
    param_fit = curve_fit(spatiotemporal_variogram, lags, emp_vario.variogram, [1.0,1.0,1.0,1.0,1.0,1.0])
    return param_fit.param[1], param_fit.param[2], param_fit.param[3]
end

# Function to calculate pairwise distances
function pairwise_distances_st(measurements)
    n = length(measurements)
    h = zeros(n, n)  # Spatial distances
    u = zeros(n, n)  # Temporal distances
    γ = zeros(n, n)  # Semi-variogram

    for i in 1:n
        for j in 1:n
            h[i, j] = sqrt((measurements[i].p[1] - measurements[j].p[1])^2 + (measurements[i].p[2] - measurements[j].p[2])^2)
            u[i, j] = abs(measurements[i].t - measurements[j].t)
            γ[i, j] = 0.5 * (measurements[i].y - measurements[j].y)^2
        end
    end
    return h, u, γ
end

# Function to calculate the empirical variogram
function empirical_variogram_st(h, u, γ, n_bins=50)
    h_max = maximum(h)
    u_max = maximum(u)
    h_edges = range(0, h_max, length=n_bins+1)
    u_edges = range(0, u_max, length=n_bins+1)
    
    bins = DataFrame(h_mid = Float64[], u_mid = Float64[], variogram = Float64[], count = Int64[])
    
    for i in 1:n_bins
        for j in 1:n_bins
            h_bin = h_edges[i] .<= h .< h_edges[i+1]
            u_bin = u_edges[j] .<= u .< u_edges[j+1]
            mask = h_bin .& u_bin
            bin_count = count(mask)
            if bin_count > 0
                h_mid = (h_edges[i] + h_edges[i+1]) / 2
                u_mid = (u_edges[j] + u_edges[j+1]) / 2
                variogram_value = mean(γ[mask])
                push!(bins, (h_mid, u_mid, variogram_value, bin_count))
            end
        end
    end
    return bins
end

function matern12_variogram(h, u, σ_sq, λₓ, λₜ)
    return σ_sq .* (1 .- exp.(-h./λₓ) .* exp.(-u./λₜ))
end

function matern12_variogram(lags, p)
    h = lags[:,1]
    u = lags[:,2]
    return p[1] .* (1 .- exp.(-h./p[2]).* exp.(-u./p[3]))
end


end


