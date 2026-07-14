module Variograms

using LinearAlgebra, LsqFit, StatsBase, DataFrames, Distances

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

# Function to calculate the empirical variogram
function empirical_variogram(h, γ, n_bins=50)
    h_max = maximum(h)
    h_edges = range(0, h_max, length=n_bins+1)
    
    # Pre-allocate flat vectors for performance
    bin_sums = zeros(eltype(γ), n_bins)
    bin_counts = zeros(Int64, n_bins)
    
    # Calculate the width of each bin for fast lookup mapping
    bin_width = h_max / n_bins
    
    n = size(h, 1)
    
    # Accumulate sums and counts in a single pass over the upper triangle
    for i in 1:n
        for j in (i+1):n
            h_val = h[i, j]
            
            # Fast mathematical indexing to find which bin h_val belongs to
            if h_val >= h_max
                bin_idx = n_bins
            else
                bin_idx = floor(Int, h_val / bin_width) + 1
            end
            
            if 1 <= bin_idx <= n_bins
                bin_sums[bin_idx] += γ[i, j]
                bin_counts[bin_idx] += 1
            end
        end
    end
    
    # Construct the final vectors for rows that actually have data
    h_mids = Float64[]
    variograms = Float64[]
    counts = Int64[]
    
    # Filter out empty bins and compute means
    for i in 1:n_bins
        if bin_counts[i] > 0
            h_mid = (h_edges[i] + h_edges[i+1]) / 2
            push!(h_mids, h_mid)
            push!(variograms, bin_sums[i] / bin_counts[i])
            push!(counts, bin_counts[i])
        end
    end
    
    # Construct the DataFrame all at once (fast!)
    return DataFrame(h_mid = h_mids, variogram = variograms, count = counts)
end

function covariance(x1, x2, λₓ, σ)
    return (σ^2) .* exp(-norm.(x1.-x2)/λₓ);
end

function matern12_variogram(h, σ_sq, λₓ)
    return σ_sq .* (1 .- exp.(-h./λₓ))
end

function matern12_variogram(h, p)
    return p[1] .* (1 .- exp.(-h./p[2]))
end

function matern12_log(h, σ_sq, λₓ)
    return log(σ_sq).-(h./λₓ)
    # return -(1/λₓ).*h
end

function matern12_lin(h, σ_sq, λₓ)
    return -1/λₓ
end

function spatio_variogram(h,p)
    return (p[1]^2) .* (1 .- exp.(-h./p[2]))
end

function temporal_variogram(beta_0,beta_1,beta_2,u)
    return beta_0 - beta_1.*u + beta_2.*(cos.(pi.*u./12.5) - 1)
end

function spatiotemporal_variogram(lambda,beta_0,beta_1,beta_2,h,u,l)
    return ((lambda^2) .* exp(-h./l)) .* (beta_0 - beta_1.*u + beta_2.*(cos.(pi.*u./12.5) - 1))
end

"Fit hyperparameters to measurements"
function hp_fit(lon, lat, y)
    h, γ = pairwise_distances(lon, lat, y);
    emp_vario = empirical_variogram(h, γ);
    param_fit = curve_fit(spatio_variogram, emp_vario.h_mid, emp_vario.variogram, [1.0,1.0]);
    return param_fit.param[1], param_fit.param[2]
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


