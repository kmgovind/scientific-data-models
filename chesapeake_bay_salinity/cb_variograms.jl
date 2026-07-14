module Variograms

using LinearAlgebra, LsqFit, StatsBase, DataFrames, Distances

# Function to calculate pairwise distances
function pairwise_distances(lon, lat, y)
    n = length(y)

    h = zeros(n, n)
    γ = zeros(n, n)

    for i in 1:n
        for j in 1:n
            h[i, j] = haversine([lon[i], lat[i]], [lon[j], lat[j]])/1000
            γ[i, j] = 0.5 * (y[i] - y[j])^2
        end
    end

    return h, γ
end

# Function to calculate the empirical variogram
function empirical_variogram(h, γ, n_bins=50)
    # h_max = maximum(h)
    h_max = maximum(h)
    h_edges = range(0, h_max, length=n_bins+1)
    
    bins = DataFrame(h_mid = Float64[], variogram = Float64[], count = Int64[])
    
    for i in 1:n_bins
        h_bin = h_edges[i] .<= h .< h_edges[i+1]
        bin_count = count(h_bin)
        if bin_count > 0
            h_mid = (h_edges[i] + h_edges[i+1]) / 2
            variogram_value = mean(γ[h_bin])
            push!(bins, (h_mid, variogram_value, bin_count))
            # push the data all together at the end rather than on each loop
        end
    end

    return bins
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


