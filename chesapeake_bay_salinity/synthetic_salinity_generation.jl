using LinearAlgebra, Random, DataFrames, Dates

"""
    simulate_large_grid(grid_x, grid_y, grid_t, params; mean_val=0.0, nugget=1e-6)

Generates high-dimensional Gaussian Random Fields on a regular grid using Kronecker factorization.

# Arguments
- `grid_x`: Vector of X-coordinates (e.g., longitudes/easting)
- `grid_y`: Vector of Y-coordinates (e.g., latitudes/northing)
- `grid_t`: Vector of timestamps (e.g., hours or Float64 steps)
- `params`: Vector/Tuple `[σ_sq, λ_space, λ_time]`
"""
function simulate_large_grid(
    grid_x::AbstractVector{Float64},
    grid_y::AbstractVector{Float64},
    grid_t::AbstractVector{Float64},
    params::Vector{Float64};
    mean_val::Float64 = 0.0,
    nugget::Float64 = 1e-6,
    rng::AbstractRNG = Random.GLOBAL_RNG
)
    σ_sq, λ_s, λ_t = params
    
    Nx = length(grid_x)
    Ny = length(grid_y)
    Ns = Nx * Ny
    Nt = length(grid_t)

    # 1. Generate 2D Spatial Grid Coordinates
    spatial_coords = Vector{Tuple{Float64, Float64}}(undef, Ns)
    idx = 1
    for x in grid_x, y in grid_y
        spatial_coords[idx] = (x, y)
        idx += 1
    end

    # 2. Compute Spatial Covariance Matrix Σ_s (Ns x Ns)
    Σ_s = Matrix{Float64}(undef, Ns, Ns)
    for i in 1:Ns
        p1 = spatial_coords[i]
        for j in i:Ns
            p2 = spatial_coords[j]
            h = hypot(p1[1] - p2[1], p1[2] - p2[2]) # Spatial distance
            val = exp(-h / λ_s)                    # Spatial correlation
            Σ_s[i, j] = val
            Σ_s[j, i] = val
        end
    end

    # 3. Compute Temporal Covariance Matrix Σ_t (Nt x Nt)
    Σ_t = Matrix{Float64}(undef, Nt, Nt)
    for i in 1:Nt
        for j in i:Nt
            u = abs(grid_t[i] - grid_t[j])        # Temporal distance
            val = exp(-u / λ_t)                   # Temporal correlation
            Σ_t[i, j] = val
            Σ_t[j, i] = val
        end
    end

    # 4. Diagonal Ridge Regularization (Nugget effect)
    for i in 1:Ns; Σ_s[i, i] += nugget; end
    for i in 1:Nt; Σ_t[i, i] += nugget; end

    # 5. Factorize small spatial & temporal matrices independently
    L_s = cholesky(Symmetric(Σ_s)).L
    L_t = cholesky(Symmetric(Σ_t)).L

    # 6. Kronecker Sampling Matrix Trick: Y = μ + sqrt(σ²) * (L_s * Z * L_t')
    Z = randn(rng, Ns, Nt)
    Y = mean_val .+ sqrt(σ_sq) .* (L_s * Z * L_t')

    return spatial_coords, grid_t, Y
end

"""
    grid_to_dataframe(grid_x, grid_y, grid_t, Y)

Flattens the 3D grid matrix output into a standard DataFrame for downstream tasks.
"""
function grid_to_dataframe(grid_x, grid_y, grid_t, Y)
    Nx, Ny, Nt = length(grid_x), length(grid_y), length(grid_t)
    total_rows = Nx * Ny * Nt

    df_x = Vector{Float64}(undef, total_rows)
    df_y = Vector{Float64}(undef, total_rows)
    df_t = Vector{Float64}(undef, total_rows)
    df_val = Vector{Float64}(undef, total_rows)

    idx = 1
    s_idx = 1
    for x in grid_x
        for y in grid_y
            for t_idx in 1:Nt
                df_x[idx] = x
                df_y[idx] = y
                df_t[idx] = grid_t[t_idx]
                df_val[idx] = Y[s_idx, t_idx]
                idx += 1
            end
            s_idx += 1
        end
    end

    return DataFrame(x = df_x, y = df_y, time = df_t, y_val = df_val)
end

"""
    gaussian_plume_mean(x, y, t, x0, y0, t0, Q, vx, vy, D, S_ambient)

Computes the 2D advection-diffusion analytical mean plume concentration from a continuous source.
"""
function gaussian_plume_mean(x, y, t, x0, y0, t0, Q, vx, vy, D, S_ambient)
    dt = t - t0
    if dt <= 0.0
        return S_ambient
    end
    
    # Distance from moving advective center
    cx = x0 + vx * dt
    cy = y0 + vy * dt
    r_sq = (x - cx)^2 + (y - cy)^2
    
    # 2D Advection-Diffusion analytical solution
    plume_conc = (Q / (4 * π * D * dt)) * exp(-r_sq / (4 * D * dt))
    return S_ambient + plume_conc
end

"""
    simulate_salinity_plume(grid_x, grid_y, grid_t, params, plume_params)

Generates anisotropic salinity plume data using a spatio-temporal Gaussian Random Field superimposed 
on an advection-diffusion mean trend.
"""
function simulate_salinity_plume(
    grid_x::AbstractVector{Float64},
    grid_y::AbstractVector{Float64},
    grid_t::AbstractVector{Float64},
    params::Vector{Float64},        # [σ_sq, λ_along, λ_trans, λ_time]
    plume_params::NamedTuple;      # (source_x, source_y, t0, release_rate, vx, vy, diffusion, ambient_salinity)
    nugget::Float64 = 1e-6,
    rng::AbstractRNG = Random.GLOBAL_RNG
)
    σ_sq, λ_along, λ_trans, λ_t = params
    
    Nx, Ny, Nt = length(grid_x), length(grid_y), length(grid_t)
    Ns = Nx * Ny

    # 1. Flow Angle for Anisotropic Rotation
    θ = atan(plume_params.vy, plume_params.vx)
    cos_θ, sin_θ = cos(θ), sin(θ)

    # 2. Build Anisotropic Spatial Covariance Matrix
    spatial_coords = [(x, y) for x in grid_x for y in grid_y]
    Σ_s = Matrix{Float64}(undef, Ns, Ns)

    for i in 1:Ns
        x1, y1 = spatial_coords[i]
        for j in i:Ns
            x2, y2 = spatial_coords[j]
            
            dx = x2 - x1
            dy = y2 - y1
            
            # Rotate coordinates into principal flow direction
            d_along =  dx * cos_θ + dy * sin_θ
            d_trans = -dx * sin_θ + dy * cos_θ
            
            # Anisotropic distance metric
            h_aniso = sqrt((d_along / λ_along)^2 + (d_trans / λ_trans)^2)
            
            val = exp(-h_aniso)
            Σ_s[i, j] = val
            Σ_s[j, i] = val
        end
    end

    # 3. Build Temporal Covariance Matrix
    Σ_t = Matrix{Float64}(undef, Nt, Nt)
    for i in 1:Nt, j in i:Nt
        u = abs(grid_t[i] - grid_t[j])
        val = exp(-u / λ_t)
        Σ_t[i, j] = val
        Σ_t[j, i] = val
    end

    # 4. Cholesky Factorization
    for i in 1:Ns; Σ_s[i, i] += nugget; end
    for i in 1:Nt; Σ_t[i, i] += nugget; end

    L_s = cholesky(Symmetric(Σ_s)).L
    L_t = cholesky(Symmetric(Σ_t)).L

    # 5. Generate Stochastic Residual Field Z ~ N(0, Σ)
    Z = sqrt(σ_sq) .* (L_s * randn(rng, Ns, Nt) * L_t')

    # 6. Combine Deterministic Mean Trend + Stochastic Residual Field
    Salinity = Matrix{Float64}(undef, Ns, Nt)
    
    for (s_idx, (x, y)) in enumerate(spatial_coords)
        for (t_idx, t) in enumerate(grid_t)
            
            # Evaluate deterministic plume background
            μ = gaussian_plume_mean(
                x, y, t, 
                plume_params.source_x, plume_params.source_y, plume_params.t0,
                plume_params.release_rate, plume_params.vx, plume_params.vy, 
                plume_params.diffusion, plume_params.ambient_salinity
            )

            # Add stochastic turbulence and enforce physical lower bound
            val = μ + Z[s_idx, t_idx]
            Salinity[s_idx, t_idx] = max(plume_params.ambient_salinity, val)
        end
    end

    return spatial_coords, grid_t, Salinity
end

"""
    gaussian_plume_continuous_mean(x, y, t, plume_params)

Calculates continuous plume accumulation from a continuous release source at point (x0, y0).
"""
function gaussian_plume_continuous_mean(x, y, t, p)
    if t <= p.t0
        return p.ambient_salinity
    end
    
    # Distance from source along current vector
    dx = x - p.source_x
    dy = y - p.source_y
    
    # Rotate into current direction (vx, vy)
    U = hypot(p.vx, p.vy)
    if U == 0.0
        return p.ambient_salinity
    end
    
    x_along = (dx * p.vx + dy * p.vy) / U
    y_trans = (-dx * p.vy + dy * p.vx) / U
    
    # Steady-state / Continuous release advection-diffusion approximation
    r = hypot(x_along, y_trans)
    if r == 0.0
        return p.ambient_salinity + (p.release_rate / (2 * π * p.diffusion))
    end
    
    # Continuous plume decay profile
    conc = (p.release_rate / (2 * π * p.diffusion)) * exp((x_along * U - r * U) / (2 * p.diffusion))
    return p.ambient_salinity + max(0.0, conc)
end

"""
    simulate_longterm_plume(grid_x, grid_y, t_start, t_end, dt, cov_params, plume_params; output_file=nothing)

Sequentially generates synthetic salinity plume frames over long time steps using O(1) memory per step.
"""
function simulate_longterm_plume(
    grid_x::AbstractVector{Float64},
    grid_y::AbstractVector{Float64},
    t_start::Float64,
    t_end::Float64,
    dt::Float64,                      # Time step size (e.g. 0.5 hours)
    cov_params::Vector{Float64},      # [σ_sq, λ_along, λ_trans, λ_time]
    plume_params::NamedTuple;
    nugget::Float64 = 1e-6,
    rng::AbstractRNG = Random.GLOBAL_RNG
)
    σ_sq, λ_along, λ_trans, λ_t = cov_params
    Nx, Ny = length(grid_x), length(grid_y)
    Ns = Nx * Ny

    # 1. Flow Angle for Spatial Anisotropy
    θ = atan(plume_params.vy, plume_params.vx)
    cos_θ, sin_θ = cos(θ), sin(θ)

    # 2. Build and Factorize Spatial Covariance Matrix ONCE
    spatial_coords = [(x, y) for x in grid_x for y in grid_y]
    Σ_s = Matrix{Float64}(undef, Ns, Ns)

    for i in 1:Ns
        x1, y1 = spatial_coords[i]
        for j in i:Ns
            x2, y2 = spatial_coords[j]
            dx, dy = x2 - x1, y2 - y1
            
            d_along =  dx * cos_θ + dy * sin_θ
            d_trans = -dx * sin_θ + dy * cos_θ
            
            h_aniso = sqrt((d_along / λ_along)^2 + (d_trans / λ_trans)^2)
            val = exp(-h_aniso)
            Σ_s[i, j] = val
            Σ_s[j, i] = val
        end
    end
    
    for i in 1:Ns; Σ_s[i, i] += nugget; end
    L_s = cholesky(Symmetric(Σ_s)).L

    # 3. Setup Time Stepping
    times = collect(t_start:dt:t_end)
    Nt = length(times)
    
    # Pre-calculate AR(1) recurrence constant for constant dt
    α = exp(-dt / λ_t)
    β = sqrt(1.0 - α^2)

    # 4. Initialize Residual Field Z_0 ~ N(0, σ² Σ_s)
    Z_current = sqrt(σ_sq) .* (L_s * randn(rng, Ns))
    
    # 3D Array to hold salinity (Ns x Nt) - or write directly to disk/stream
    Salinity_3D = Array{Float64, 3}(undef, Nx, Ny, Nt)

    # 5. Continuous Sequential Time Loop
    for (t_idx, t) in enumerate(times)
        if t_idx > 1
            # Exact Continuous AR(1) State Update
            w = randn(rng, Ns)
            Z_current = α .* Z_current .+ (β * sqrt(σ_sq)) .* (L_s * w)
        end

        # Combine deterministic mean trend + stochastic state
        s_idx = 1
        for (i, x) in enumerate(grid_x)
            for (j, y) in enumerate(grid_y)
                μ = gaussian_plume_continuous_mean(x, y, t, plume_params)
                val = μ + Z_current[s_idx]
                
                Salinity_3D[i, j, t_idx] = max(plume_params.ambient_salinity, val)
                s_idx += 1
            end
        end
    end

    return spatial_coords, times, Salinity_3D
end