using NCDatasets
using CairoMakie
using Dates

function animate_hatteras_currents(output_gif::String)
    # Correct NOAA ERDDAP dataset ID for East Coast HF Radar Surface Currents (6km)
    erddap_url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdHFR364km"
    
    println("Connecting to NOAA HF Radar stream ($erddap_url)...")
    ds = NCDataset(erddap_url)

    # 1. Fetch dimension coordinate vectors
    lons_all = ds["longitude"][:]
    lats_all = ds["latitude"][:]
    times_all = ds["time"][:]

    # Filter bounding box around North Carolina / Cape Hatteras
    # Longitude: -77.5°W to -71.5°W, Latitude: 33.0°N to 37.5°N
    lon_indices = findall(x -> -77.5 <= x <= -71.5, lons_all)
    lat_indices = findall(x -> 33.0 <= x <= 37.5, lats_all)

    lons = lons_all[lon_indices]
    lats = lats_all[lat_indices]

    # Select the most recent 48 hourly time steps
    num_frames = 48
    total_times = length(times_all)
    time_indices = (total_times - num_frames + 1):total_times
    times = times_all[time_indices]

    println("Fetching surface current grids ($num_frames hourly frames)...")

    # Load u (eastward) and v (northward) surface velocity components (m/s)
    # Array dimensions in ERDDAP griddap: (longitude, latitude, depth, time)
    u_raw = ds["u"][lon_indices, lat_indices, 1, time_indices]
    v_raw = ds["v"][lon_indices, lat_indices, 1, time_indices]

    # Clean missing values/NaNs into standard Float64 matrices
    u_data = Float64.(coalesce.(u_raw, NaN))
    v_data = Float64.(coalesce.(v_raw, NaN))

    # 2. Setup CairoMakie Figure & Dark Theme
    fig = Figure(size = (900, 750), backgroundcolor = :black)
    
    ax = Axis(
        fig[1, 1],
        title = "Cape Hatteras Surface Currents (US CoastWatch HF Radar)",
        titlesize = 17,
        titlecolor = :white,
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)",
        xticklabelcolor = :white,
        yticklabelcolor = :white,
        xlabelcolor = :white,
        ylabelcolor = :white,
        backgroundcolor = "#080c16",
        limits = ((-77.5, -71.5), (33.0, 37.5))
    )

    # Reactive observables for frame playback
    frame_idx = Observable(1)
    
    speed_matrix = @lift begin
        u_slice = u_data[:, :, 1, $frame_idx]
        v_slice = v_data[:, :, 1, $frame_idx]
        
        # Velocity magnitude: sqrt(u^2 + v^2)
        spd = sqrt.(u_slice.^2 .+ v_slice.^2)
        return spd
    end

    time_label = @lift begin
        t = times[$frame_idx]
        try
            return "Time (UTC): " * Dates.format(DateTime(t), "yyyy-mm-dd HH:MM")
        catch
            return "Time (UTC): " * string(t)
        end
    end

    # Plot surface current speed heatmap
    hm = heatmap!(
        ax, 
        lons, 
        lats, 
        speed_matrix, 
        colormap = :inferno, 
        colorrange = (0.0, 2.0), # Currents range up to ~2 m/s (~4 knots)
        nan_color = :black
    )

    Colorbar(
        fig[1, 2], 
        hm, 
        label = "Surface Current Speed (m s⁻¹)", 
        labelcolor = :white,
        ticklabelcolor = :white
    )

    text!(
        ax, 
        time_label, 
        position = (-77.0, 37.0), 
        color = :white, 
        fontsize = 14
    )

    # 3. Render and save GIF
    println("Rendering animation to $output_gif...")
    record(fig, output_gif, 1:num_frames; framerate = 10) do i
        frame_idx[] = i
    end

    close(ds)
    println("Done! Successfully saved to $output_gif")
end

# Run function
animate_hatteras_currents("cape_hatteras_gulf_stream.gif")