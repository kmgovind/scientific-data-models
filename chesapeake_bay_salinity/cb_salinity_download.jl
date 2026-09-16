using NCDatasets, Dates, Downloads

# Define the date range
num_days = 2;
start_date = Date(2026, 6, 18);
dates = [start_date + Day(i) for i in 0:(num_days-1)];
cycles = ["00", "06", "12", "18"]; # Define the cycles to download

# Print a warning about the amount of data to be downloaded
println("Attempting to download all available CBOFS files for $num_days days...")
println("WARNING: This may download roughly $(num_days * 4) GB of data.")

# Define counters for successful and failed downloads
total_ok = 0
total_fail = 0

# Download the data for each date and cycle
for date in dates

    # Define the datestring
    yyyy = string(year(date))
    mm = lpad(string(month(date)), 2, "0")
    dd = lpad(string(day(date)), 2, "0")
    yyyymmdd = yyyy * mm * dd

    # Define the download directory for datafiles
    destination_dir = "datafiles/"

    for cycle in cycles
        # 2. Get the daily THREDDS catalog for this specific date
        catalog_url = "https://opendap.co-ops.nos.noaa.gov/thredds/catalog/NOAA/CBOFS/MODELS/$yyyy/$mm/$dd/catalog.html"

        f_hours = String[]
        try
            # Download catalog HTML into memory
            catalog_tmp = IOBuffer()
            Downloads.download(catalog_url, catalog_tmp)
            html = String(take!(catalog_tmp))

            # Match: cbofs.t18z.20260216.fields.n001.nc
            pat = Regex("cbofs\\.t$(cycle)z\\.$yyyymmdd\\.fields\\.n(\\d{3})\\.nc")
            f_hours = unique([m.captures[1] for m in eachmatch(pat, html)])
            sort!(f_hours)

            if isempty(f_hours)
                println("No files found in catalog for $yyyymmdd cycle t$(cycle)z.")
                continue
            end

            println("\nFound $(length(f_hours)) files for $yyyymmdd (t$(cycle)z): ", join(f_hours, ", "))
        catch e
            println("Could not read catalog for $yyyymmdd. Error: ", e)
            continue
        end

        # 3. Download every available f_hour for that cycle
        for f_hour in f_hours
            thredds_url = "https://opendap.co-ops.nos.noaa.gov/thredds/fileServer/NOAA/CBOFS/MODELS/$yyyy/$mm/$dd/cbofs.t$(cycle)z.$yyyymmdd.fields.n$f_hour.nc"
            local_file = "$destination_dir/chesapeake_salinity_$(yyyymmdd)_t$(cycle)z_n$(f_hour).nc"

            # Skip if we already downloaded it (useful if script crashes and you restart)
            if isfile(local_file)
                println("  Skipping $local_file (Already exists)")
                continue
            end

            try
                print("  Downloading $yyyymmdd t$(cycle)z n$f_hour... ")
                Downloads.download(thredds_url, local_file)
                println("Success!")
                global total_ok += 1
            catch e
                println("Failed. Error: ", e)
                global total_fail += 1
            end
        end
    end
end

println("\nDone. Successful new downloads: $total_ok | Failed: $total_fail")