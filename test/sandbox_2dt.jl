using CartesianGeometry
using Test
const T=Float64
using StaticArrays

x = collect(0.0:0.05:1.0)
y = collect(0.0:0.05:1.0)
z = [0.0, 1.0]
xyz = (x, y, z)

# define level set
# Define initial and final circle parameters
const R_init = 0.1    # Initial radius
const R_final = 0.2   # Final radius 
const a, b = 0.5, 0.5 # Circle center

levelset1 = (x, y, z) -> (x-a)^2 + (y-b)^2 - R_init^2
levelset2 = (x, y, z) -> (x-a)^2 + (y-b)^2 - R_final^2

tn,tn1 = 0.0, 1.0
dt = tn1 - tn
levelset = (x, y, t) -> (t-tn)/dt * levelset2(x, y, 0) + (1 - (t-tn)/dt) * levelset1(x, y, 0)

V, bary, interface_length, cell_types = integrate(Tuple{0}, levelset, xyz, T, nan)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)

using CSV
using DataFrames

"""
    write_array_to_csv(array, filename; delimiter=',', header=nothing, index=nothing)

Write a Julia array to a CSV file.
Handles 1D, 2D and 3D arrays. For 3D arrays, it creates separate files for each slice.

# Arguments
- `array`: Array to write to the file
- `filename`: Name of the output file. If the path doesn't exist, it will be created.
- `delimiter`: Delimiter to use in the CSV file. Default: ','
- `header`: List of column names. If nothing, no header will be written.
- `index`: List of row names. If nothing, no indices will be written.

# Returns
- The absolute path of the created file, or a vector of paths for 3D arrays.
"""
function write_array_to_csv(array, filename; delimiter=',', header=nothing, index=nothing)
    # Ensure directory exists
    mkpath(dirname(abspath(filename)))
    
    # Handle 3D arrays
    if ndims(array) == 3
        filenames = String[]
        basefilename = splitext(filename)
        
        for k in 1:size(array, 3)
            slice_filename = "$(basefilename[1])_slice$(k)$(basefilename[2])"
            slice_path = write_array_to_csv(array[:,:,k], slice_filename, 
                                          delimiter=delimiter, header=header, index=index)
            push!(filenames, slice_path)
        end
        
        println("3D Array written to $(length(filenames)) files")
        return filenames
    end
    
    # Convert to DataFrame if header or index is provided
    if !isnothing(header) || !isnothing(index)
        df = DataFrame(array, :auto)
        
        if !isnothing(header)
            if length(header) != size(array, 2)
                error("Header length ($(length(header))) must match the number of columns ($(size(array, 2)))")
            end
            rename!(df, Symbol.(header))
        end
        
        if !isnothing(index)
            if length(index) != size(array, 1)
                error("Index length ($(length(index))) must match the number of rows ($(size(array, 1)))")
            end
            insertcols!(df, 1, :Index => index)
        end
        
        CSV.write(filename, df, delim=delimiter)
    else
        # Simple write without headers
        open(filename, "w") do io
            for i in 1:size(array, 1)
                println(io, join(array[i, :], delimiter))
            end
        end
    end
    
    println("Array written to file: $filename")
    return abspath(filename)
end

"""
    write_dict_of_arrays_to_csv(data_dict, base_filename; delimiter=',')

Write a dictionary of Julia arrays to multiple CSV files.
Handles 1D, 2D and 3D arrays.

# Arguments
- `data_dict`: Dictionary with keys as names and values as arrays
- `base_filename`: Base name for the output files (without extension).
  Files will be named {base_filename}_{key}.csv
- `delimiter`: Delimiter to use in the CSV file. Default: ','

# Returns
- List of absolute paths of the created files.
"""
function write_dict_of_arrays_to_csv(data_dict, base_filename; delimiter=',')
    files_created = String[]
    
    # Ensure directory exists
    output_dir = dirname(abspath(base_filename))
    mkpath(output_dir)
    
    # Base name without extension or path
    base_name = basename(base_filename)
    if contains(base_name, ".")
        base_name = split(base_name, ".")[1]
    end
    
    # Write each array to a separate file
    for (key, array) in data_dict
        filename = joinpath(output_dir, "$(base_name)_$(key).csv")
        
        # Handle different array dimensions
        if ndims(array) == 1
            # 1D array: reshape to column vector for CSV writing
            reshaped_array = reshape(array, :, 1)
            path = write_array_to_csv(reshaped_array, filename, delimiter=delimiter)
            push!(files_created, path)
        elseif ndims(array) == 2
            # 2D array: write directly
            path = write_array_to_csv(array, filename, delimiter=delimiter)
            push!(files_created, path)
        elseif ndims(array) == 3
            # 3D array: write slices to separate files
            paths = write_array_to_csv(array, filename, delimiter=delimiter)
            append!(files_created, paths)
        else
            error("Arrays with more than 3 dimensions are not supported")
        end
        
        println("Array '$key' written to CSV")
    end
    
    return files_created
end

"""
    save_capacity_results(V, bary, interface_length, cell_types, As, Ws, Bs, output_dir="./outputs")

Save capacity results to CSV files, handling 2D and 3D cases.

# Arguments
- `V`: Volume/Area array
- `bary`: Barycenter array
- `interface_length`: Interface length array
- `cell_types`: Cell types array
- `As`: Face apertures
- `Ws`: Staggered volumes
- `Bs`: Boundary fractions
- `output_dir`: Directory to save files

# Returns
- List of absolute paths of the created files.
"""
function save_capacity_results(V, bary, interface_length, cell_types, As, Ws, Bs, output_dir="./outputs")
    # Create dictionaries for the different types of data
    volumes = Dict("V" => V, "cell_types" => cell_types, "interface_length" => interface_length)
    
    # For bary (contains vectors), extract components
    if eltype(bary) <: StaticArrays.SVector
        # Determine dimensionality from bary
        dims = length(first(bary))
        
        # Extract each component
        for d in 1:dims
            volumes["bary_$(d)"] = getindex.(bary, d)
        end
    else
        volumes["bary"] = bary
    end
    
    # Add As, Ws, Bs components
    if As isa Tuple
        for (i, a) in enumerate(As)
            volumes["As_$(i)"] = a
        end
    else
        volumes["As"] = As
    end
    
    if Ws isa Tuple
        for (i, w) in enumerate(Ws)
            volumes["Ws_$(i)"] = w
        end
    else
        volumes["Ws"] = Ws
    end
    
    if Bs isa Tuple
        for (i, b) in enumerate(Bs)
            volumes["Bs_$(i)"] = b
        end
    else
        volumes["Bs"] = Bs
    end
    
    # Create base filename
    base_filename = joinpath(output_dir, "capacity_results")
    
    # Save all results to CSV files
    files = write_dict_of_arrays_to_csv(volumes, base_filename)
    println("Saved $(length(files)) capacity result files to $output_dir")
    
    return files
end

# Optional: Add support for visualizing 3D data
"""
    visualize_capacities_3d(capacities, time_index=1; title="3D Mesh Capacities", 
                            colormap=:viridis, figsize=(800, 600))

Visualize a time slice of 3D capacity data.

# Arguments
- `capacities`: Dictionary of 3D capacities with keys like "V", "Ax", etc.
- `time_index`: The time slice to visualize (default: 1)
- `title`: Main title of the figure
- `colormap`: Colormap to use for visualizations
- `figsize`: Size of the figure (width, height) in pixels
"""
function visualize_capacities_3d(capacities, time_index=1; 
                                title="3D Mesh Capacities", 
                                colormap=:viridis, 
                                figsize=(800, 600))
    n_capacities = length(capacities)
    if n_capacities == 0
        println("No capacities to visualize")
        return
    end
    
    # Determine subplot layout
    n_cols = min(3, n_capacities)
    n_rows = ceil(Int, n_capacities / n_cols)
    
    # Create the plot
    p = plot(layout=(n_rows, n_cols), size=figsize, 
             title="$title (Time slice: $time_index)", 
             titleloc=:center, titlefont=font(14))
    
    # Find max value for normalization
    vmax = maximum(maximum(arr[:,:,time_index]) for (_, arr) in capacities 
                   if ndims(arr) == 3 && time_index <= size(arr, 3) && !isempty(arr))
    
    i = 1
    for (name, array) in capacities
        if i <= n_rows * n_cols
            if ndims(array) == 3 && time_index <= size(array, 3)
                heatmap!(p, array[:,:,time_index], subplot=i, title=name, colorbar=true, 
                         color=colormap, clim=(0, vmax), xlabel="i", ylabel="j")
            elseif ndims(array) == 2
                heatmap!(p, array, subplot=i, title=name, colorbar=true, 
                         color=colormap, xlabel="i", ylabel="j")
            elseif ndims(array) == 1
                plot!(p, array, subplot=i, title=name, 
                      xlabel="Index", ylabel="Value", legend=false)
            end
            i += 1
        end
    end
    
    return p
end

# For your 3D case (2D + time)
files = save_capacity_results(V, bary, interface_length, cell_types, As, Ws, Bs, "./outputs_3d")

