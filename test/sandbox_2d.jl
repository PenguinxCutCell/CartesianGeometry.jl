using CartesianGeometry
using Test
const T=Float64
using StaticArrays
x = [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
     0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
y = [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
        0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
xyz = (x, y)

# define level set
const R = 0.4
const a, b = 0.5, 0.5

levelset = HyperSphere(R, (a, b))

V, bary, interface_length, cell_types = integrate(Tuple{0}, levelset, xyz, T, nan)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)

println(Bs)
readline()

using CSV
using DataFrames

"""
    write_array_to_csv(array, filename; delimiter=',', header=nothing, index=nothing)

Write a Julia array to a CSV file.

# Arguments
- `array`: Array to write to the file
- `filename`: Name of the output file. If the path doesn't exist, it will be created.
- `delimiter`: Delimiter to use in the CSV file. Default: ','
- `header`: List of column names. If nothing, no header will be written.
- `index`: List of row names. If nothing, no indices will be written.

# Returns
- The absolute path of the created file.
"""
function write_array_to_csv(array, filename; delimiter=',', header=nothing, index=nothing)
    # Ensure directory exists
    mkpath(dirname(abspath(filename)))
    
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
        
        # For 1D arrays, reshape to column vector for CSV writing
        if ndims(array) == 1
            reshaped_array = reshape(array, :, 1)
            write_array_to_csv(reshaped_array, filename, delimiter=delimiter)
        else
            write_array_to_csv(array, filename, delimiter=delimiter)
        end
        
        push!(files_created, abspath(filename))
        println("Array '$key' written to file: $filename")
    end
    
    return files_created
end

"""
    visualize_capacities(capacities; title="Mesh Capacities", colormap=:viridis, figsize=(800, 600))

Visualize different mesh capacities using subplots.

# Arguments
- `capacities`: Dictionary of capacities with keys like "V", "Ax", etc.
- `title`: Main title of the figure
- `colormap`: Colormap to use for visualizations
- `figsize`: Size of the figure (width, height) in pixels
"""
function visualize_capacities(capacities; title="Mesh Capacities", colormap=:viridis, figsize=(800, 600))
    n_capacities = length(capacities)
    if n_capacities == 0
        println("No capacities to visualize")
        return
    end
    
    # Determine subplot layout
    n_cols = min(3, n_capacities)
    n_rows = ceil(Int, n_capacities / n_cols)
    
    # Create the plot
    p = plot(layout=(n_rows, n_cols), size=figsize, title=title, titleloc=:center, titlefont=font(14))
    
    # Find max value for normalization
    vmax = maximum(maximum(arr) for (_, arr) in capacities if !isempty(arr))
    
    i = 1
    for (name, array) in capacities
        if i <= n_rows * n_cols
            if ndims(array) == 2
                heatmap!(p, array, subplot=i, title=name, colorbar=true, 
                         color=colormap, clim=(0, vmax), xlabel="i", ylabel="j")
            else
                # For 1D arrays, plot them as line plots
                plot!(p, array, subplot=i, title=name, 
                      xlabel="Index", ylabel="Value", legend=false)
            end
            i += 1
        end
    end
    
    return p
end

# Example usage with your CartesianGeometry results
function save_capacity_results(V, bary, interface_length, cell_types, As, Ws, Bs, output_dir="./outputs")
    # Create dictionaries for the different types of data
    volumes = Dict("V" => V, "cell_types" => cell_types, "interface_length" => interface_length)
    
    # For bary (contains vectors), extract components
    if eltype(bary) <: StaticArrays.SVector
        bary_x = getindex.(bary, 1)
        bary_y = getindex.(bary, 2)
        volumes["bary_x"] = bary_x
        volumes["bary_y"] = bary_y
    else
        volumes["bary"] = bary
    end
    
    # Add As, Ws, Bs components
    if As isa Tuple
        for (i, a) in enumerate(As)
            volumes["As_$(i)"] = a
        end
    end
    
    if Ws isa Tuple
        for (i, w) in enumerate(Ws)
            volumes["Ws_$(i)"] = w
        end
    end
    
    if Bs isa Tuple
        for (i, b) in enumerate(Bs)
            volumes["Bs_$(i)"] = b
        end
    end
    
    # Create base filename
    base_filename = joinpath(output_dir, "capacity_results")
    
    # Save all results to CSV files
    files = write_dict_of_arrays_to_csv(volumes, base_filename)
    println("Saved $(length(files)) capacity result files to $output_dir")
    
    return files
end

# After running your simulation
files = save_capacity_results(V, bary, interface_length, cell_types, As, Ws, Bs)