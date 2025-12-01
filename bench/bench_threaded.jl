"""
Benchmark script to compare serial vs threaded integration performance.

Run with multiple threads to see performance improvement:
    julia --threads=4 --project=. bench/bench_threaded.jl

Or with all available cores:
    julia --threads=auto --project=. bench/bench_threaded.jl
"""

using CartesianGeometry
using Printf

const T = Float64

# Helper to create 2D mesh
function create_mesh_2d(n)
    x = collect(range(0.0, 1.0, length=n+1))
    y = collect(range(0.0, 1.0, length=n+1))
    return (x, y)
end

# Helper to create 3D mesh  
function create_mesh_3d(n)
    x = collect(range(0.0, 1.0, length=n+1))
    y = collect(range(0.0, 1.0, length=n+1))
    z = collect(range(0.0, 1.0, length=n+1))
    return (x, y, z)
end

# Level set functions
const R = 0.4
const center = (0.5, 0.5, 0.5)
levelset_2d = HyperSphere(R, (0.5, 0.5))
levelset_3d = HyperSphere(R, (0.5, 0.5, 0.5))

function run_benchmark()
    println("="^70)
    println("CartesianGeometry.jl - Serial vs Threaded Performance Comparison")
    println("="^70)
    println()
    println("Number of threads: $(Threads.nthreads())")
    println()
    
    # Warm-up runs
    println("Warming up...")
    mesh_2d = create_mesh_2d(10)
    mesh_3d = create_mesh_3d(10)
    integrate(Tuple{0}, levelset_2d, mesh_2d, T, zero)
    integrate_threaded(Tuple{0}, levelset_2d, mesh_2d, T, zero)
    integrate(Tuple{0}, levelset_3d, mesh_3d, T, zero)
    integrate_threaded(Tuple{0}, levelset_3d, mesh_3d, T, zero)
    integrate(Tuple{1}, levelset_2d, mesh_2d, T, zero)
    integrate_threaded(Tuple{1}, levelset_2d, mesh_2d, T, zero)
    integrate(Tuple{1}, levelset_3d, mesh_3d, T, zero)
    integrate_threaded(Tuple{1}, levelset_3d, mesh_3d, T, zero)
    println("Warm-up complete.\n")
    
    println("-"^70)
    println("2D Volume Integration (Tuple{0})")
    println("-"^70)
    println(@sprintf("%-12s %12s %12s %12s %12s", "Resolution", "Serial (s)", "Threaded (s)", "Speedup", "Cells"))
    
    resolutions_2d = [20, 40, 60, 80, 100, 150, 200]
    
    for n in resolutions_2d
        mesh = create_mesh_2d(n)
        num_cells = n * n
        
        # Serial timing (multiple runs)
        GC.gc()
        t_serial = @elapsed for _ in 1:3
            integrate(Tuple{0}, levelset_2d, mesh, T, zero)
        end
        t_serial /= 3
        
        # Threaded timing (multiple runs)
        GC.gc()
        t_threaded = @elapsed for _ in 1:3
            integrate_threaded(Tuple{0}, levelset_2d, mesh, T, zero)
        end
        t_threaded /= 3
        
        speedup = t_serial / t_threaded
        println(@sprintf("%-12s %12.4f %12.4f %12.2fx %12d", "$(n)×$(n)", t_serial, t_threaded, speedup, num_cells))
    end
    
    println()
    println("-"^70)
    println("3D Volume Integration (Tuple{0})")
    println("-"^70)
    println(@sprintf("%-12s %12s %12s %12s %12s", "Resolution", "Serial (s)", "Threaded (s)", "Speedup", "Cells"))
    
    resolutions_3d = [10, 15, 20, 25, 30, 40, 50]
    
    for n in resolutions_3d
        mesh = create_mesh_3d(n)
        num_cells = n * n * n
        
        # Serial timing
        GC.gc()
        t_serial = @elapsed integrate(Tuple{0}, levelset_3d, mesh, T, zero)
        
        # Threaded timing
        GC.gc()
        t_threaded = @elapsed integrate_threaded(Tuple{0}, levelset_3d, mesh, T, zero)
        
        speedup = t_serial / t_threaded
        println(@sprintf("%-12s %12.4f %12.4f %12.2fx %12d", "$(n)×$(n)×$(n)", t_serial, t_threaded, speedup, num_cells))
    end
    
    println()
    println("-"^70)
    println("2D Surface Integration (Tuple{1})")
    println("-"^70)
    println(@sprintf("%-12s %12s %12s %12s %12s", "Resolution", "Serial (s)", "Threaded (s)", "Speedup", "Cells"))
    
    for n in resolutions_2d
        mesh = create_mesh_2d(n)
        num_cells = n * n
        
        # Serial timing
        GC.gc()
        t_serial = @elapsed for _ in 1:3
            integrate(Tuple{1}, levelset_2d, mesh, T, zero)
        end
        t_serial /= 3
        
        # Threaded timing
        GC.gc()
        t_threaded = @elapsed for _ in 1:3
            integrate_threaded(Tuple{1}, levelset_2d, mesh, T, zero)
        end
        t_threaded /= 3
        
        speedup = t_serial / t_threaded
        println(@sprintf("%-12s %12.4f %12.4f %12.2fx %12d", "$(n)×$(n)", t_serial, t_threaded, speedup, num_cells))
    end
    
    println()
    println("-"^70)
    println("3D Surface Integration (Tuple{1})")
    println("-"^70)
    println(@sprintf("%-12s %12s %12s %12s %12s", "Resolution", "Serial (s)", "Threaded (s)", "Speedup", "Cells"))
    
    for n in resolutions_3d
        mesh = create_mesh_3d(n)
        num_cells = n * n * n
        
        # Serial timing
        GC.gc()
        t_serial = @elapsed integrate(Tuple{1}, levelset_3d, mesh, T, zero)
        
        # Threaded timing
        GC.gc()
        t_threaded = @elapsed integrate_threaded(Tuple{1}, levelset_3d, mesh, T, zero)
        
        speedup = t_serial / t_threaded
        println(@sprintf("%-12s %12.4f %12.4f %12.2fx %12d", "$(n)×$(n)×$(n)", t_serial, t_threaded, speedup, num_cells))
    end
    
    println()
    println("-"^70)
    println("Second Kind Integration (3D)")
    println("-"^70)
    println(@sprintf("%-12s %12s %12s %12s %12s", "Resolution", "Serial (s)", "Threaded (s)", "Speedup", "Cells"))
    
    for n in resolutions_3d[1:5]  # Use smaller resolutions for second kind
        mesh = create_mesh_3d(n)
        num_cells = n * n * n
        
        # First compute barycenters
        _, bary, _, _, _ = integrate(Tuple{0}, levelset_3d, mesh, T, zero)
        
        # Serial timing for second kind
        GC.gc()
        t_serial = @elapsed integrate(Tuple{0}, levelset_3d, mesh, T, zero, bary)
        
        # Threaded timing for second kind
        GC.gc()
        t_threaded = @elapsed integrate_threaded(Tuple{0}, levelset_3d, mesh, T, zero, bary)
        
        speedup = t_serial / t_threaded
        println(@sprintf("%-12s %12.4f %12.4f %12.2fx %12d", "$(n)×$(n)×$(n)", t_serial, t_threaded, speedup, num_cells))
    end
    
    println()
    println("="^70)
    println("Benchmark complete!")
    println("="^70)
end

# Run the benchmark
run_benchmark()
