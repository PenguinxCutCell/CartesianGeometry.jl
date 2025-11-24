# Utilities
isfull(val, full_val) = isapprox(val, full_val; atol=1e-8)
isempty(val) = isapprox(val, 0.0; atol=1e-8)

########################
# 1D Implementation
########################

function implicit_integration(mesh::Tuple{Vector}, Φ; tol=1e-6)
    x_coords = mesh[1]
    nx = length(x_coords)-1
    dx = x_coords[2] - x_coords[1]

    # Compute volumes (1D length segments)
    V = zeros(nx+1)
    for i in 1:nx
        a = (x_coords[i],)
        b = (x_coords[i+1],)
        V[i] = ImplicitIntegration.integrate((x)->1, Φ, a, b;tol=tol).val
    end

    # Classify cells
    cell_types = zeros(nx+1)
    for i in 1:nx
        if isempty(V[i])
            cell_types[i] = 0
        elseif isfull(V[i], dx)
            cell_types[i] = 1
        else
            cell_types[i] = -1
        end
    end

    # Compute cell centroids
    C_ω = zeros(nx+1) # cell centroids
    for i in 1:nx
        a = (x_coords[i],)
        b = (x_coords[i+1],)
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b;tol=tol).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b;tol=tol).val / area
            C_ω[i] = isnan(x_c) ? 0.5*(x_coords[i]+x_coords[i+1]) : x_c
        else
            C_ω[i] = 0.5*(x_coords[i]+x_coords[i+1])
        end
    end

    # Compute interface centroids
    C_γ = zeros(nx+1)
    Γ = zeros(nx+1)
    for i in 1:nx
        a = (x_coords[i],)
        b = (x_coords[i+1],)
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b, surface=true;tol=tol).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b, surface=true;tol=tol).val / area
            C_γ[i] = x_c
            Γ[i] = area
        else
            C_γ[i] = NaN
            Γ[i] = 0.0
        end
    end
    # Compute W (staggered volumes) based on cell centroids
    Wx = zeros(nx+1)
    for i in 1:nx+1
        xi = C_ω[max(i-1,1)]
        xip1 = C_ω[min(i,nx)]
        Wx[i] = ImplicitIntegration.integrate((x)->1, Φ, (xi,), (xip1,);tol=tol).val
    end
    W = (Wx,)

    # Compute A (faces capacity)
    Ax = zeros(nx+1)
    for i in 1:nx+1
        Ax[i] = Φ((x_coords[i],)) ≤ 0.0 ? 1.0 : 0.0
    end
    A = (Ax,)

    # Compute B (values at cell centroids):
    Bx = zeros(nx+1)
    for i in 1:nx
        xi = C_ω[i]
        Bx[i] = Φ((xi,)) ≤ 0.0 ? 1.0 : 0.0
    end
    B = (Bx,)

    # Convert C_ω to Vector of SVector
    C_ω = [SVector{1,Float64}([x]) for x in C_ω]
    C_γ = [SVector{1,Float64}([x]) for x in C_γ]

    return V, cell_types, C_ω, C_γ, Γ, W, A, B
end

########################
# 2D Implementation
########################

function implicit_integration(mesh::Tuple{Vector,Vector}, Φ; tol=1e-6)
    x_coords, y_coords = mesh
    nx = length(x_coords)-1
    ny = length(y_coords)-1
    dx = x_coords[2]-x_coords[1]
    dy = y_coords[2]-y_coords[1]
    dim = 2

    # Total number of cells in flattened array
    total_cells = (nx+1)*(ny+1)
    
    # Linear indexing function for consistent use throughout
    linear_idx(i,j) = i + (j-1)*(nx+1)

    # Volumes (areas)
    V = zeros(total_cells)
    for i in 1:nx
        for j in 1:ny
            idx = linear_idx(i, j)
            a = (x_coords[i], y_coords[j])
            b = (x_coords[i+1], y_coords[j+1])
            V[idx] = ImplicitIntegration.integrate((x)->1, Φ, a, b;tol=tol).val
        end
    end

    # Cell types
    cell_types = zeros(Int, total_cells)
    for i in 1:nx
        for j in 1:ny
            idx = linear_idx(i, j)
            vol = V[idx]
            if isempty(vol)
                cell_types[idx] = 0
            elseif isfull(vol, dx*dy)
                cell_types[idx] = 1
            else
                cell_types[idx] = -1
            end
        end
    end

    # Cell centroids
    C_ω = [SVector{dim,Float64}(zeros(dim)) for _ in 1:total_cells]
    for i in 1:nx
        for j in 1:ny
            idx = linear_idx(i, j)
            a = (x_coords[i], y_coords[j])
            b = (x_coords[i+1], y_coords[j+1])
            area = V[idx]  # Reuse the already computed volume
            if area > 0
                x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b;tol=tol).val / area
                y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b;tol=tol).val / area
                x_c = isnan(x_c) ? 0.5*(x_coords[i]+x_coords[i+1]) : x_c
                y_c = isnan(y_c) ? 0.5*(y_coords[j]+y_coords[j+1]) : y_c
                C_ω[idx] = SVector{dim,Float64}(x_c, y_c)
            else
                C_ω[idx] = SVector{dim,Float64}(0.5*(x_coords[i]+x_coords[i+1]), 0.5*(y_coords[j]+y_coords[j+1]))
            end
        end
    end

    # Interface centroids and lengths
    C_γ = [SVector{dim,Float64}(zeros(dim)) for _ in 1:total_cells]
    Γ = zeros(total_cells)
    for i in 1:nx
        for j in 1:ny
            idx = linear_idx(i, j)
            a = (x_coords[i], y_coords[j])
            b = (x_coords[i+1], y_coords[j+1])
            area = ImplicitIntegration.integrate((x)->1, Φ, a, b, surface=true;tol=tol).val
            if area > 0
                x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b, surface=true;tol=tol).val / area
                y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b, surface=true;tol=tol).val / area
                C_γ[idx] = SVector{dim,Float64}(x_c, y_c)
                Γ[idx] = area
            else
                C_γ[idx] = SVector{dim,Float64}(NaN, NaN)
                Γ[idx] = 0.0
            end
        end
    end

    # Face capacities A = (Ax, Ay, Az)
    Ax = zeros(total_cells)
    Ay = zeros(total_cells)
    
    for i in 1:nx+1, j in 1:ny
        idx = linear_idx(i, j)
        xi = x_coords[i]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        
        # Creating a 1D level set function where x is fixed
        ϕ_2d = (y) -> Φ((xi, y[1]))
        Ax[idx] = ImplicitIntegration.integrate((y)->1, ϕ_2d, (yj,), (yjp1,);tol=tol).val
    end
    
    for i in 1:nx, j in 1:ny+1
        idx = linear_idx(i, j)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        
        # Creating a 2D level set function where y is fixed
        ϕ_2d = (x) -> Φ((x[1], yj))
        Ay[idx] = ImplicitIntegration.integrate((x)->1, ϕ_2d, (xi,), (xip1,);tol=tol).val
    end
    
    A = (Ax, Ay)

    # Cell boundary fractions B = (Bx, By, Bz)
    Bx = zeros(total_cells)
    By = zeros(total_cells)
    
    for i in 1:nx, j in 1:ny
        idx = linear_idx(i, j)
        c_x = C_ω[idx][1]  # x coordinate of cell centroid
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        
        # Fixed x at centroid
        ϕ_2d = (y) -> Φ((c_x, y[1]))
        Bx[idx] = ImplicitIntegration.integrate((y)->1, ϕ_2d, (yj,), (yjp1,);tol=tol).val
    end
    
    for i in 1:nx, j in 1:ny
        idx = linear_idx(i, j)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        c_y = C_ω[idx][2]  # y coordinate of cell centroid
        
        # Fixed y at centroid
        ϕ_2d = (x) -> Φ((x[1], c_y))
        By[idx] = ImplicitIntegration.integrate((x)->1, ϕ_2d, (xi,), (xip1,);tol=tol).val
    end

    B = (Bx, By)

    # Staggered volumes W = (Wx, Wy)
    Wx = zeros((nx+1)*(ny+1))
    Wy = zeros((nx+1)*(ny+1))
    
    # For Wx (staggered volume in x-direction)
    for i in 1:nx+1  # Include ghost cell indices for boundaries
        for j in 1:ny
            idx = linear_idx(i, j)
            
            # Get neighboring cell centroids
            left_i = max(1, i-1)
            right_i = min(i, nx)
            
            left_idx = linear_idx(left_i, j)
            right_idx = linear_idx(right_i, j)
            
            # Use centroids for x-coordinates, grid points for y-coordinates
            a = (C_ω[left_idx][1], y_coords[j])
            b = (C_ω[right_idx][1], y_coords[j+1])
            
            # Integrate volume between centroids
            Wx[idx] = ImplicitIntegration.integrate((x)->1.0, Φ, a, b;tol=tol).val
        end
    end
    
    # For Wy (staggered volume in y-direction)
    for i in 1:nx
        for j in 1:ny+1  # Include ghost cell indices for boundaries
            idx = linear_idx(i, j)
            
            # Get neighboring cell centroids
            bottom_j = max(1, j-1)
            top_j = min(j, ny)
            
            bottom_idx = linear_idx(i, bottom_j)
            top_idx = linear_idx(i, top_j)
            
            # Use centroids for y-coordinates, grid points for x-coordinates
            a = (x_coords[i], C_ω[bottom_idx][2])
            b = (x_coords[i+1], C_ω[top_idx][2])
            
            # Integrate volume between centroids
            Wy[idx] = ImplicitIntegration.integrate((x)->1.0, Φ, a, b;tol=tol).val
        end
    end
    
    W = (Wx, Wy)

    return V, cell_types, C_ω, C_γ, Γ, W, A, B
end

########################
# 3D Implementation
########################

function implicit_integration(mesh::Tuple{Vector,Vector,Vector}, Φ;tol=1e-6)
    x_coords, y_coords, z_coords = mesh
    nx = length(x_coords)-1
    ny = length(y_coords)-1
    nz = length(z_coords)-1
    dx = x_coords[2]-x_coords[1]
    dy = y_coords[2]-y_coords[1]
    dz = z_coords[2]-z_coords[1]
    dim = 3
    
    # Total number of cells in flattened array
    total_cells = (nx+1)*(ny+1)*(nz+1)
    
    # Linear indexing function
    function linear_idx(i, j, k)
        i + (j-1)*(nx+1) + (k-1)*(nx+1)*(ny+1)
    end

    # Volume
    V = zeros(total_cells)
    for i in 1:nx, j in 1:ny, k in 1:nz
        idx = linear_idx(i, j, k)
        a = (x_coords[i], y_coords[j], z_coords[k])
        b = (x_coords[i+1], y_coords[j+1], z_coords[k+1])
        V[idx] = ImplicitIntegration.integrate((x)->1, Φ, a, b;tol=tol).val
    end

    # Cell types
    cell_types = zeros(total_cells)
    for i in 1:nx, j in 1:ny, k in 1:nz
        idx = linear_idx(i, j, k)
        vol = V[idx]
        if isempty(vol)
            cell_types[idx] = 0
        elseif isfull(vol, dx*dy*dz)
            cell_types[idx] = 1
        else
            cell_types[idx] = -1
        end
    end

    # Cell centroids
    C_ω = [SVector{dim,Float64}(zeros(dim)) for _ in 1:total_cells]
    for i in 1:nx, j in 1:ny, k in 1:nz
        idx = linear_idx(i, j, k)
        a = (x_coords[i], y_coords[j], z_coords[k])
        b = (x_coords[i+1], y_coords[j+1], z_coords[k+1])
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b;tol=tol).val / area
            y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b;tol=tol).val / area
            z_c = ImplicitIntegration.integrate((x)->x[3], Φ, a, b;tol=tol).val / area
            x_c = isnan(x_c) ? 0.5*(x_coords[i]+x_coords[i+1]) : x_c
            y_c = isnan(y_c) ? 0.5*(y_coords[j]+y_coords[j+1]) : y_c
            z_c = isnan(z_c) ? 0.5*(z_coords[k]+z_coords[k+1]) : z_c
            C_ω[idx] = SVector{dim,Float64}(x_c, y_c, z_c)
        else
            C_ω[idx] = SVector{dim,Float64}(
                0.5*(x_coords[i]+x_coords[i+1]),
                0.5*(y_coords[j]+y_coords[j+1]),
                0.5*(z_coords[k]+z_coords[k+1])
            )
        end
    end

    # Interface centroids
    C_γ = [SVector{dim,Float64}(zeros(dim)) for _ in 1:total_cells]
    Γ = zeros(total_cells)
    for i in 1:nx, j in 1:ny, k in 1:nz
        idx = linear_idx(i, j, k)
        a = (x_coords[i], y_coords[j], z_coords[k])
        b = (x_coords[i+1], y_coords[j+1], z_coords[k+1])
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b, surface=true;tol=tol).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b, surface=true;tol=tol).val / area
            y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b, surface=true;tol=tol).val / area
            z_c = ImplicitIntegration.integrate((x)->x[3], Φ, a, b, surface=true;tol=tol).val / area
            C_γ[idx] = SVector{dim,Float64}(x_c, y_c, z_c)
            Γ[idx] = area
        else
            C_γ[idx] = SVector{dim,Float64}(NaN, NaN, NaN)
            Γ[idx] = 0.0
        end
    end

    # Staggered volumes W = (Wx, Wy, Wz)
    Wx = zeros(total_cells)
    Wy = zeros(total_cells)
    Wz = zeros(total_cells)
    
    for i in 1:nx+1, j in 1:ny, k in 1:nz
        idx = linear_idx(i, j, k)
        prev_i = max(i-1, 1)
        next_i = min(i, nx)
        
        prev_idx = linear_idx(prev_i, j, k)
        next_idx = linear_idx(next_i, j, k)
        
        xi = C_ω[prev_idx][1]
        xip1 = C_ω[next_idx][1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        
        Wx[idx] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zk), (xip1,yjp1,zkp1);tol=tol).val
    end
    
    for i in 1:nx, j in 1:ny+1, k in 1:nz
        idx = linear_idx(i, j, k)
        prev_j = max(j-1, 1)
        next_j = min(j, ny)
        
        prev_idx = linear_idx(i, prev_j, k)
        next_idx = linear_idx(i, next_j, k)
        
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = C_ω[prev_idx][2]
        yjp1 = C_ω[next_idx][2]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        
        Wy[idx] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zk), (xip1,yjp1,zkp1);tol=tol).val
    end
    
    for i in 1:nx, j in 1:ny, k in 1:nz+1
        idx = linear_idx(i, j, k)
        prev_k = max(k-1, 1)
        next_k = min(k, nz)
        
        prev_idx = linear_idx(i, j, prev_k)
        next_idx = linear_idx(i, j, next_k)
        
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = C_ω[prev_idx][3]
        zkp1 = C_ω[next_idx][3]
        
        Wz[idx] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zk), (xip1,yjp1,zkp1);tol=tol).val
    end
    
    W = (Wx, Wy, Wz)

    # Face capacities A = (Ax, Ay, Az)
    Ax = zeros(total_cells)
    Ay = zeros(total_cells)
    Az = zeros(total_cells)
    
    for i in 1:nx+1, j in 1:ny, k in 1:nz
        idx = linear_idx(i, j, k)
        xi = x_coords[i]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        
        # Creating a 2D level set function where x is fixed
        ϕ_2d = (y) -> Φ((xi, y[1], y[2]))
        Ax[idx] = ImplicitIntegration.integrate((y)->1, ϕ_2d, (yj,zk), (yjp1,zkp1);tol=tol).val
    end
    
    for i in 1:nx, j in 1:ny+1, k in 1:nz
        idx = linear_idx(i, j, k)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        
        # Creating a 2D level set function where y is fixed
        ϕ_2d = (x) -> Φ((x[1], yj, x[2]))
        Ay[idx] = ImplicitIntegration.integrate((x)->1, ϕ_2d, (xi,zk), (xip1,zkp1);tol=tol).val
    end
    
    for i in 1:nx, j in 1:ny, k in 1:nz+1
        idx = linear_idx(i, j, k)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        
        # Creating a 2D level set function where z is fixed
        ϕ_2d = (x) -> Φ((x[1], x[2], zk))
        Az[idx] = ImplicitIntegration.integrate((x)->1, ϕ_2d, (xi,yj), (xip1,yjp1);tol=tol).val
    end
    
    A = (Ax, Ay, Az)

    # Cell boundary fractions B = (Bx, By, Bz)
    Bx = zeros(total_cells)
    By = zeros(total_cells)
    Bz = zeros(total_cells)
    
    for i in 1:nx, j in 1:ny, k in 1:nz
        idx = linear_idx(i, j, k)
        c_x = C_ω[idx][1]  # x coordinate of cell centroid
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        
        # Fixed x at centroid
        ϕ_2d = (y) -> Φ((c_x, y[1], y[2]))
        Bx[idx] = ImplicitIntegration.integrate((y)->1, ϕ_2d, (yj,zk), (yjp1,zkp1);tol=tol).val
    end
    
    for i in 1:nx, j in 1:ny, k in 1:nz
        idx = linear_idx(i, j, k)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        c_y = C_ω[idx][2]  # y coordinate of cell centroid
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        
        # Fixed y at centroid
        ϕ_2d = (x) -> Φ((x[1], c_y, x[2]))
        By[idx] = ImplicitIntegration.integrate((x)->1, ϕ_2d, (xi,zk), (xip1,zkp1);tol=tol).val
    end
    
    for i in 1:nx, j in 1:ny, k in 1:nz
        idx = linear_idx(i, j, k)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        c_z = C_ω[idx][3]  # z coordinate of cell centroid
        
        # Fixed z at centroid
        ϕ_2d = (x) -> Φ((x[1], x[2], c_z))
        Bz[idx] = ImplicitIntegration.integrate((x)->1, ϕ_2d, (xi,yj), (xip1,yjp1);tol=tol).val
    end
    
    B = (Bx, By, Bz)

    return V, cell_types, C_ω, C_γ, Γ, W, A, B
end

########################
# 4D Implementation
########################

function implicit_integration(mesh::Tuple{Vector,Vector,Vector,Vector}, Φ;tol=1e-6)
    x_coords, y_coords, z_coords, t_coords = mesh
    nx = length(x_coords)-1
    ny = length(y_coords)-1
    nz = length(z_coords)-1
    nt = length(t_coords)-1
    dx = x_coords[2]-x_coords[1]
    dy = y_coords[2]-y_coords[1]
    dz = z_coords[2]-z_coords[1]
    dt = t_coords[2]-t_coords[1]
    dim = 4
    
    # Total number of cells in flattened array
    total_cells = (nx+1)*(ny+1)*(nz+1)*(nt+1)
    
    # Linear indexing function
    function linear_idx(i, j, k, l)
        i + (j-1)*(nx+1) + (k-1)*(nx+1)*(ny+1) + (l-1)*(nx+1)*(ny+1)*(nz+1)
    end

    # Volume (4D hypervolume)
    V = zeros(total_cells)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        a = (x_coords[i], y_coords[j], z_coords[k], t_coords[l])
        b = (x_coords[i+1], y_coords[j+1], z_coords[k+1], t_coords[l+1])
        V[idx] = ImplicitIntegration.integrate((x)->1, Φ, a, b;tol=tol).val
    end

    # Cell types
    cell_types = zeros(Int, total_cells)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        vol = V[idx]
        if isempty(vol)
            cell_types[idx] = 0
        elseif isfull(vol, dx*dy*dz*dt)
            cell_types[idx] = 1
        else
            cell_types[idx] = -1
        end
    end

    # Cell centroids
    C_ω = [SVector{dim,Float64}(zeros(dim)) for _ in 1:total_cells]
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        a = (x_coords[i], y_coords[j], z_coords[k], t_coords[l])
        b = (x_coords[i+1], y_coords[j+1], z_coords[k+1], t_coords[l+1])
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b;tol=tol).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b;tol=tol).val / area
            y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b;tol=tol).val / area
            z_c = ImplicitIntegration.integrate((x)->x[3], Φ, a, b;tol=tol).val / area
            t_c = ImplicitIntegration.integrate((x)->x[4], Φ, a, b;tol=tol).val / area
            x_c = isnan(x_c) ? 0.5*(x_coords[i]+x_coords[i+1]) : x_c
            y_c = isnan(y_c) ? 0.5*(y_coords[j]+y_coords[j+1]) : y_c
            z_c = isnan(z_c) ? 0.5*(z_coords[k]+z_coords[k+1]) : z_c
            t_c = isnan(t_c) ? 0.5*(t_coords[l]+t_coords[l+1]) : t_c
            C_ω[idx] = SVector{dim,Float64}(x_c, y_c, z_c, t_c)
        else
            C_ω[idx] = SVector{dim,Float64}(
                0.5*(x_coords[i]+x_coords[i+1]),
                0.5*(y_coords[j]+y_coords[j+1]),
                0.5*(z_coords[k]+z_coords[k+1]),
                0.5*(t_coords[l]+t_coords[l+1])
            )
        end
    end

    # Interface centroids (hypervolume boundary)
    C_γ = [SVector{dim,Float64}(zeros(dim)) for _ in 1:total_cells]
    Γ = zeros(total_cells)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        a = (x_coords[i], y_coords[j], z_coords[k], t_coords[l])
        b = (x_coords[i+1], y_coords[j+1], z_coords[k+1], t_coords[l+1])
        area = ImplicitIntegration.integrate((x)->1, Φ, a, b, surface=true;tol=tol).val
        if area > 0
            x_c = ImplicitIntegration.integrate((x)->x[1], Φ, a, b, surface=true;tol=tol).val / area
            y_c = ImplicitIntegration.integrate((x)->x[2], Φ, a, b, surface=true;tol=tol).val / area
            z_c = ImplicitIntegration.integrate((x)->x[3], Φ, a, b, surface=true;tol=tol).val / area
            t_c = ImplicitIntegration.integrate((x)->x[4], Φ, a, b, surface=true;tol=tol).val / area
            C_γ[idx] = SVector{dim,Float64}(x_c, y_c, z_c, t_c)
            Γ[idx] = area
        else
            C_γ[idx] = SVector{dim,Float64}(NaN, NaN, NaN, NaN)
            Γ[idx] = 0.0
        end
    end

    # Staggered volumes W = (Wx, Wy, Wz, Wt)
    Wx = zeros(total_cells)
    Wy = zeros(total_cells)
    Wz = zeros(total_cells)
    Wt = zeros(total_cells)
    
    # For Wx (staggered volume in x-direction)
    for i in 1:nx+1, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        prev_i = max(i-1, 1)
        next_i = min(i, nx)
        
        prev_idx = linear_idx(prev_i, j, k, l)
        next_idx = linear_idx(next_i, j, k, l)
        
        xi = C_ω[prev_idx][1]
        xip1 = C_ω[next_idx][1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        tl = t_coords[l]
        tlp1 = t_coords[l+1]
        
        Wx[idx] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zk,tl), (xip1,yjp1,zkp1,tlp1);tol=tol).val
    end
    
    # For Wy (staggered volume in y-direction)
    for i in 1:nx, j in 1:ny+1, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        prev_j = max(j-1, 1)
        next_j = min(j, ny)
        
        prev_idx = linear_idx(i, prev_j, k, l)
        next_idx = linear_idx(i, next_j, k, l)
        
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = C_ω[prev_idx][2]
        yjp1 = C_ω[next_idx][2]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        tl = t_coords[l]
        tlp1 = t_coords[l+1]
        
        Wy[idx] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zk,tl), (xip1,yjp1,zkp1,tlp1);tol=tol).val
    end
    
    # For Wz (staggered volume in z-direction)
    for i in 1:nx, j in 1:ny, k in 1:nz+1, l in 1:nt
        idx = linear_idx(i, j, k, l)
        prev_k = max(k-1, 1)
        next_k = min(k, nz)
        
        prev_idx = linear_idx(i, j, prev_k, l)
        next_idx = linear_idx(i, j, next_k, l)
        
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = C_ω[prev_idx][3]
        zkp1 = C_ω[next_idx][3]
        tl = t_coords[l]
        tlp1 = t_coords[l+1]
        
        Wz[idx] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zk,tl), (xip1,yjp1,zkp1,tlp1);tol=tol).val
    end
    
    # For Wt (staggered volume in t-direction)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt+1
        idx = linear_idx(i, j, k, l)
        prev_l = max(l-1, 1)
        next_l = min(l, nt)
        
        prev_idx = linear_idx(i, j, k, prev_l)
        next_idx = linear_idx(i, j, k, next_l)
        
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        tl = C_ω[prev_idx][4]
        tlp1 = C_ω[next_idx][4]
        
        Wt[idx] = ImplicitIntegration.integrate((x)->1, Φ, (xi,yj,zk,tl), (xip1,yjp1,zkp1,tlp1);tol=tol).val
    end
    
    W = (Wx, Wy, Wz, Wt)

    # Face capacities A = (Ax, Ay, Az, At)
    Ax = zeros(total_cells)
    Ay = zeros(total_cells)
    Az = zeros(total_cells)
    At = zeros(total_cells)
    
    # For Ax (x-face capacity - YZT hyperplane at constant x)
    for i in 1:nx+1, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        xi = x_coords[i]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        tl = t_coords[l]
        tlp1 = t_coords[l+1]
        
        # Creating a 3D level set function where x is fixed
        ϕ_3d = (y) -> Φ((xi, y[1], y[2], y[3]))
        Ax[idx] = ImplicitIntegration.integrate((y)->1, ϕ_3d, (yj,zk,tl), (yjp1,zkp1,tlp1);tol=tol).val
    end
    
    # For Ay (y-face capacity - XZT hyperplane at constant y)
    for i in 1:nx, j in 1:ny+1, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        tl = t_coords[l]
        tlp1 = t_coords[l+1]
        
        # Creating a 3D level set function where y is fixed
        ϕ_3d = (x) -> Φ((x[1], yj, x[2], x[3]))
        Ay[idx] = ImplicitIntegration.integrate((x)->1, ϕ_3d, (xi,zk,tl), (xip1,zkp1,tlp1);tol=tol).val
    end
    
    # For Az (z-face capacity - XYT hyperplane at constant z)
    for i in 1:nx, j in 1:ny, k in 1:nz+1, l in 1:nt
        idx = linear_idx(i, j, k, l)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        tl = t_coords[l]
        tlp1 = t_coords[l+1]
        
        # Creating a 3D level set function where z is fixed
        ϕ_3d = (x) -> Φ((x[1], x[2], zk, x[3]))
        Az[idx] = ImplicitIntegration.integrate((x)->1, ϕ_3d, (xi,yj,tl), (xip1,yjp1,tlp1);tol=tol).val
    end
    
    # For At (t-face capacity - XYZ hyperplane at constant t)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt+1
        idx = linear_idx(i, j, k, l)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        tl = t_coords[l]
        
        # Creating a 3D level set function where t is fixed
        ϕ_3d = (x) -> Φ((x[1], x[2], x[3], tl))
        At[idx] = ImplicitIntegration.integrate((x)->1, ϕ_3d, (xi,yj,zk), (xip1,yjp1,zkp1);tol=tol).val
    end
    
    A = (Ax, Ay, Az, At)

    # Cell boundary fractions B = (Bx, By, Bz, Bt)
    Bx = zeros(total_cells)
    By = zeros(total_cells)
    Bz = zeros(total_cells)
    Bt = zeros(total_cells)
    
    # For Bx (x-boundary fraction at centroid)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        c_x = C_ω[idx][1]  # x coordinate of cell centroid
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        tl = t_coords[l]
        tlp1 = t_coords[l+1]
        
        # Fixed x at centroid
        ϕ_3d = (y) -> Φ((c_x, y[1], y[2], y[3]))
        Bx[idx] = ImplicitIntegration.integrate((y)->1, ϕ_3d, (yj,zk,tl), (yjp1,zkp1,tlp1);tol=tol).val
    end
    
    # For By (y-boundary fraction at centroid)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        c_y = C_ω[idx][2]  # y coordinate of cell centroid
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        tl = t_coords[l]
        tlp1 = t_coords[l+1]
        
        # Fixed y at centroid
        ϕ_3d = (x) -> Φ((x[1], c_y, x[2], x[3]))
        By[idx] = ImplicitIntegration.integrate((x)->1, ϕ_3d, (xi,zk,tl), (xip1,zkp1,tlp1);tol=tol).val
    end
    
    # For Bz (z-boundary fraction at centroid)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        c_z = C_ω[idx][3]  # z coordinate of cell centroid
        tl = t_coords[l]
        tlp1 = t_coords[l+1]
        
        # Fixed z at centroid
        ϕ_3d = (x) -> Φ((x[1], x[2], c_z, x[3]))
        Bz[idx] = ImplicitIntegration.integrate((x)->1, ϕ_3d, (xi,yj,tl), (xip1,yjp1,tlp1);tol=tol).val
    end
    
    # For Bt (t-boundary fraction at centroid)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nt
        idx = linear_idx(i, j, k, l)
        xi = x_coords[i]
        xip1 = x_coords[i+1]
        yj = y_coords[j]
        yjp1 = y_coords[j+1]
        zk = z_coords[k]
        zkp1 = z_coords[k+1]
        c_t = C_ω[idx][4]  # t coordinate of cell centroid
        
        # Fixed t at centroid
        ϕ_3d = (x) -> Φ((x[1], x[2], x[3], c_t))
        Bt[idx] = ImplicitIntegration.integrate((x)->1, ϕ_3d, (xi,yj,zk), (xip1,yjp1,zkp1);tol=tol).val
    end
    
    B = (Bx, By, Bz, Bt)

    return V, cell_types, C_ω, C_γ, Γ, W, A, B
end