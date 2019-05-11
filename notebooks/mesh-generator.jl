using LinearAlgebra, Triangle

function delaunay(p)
    """
    Generates Delaunay triangulation of array of nodes
    
    Arguments:
        p - array of nodes
    """
    vertex_map = Array{Int64,1}(1:length(p))
    t = Triangle.basic_triangulation(convert(Array{Float64}, hcat(p...)'), vertex_map)
end

function inpolygon(p, pv)
    """
    Returns boolean to indicate whether or not a point is inside
    of a closed polygon
    
    Arguments:
        p  - point to be considered
        pv - list of points that defines a closed polygon
    """
    crossings = 0
    for i in 1:length(pv)-1
        if ((pv[i][2] <= p[2]) && (pv[i+1][2] > p[2])) || ((pv[i][2] > p[2]) && (pv[i+1][2] <= p[2]))            
            intersection = (p[2] - pv[i][2]) / (pv[i+1][2] - pv[i][2]) # intersection x coord
            if p[1] < pv[i][1] + intersection * (pv[i+1][1] - pv[i][1])
                crossings += 1
            end
        end
    end
    return crossings % 2 == 1
end

function tri_area(tri)
    """
    Calculates the area of a triangle
    
    Arguments:
        tri - list of points that defines a triangle, not closed
    """
    a, b, c = tri
    return abs((a[1] * (b[2] - c[2]) + b[1] * (c[2] - a[2]) + c[1] * (a[2] - b[2]))/2)
end

function tri_centroid(tri)
    """
    Finds the centroid of a triangle
    
    Arguments:
        tri - list of points that defines a triangle, not closed
    """
    a, b, c = tri
    return [1/3*(a[1] + b[1] + c[1]), 1/3*(a[2] + b[2] + c[2])]
end

function tri_circumcenter(tri)
    """
    Finds the circumcenter of a triangle
    
    Arguments:
        tri - list of points that defines a triangle, not closed
    """
    a, b, c = tri
    d = 2 * (a[1]*(b[2] - c[2]) + b[1]*(c[2] - a[2]) + c[1]*(a[2] - b[2]))
    ux = 1 / d * ((a[1]^2 + a[2]^2)*(b[2] - c[2]) + (b[1]^2 + b[2]^2)*(c[2] - a[2]) + (c[1]^2 + c[2]^2)*(a[2] - b[2]))
    uy = 1 / d * ((a[1]^2 + a[2]^2)*(c[1] - b[1]) + (b[1]^2 + b[2]^2)*(a[1] - c[1]) + (c[1]^2 + c[2]^2)*(b[1] - a[1]))
    return [ux, uy]
end

function dist(p1, p2)
    """
    Calculates the distance between two points
    
    Arguments:
        p1 - first point
        p2 - second point
    """
    return sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
end

function generate_nodes(p1, p2, hmax)
    """
    Generates nodes for Delaunay triangulation
    
    Arguments:
        p1   - first endpoint
        p2   - second endpoint
        hmax - parameter to determine density of nodes
    """
    distance = dist(p1, p2)
    
    # number of segments
    m = distance / hmax
    if m != m ÷ 1
        m += 1
    end
    n = m ÷ 1
    
    # distance between nodes
    d = distance / n
    
    # normed vector p1 → p2
    u = (p2 - p1) ./ distance
    
    # create nodes
    nodes = []
    for k = 1:n-1
        p3 = p1 + (k * d .* u)
        push!(nodes, p3)
    end
    
    return nodes
end

function remove_outside_triangles!(t, p, pv)
    """
    Removes triangles whose centroids are outside of the polygon
    
    Arguments:
        t  - triangulation of nodes `p`
        p  - nodes of polygon
        pv - closed polygon
    """
    idx_to_delete = []
    for tri in t
        centroid = tri_centroid(p[tri])
        if !(inpolygon(centroid, pv))
            idx = findfirst(t .== [tri])
            push!(idx_to_delete, idx)
        end
    end
    deleteat!(t, idx_to_delete)
end

function largest_area(t, p)
    """
    Determines the area and indices of triangle with largest area
    
    Arguments:
        t - triangulation of nodes `p`
        p - nodes of polygon
    """
    A, A_tri = 0, []
    for tri in t
        area = tri_area(p[tri])
        if area > A
            A = area
            A_tri = tri
        end
    end
    return A, A_tri
end

function pmesh(pv, hmax)
    """
    Generates triangular mesh of a polygon
    
    Arguments:
        pv   - closed polygon
        hmax - parameter to determine density of nodes
    """
    p = []
    
    for j = 1:length(pv)-1
        push!(p, pv[j])
        p1, p2 = pv[j], pv[j+1]
        append!(p, generate_nodes(p1, p2, hmax))           
    end
    push!(p, pv[end])
    
    t = nothing
    while true
        t = delaunay(p)
        
        remove_outside_triangles!(t, p, pv)
        
        A, A_tri = largest_area(t, p)
        
        if A > hmax^2/2
            circumcenter = tri_circumcenter(p[A_tri])
            if !(circumcenter in p)
                push!(p, circumcenter)
            end
        else
            break
        end
    end
    
    return p, t
end