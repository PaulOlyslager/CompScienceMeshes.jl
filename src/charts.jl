#export Simplex
import Base.sign

# U: the dimension of the universe
# D: the dimension of the manifold
# N: the number of vertices
# T: the type of the coordinates
# C: the complimentary dimension (should always be U-D)
abstract type AbstractSimplex{U,D,C,N,T} end
struct Simplex{U,D,C,N,T} <: AbstractSimplex{U,D,C,N,T}
    vertices::SVector{N,SVector{U,T}}
    tangents::SVector{D,SVector{U,T}}
    normals::SVector{C,SVector{U,T}}
    volume::T
end
struct MirroredSimplex{U,D,C,N,T} <: AbstractSimplex{U,D,C,N,T}
    simplex::Simplex{U,D,C,N,T}
end


normal(t::Simplex{U,D,1}) where {U,D} = t.normals[1]
normal(t::MirroredSimplex) = -normal(t.simplex)
normals(t::Simplex) = t.normals
normals(t::MirroredSimplex) = -normals(t.simplex)

sign(s::Simplex{U,D,1}) where {U,D} = 1
sign(s::MirroredSimplex{<:Simplex{U,D,1}}) where {U,D} = -1

mirror(t::Simplex{U,D,1}) where {U,D} = MirroredSimplex(t)
mirror(t::MirroredSimplex{U,D,1}) where {U,D} = t.simplex

#normal(t::Simplex{3,2,1,3,<:Number}) = t.normals[1]
dimtype(splx::AbstractSimplex{U,D}) where {U,D} = Val{D}
"""
    permute_simplex(simplex,permutation)

Permutation is a Vector v which sets the v[i]-th vertex at the i-th place.

Return Simplex with permuted vertices list, tangents are recalculated, normal is kept the same
"""
function permute_vertices(s::AbstractSimplex{U,D,1,N,T},permutation::Union{Vector{P},SVector{P}}) where {U,D,P,N,T}
    vert = SVector{N,SVector{U,T}}(verticeslist(s)[permutation])
    simp = simplex(vert)
    sign(dot(normal(s),normal(simp))) == 1 && return simp
    return mirror(simp)
end

# """
#     flip_normal(simplex)

# Flips the normal of the simplex. Only on triangles embedded in 3D space
# """
# flip_normal(t::Simplex{3,2,1,3,<:Number}) = Simplex(t.vertices,t.tangents,-t.normals,t.volume)

"""
    flip_normal(simplex, sign)

Flips the normal of the simplex if sign is -1. Only on triangles embedded in 3D space
"""
function flip_normal(t::Simplex{3,2,1,3,<:Number},sign::Int)
    sign == 1 && return t
    return mirror(t)
    end
export flip_normal
"""
    boundary(simplex)

returns a list of simplices where the extra normal is defined to be pointing outward of the original simplex.
The basis defined by the tangents and normals is right-handed.
"""
function boundary(p::Simplex{U,D,C,N,T}) where {U,D,C,N,T}
    overall_normals = normals(p)
    v = verticeslist(p)
    s = [simplex(SVector((v[union(1:i-1,i+1:N)]...)),overall_normals) for i in 1:N]
    midles = [sum(v[union(1:i-1,i+1:N)])/(N-1) for i in 1:N]
    middle = sum(v)/N
    boundary = Simplex{U,D-1,C+1,N-1,T}[]
    
    for (mid,simp) in zip(midles,s)
        TV = typeof(verticeslist(simp))
        TT = typeof(simp.tangents)
        permutation = [i for i in 1:N-1]
        perm_tang = [i for i in 1:N-2]
        final_normals = SVector{U,T}[]
        TP = typeof(normals(simp))
        for n in normals(simp)
            if sign(dot(mid-middle,n)) != 0.0
                push!(final_normals,SVector{U,T}(n*sign(dot(mid-middle,n))))
            else
                push!(final_normals,n)
            end
        end
        m = hcat([simp.tangents;final_normals]...)
        if det(m) < -eps()
            permutation[1], permutation[2] = 2, 1
            perm_tang[1], perm_tang[2] = 2, 1
        end

        push!(boundary,typeof(simp)(TV(verticeslist(simp)[permutation]),TT(simp.tangents[perm_tang]),TP(final_normals),simp.volume))
    end
    return boundary
end

function boundary(p::Simplex{U,D,C,3,T}) where {U,D,C,T}
    overall_normals = normals(p)
    v = verticeslist(p)
    N=3
    s = [simplex(SVector((v[union(1:i-1,i+1:N)]...)),overall_normals) for i in 1:N]
    midles = [sum(v[union(1:i-1,i+1:N)])/(N-1) for i in 1:N]
    middle = sum(v)/N
    boundary = Simplex{U,D-1,C+1,N-1,T}[]
    
    for (mid,simp) in zip(midles,s)
        TV = typeof(verticeslist(simp))
        TT = typeof(simp.tangents)
        permutation = [i for i in 1:N-1]
        perm_tang = 1
        final_normals = SVector{U,T}[]
        TP = typeof(normals(simp))
        ns = Vector(normals(simp))
        ns[end] = SVector{U,T}(ns[end]*sign(dot(mid-middle,ns[end])))
        m = hcat([simp.tangents;ns]...)
        if det(m) < -eps()
            permutation[1], permutation[2] = 2, 1
            perm_tang = -1
        end
        push!(boundary,typeof(simp)(TV(verticeslist(simp)[permutation]),TT(simp.tangents*perm_tang),TP(final_normals),simp.volume))
    end
    return boundary
end

tangents(s::Simplex,i) = s.tangents[i]
tangents(s::MirroredSimplex,i) = tangents(s.simplex,i)
"""
    coordtype(simplex)

Return coordinate type used by simplex.
"""
coordtype(::Type{<:AbstractSimplex{U,D,C,N,T}}) where {U,D,C,N,T} = T
coordtype(p::AbstractSimplex) = coordtype(typeof(p))


"""
    volume(simplex)

Return the volume of the simplex.
"""
volume(p::Simplex) = p.volume
volume(p::MirroredSimplex) = volume(p.simplex)

"""
A tuple of points, aka an interval behaves trivially like a chart
"""
volume(x::Tuple{T,T}) where {T} = norm(x[2]-x[1])


"""
    dimension(simplex)

Return the manifold dimension of the simplex.
"""
dimension(::Type{<:AbstractSimplex{U,D,C,N,T}}) where {U,D,C,N,T} = D

"""
    dimension(simplex)

Return the manifold dimension of the simplex.
"""
dimension(p::AbstractSimplex) = dimension(typeof(p))


"""
    length(simplex)

Returns the number of vertices (equals dimension + 1)
"""
Base.length(p::AbstractSimplex) = dimension(typeof(p))+1


"""
  universedimension(p)

Return the dimension of the universe in which `p` is embedded.
"""
universedimension(::Type{<:AbstractSimplex{U,D,C,N,T}}) where {U,D,C,N,T} = U
universedimension(p::AbstractSimplex) = universedimension(typeof(p))


"""
    getindex(simplex, I)
    simplex[I]

Get the vertices at index I (scalar or array) defining the simplex
"""
getindex(p::AbstractSimplex, I::Union{Number,SVector,Array}) = verticeslist(p)[I]



"""
    simplex(vertices)
    simplex(v1, v2, ...)
    simplex(vertices, Val{D})

Build a D-dimensional simplex. The vertices can be passed in
an array (static or dynamic), or supplied separately. If the
length of the array is not part of its type, the speed of the
construction can be improved by supplying an extra Val{D}
argument. In case it is not clear from the context whether
the vertex array is dynamically or statically sized, use the
third form as it will not incur notable performance hits.

Note that D is the dimension of the simplex, i.e. the number
of vertices supplied minus one.
"""
@generated function simplex(vertices::SVector{D1,P}) where {D1,P}
    U = length(P)
    D = D1 - 1
    C = U-D
    T = eltype(P)
    xp1 =:(())
    for i in 1:D
        push!(xp1.args, :(vertices[$i]-vertices[end]))
    end
    xp2 = :(SVector{$D,P}($xp1))
    quote
        tangents = $xp2
        normals, volume = _normals(tangents, Val{$C})
        Simplex(vertices, tangents, normals, $T(volume))
    end
end

@generated function simplex(vertices::SVector{D1,P}, defined_normals::SVector{C,SVector{U,T}}) where {D1,P,C,U,T}
    @assert U == length(P)
    D = D1 - 1
    @assert C <= U-D
    CC = U-D
    @assert T == eltype(P)
    xp1 =:(())
    for i in 1:D
        push!(xp1.args, :(vertices[$i]-vertices[end]))
    end
    xp2 = :(SVector{$D,P}($xp1))
    quote

        tangents = $xp2
        normals, volume = _normals(tangents, defined_normals, Val{$CC})
        Simplex(vertices, tangents, normals, $T(volume))
    end
end

simplex(vertices...) = simplex(SVector((vertices...,)))

@generated function simplex(vertices, ::Type{Val{D}}) where D
    P = eltype(vertices)
    D1 = D + 1
    xp = :(())
    for i in 1:D1
        push!(xp.args, :(vertices[$i]))
    end
    :(simplex(SVector{$D1,$P}($xp)))
end

# normal(s::Simplex{3,2,1,3,T}) where {T} = normalize(cross(s[1]-s[3],s[2]-s[3]))


function _normals(tangents, ::Type{Val{1}})
    PT = eltype(tangents)
    D  = length(tangents)
    T  = eltype(PT)

    n = zeros(T,D+1)
    b = Array{T}(undef,D,D)

    for i in 1:D+1
        fill!(b, zero(T))
        for j in 1:D
            for k in 1:i-1
                b[k,j] = tangents[j][k]
            end
            for k in i:D
                b[k,j] = tangents[j][k+1]
            end
        end
        n[i] = (-1)^(i-1) * det(b)
    end

    n *= (-1)^D / norm(n)
    normals = SVector{1,PT}([PT(n)])

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    volume = sqrt(abs(det(metric))) /  factorial(D)

    return normals, volume

end

function _normals(tangents::SVector{2,SVector{3,T}}, ::Type{Val{1}}) where {T}

    n = tangents[1] × tangents[2]
    l = norm(n)

    P = SVector{3,T}
    SVector{1,P}(n/l), 0.5*l
end




function _normals(tangents, ::Type{Val{C}}) where C
    PT = eltype(tangents)
    D  = length(tangents)
    U = length(PT)
    T  = eltype(PT)

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    volume = sqrt(abs(det(metric))) / factorial(D)
    t = hcat(tangents...)
    length(t) != 0 ? N = Matrix{T}(I,U,U)-t*(metric^-1)*transpose(t) : N = Matrix{T}(I,U,U)
    # Fix this. This function needs to become gneerated
    #N = Matrix{T}(I,U,U)-t*(metric^-1)*transpose(t)
    n = PT[]

    for i in eachrow(qr(N).R)
        norm(i)>0.8 && push!(n,PT(i))
    end

    if C > 0
        det(hcat([tangents;n]...)) < -eps() && (n[end] = -n[end])
    end
   
    normals = SVector{C,PT}(PT[n[i] for i in 1:C])

    return normals, volume
end
function _normals(tangents, defined_normals, ::Type{Val{C}}) where C #TODO fix this function
    PT = eltype(tangents)
    D  = length(tangents)
    U = length(PT)
    T  = eltype(PT)
    

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    ndef = hcat(defined_normals...)
    length(ndef)==0 && (ndef=0)
    volume = sqrt(abs(det(metric))) / factorial(D)
    t = hcat(tangents...)
    # Fix this. This function needs to become gneerated
    N = I-t*(metric^-1)*transpose(t).-ndef*transpose(ndef)
    n = PT[]
    for i in defined_normals; push!(n,i); end
    
    for i in eachrow(qr(transpose(N)).R)
        norm(i)>0.8 && push!(n,PT(Vector(i)))
    end
    det(hcat([tangents;n]...)) < -eps() && (n[end] = -n[end])
    normals = SVector{C,PT}((n)...)

    return normals, volume
end


function _normals(tangents::SVector{1,SVector{3,T}} where {T}, ::Type{Val{2}})
    P = eltype(tangents)
    normals = SVector{2,P}(zero(P), zero(P))
    # volume = dot(cross(tangents[1], tangents[2]), tangents[3]) / 6
    volume = norm(tangents[1])
    return normals, volume
end

function _normals(tangents::SVector{3,SVector{3,T}} where {T}, ::Type{Val{0}})
    P = eltype(tangents)
    normals = SVector{0,P}()
    volume = dot(cross(tangents[1], tangents[2]), tangents[3]) / 6
    return normals, volume
end



"""
    barytocart(simplex, uv)

Returns the point in the simplex with barycentric coordinates uv
"""
function barytocart(mani::AbstractSimplex, u)
    r = last(verticeslist(mani))
    for i in 1 : dimension(mani)
        ti = tangents(mani,i)
        ui = u[i]
        #r += mani.tangents[i] * u[i]
        r += ti * ui
    end
    return r
end


"""
    carttobary(simplex, point) -> barycoords

Compute the barycentric coordinates on 'simplex' of 'point'.
"""
function carttobary(p::Simplex{U,D,C,N,T}, cart) where {U,D,C,N,T}

    G = [dot(tangents(p,i), tangents(p,j)) for i in 1:D, j in 1:D]
    #w = [dot(p.tangents[i], cart - p.vertices[end]) for i in 1:D]

    o = verticeslist(p)[end]
    w = [sum(t[j]*(cart[j]-o[j]) for j in 1:length(cart)) for t in tangentlist(p)]

    u = G \ w

    return SVector{D}(u)
end

const _combs24 = [SVector{2}(c) for c in combinations(@SVector[1,2,3,4],2)]
const _edgeidx24 = [relorientation(c, SVector(1,2,3,4)) for c in _combs24]
const _combs24_pos = [(p > 0 ? _combs24[i] : reverse(_combs24[i])) for (i,p) in enumerate(_edgeidx24)]
const _edgeidx24_abs = abs.(_edgeidx24)

function edges(s::AbstractSimplex{3,3})
    T = eltype(eltype(verticeslist(s)))
    P = Simplex{3,1,2,2,T}
    Edges = Vector{P}(undef, length(_combs24_pos))
    for (i,c) in zip(_edgeidx24_abs, _combs24_pos)
        edge = simplex(verticeslist(s)[c[1]], verticeslist(s)[c[2]])
        Edges[i] = edge
    end
    return Edges
end

# function edges(s::Simplex{3,3})
#     C = [SVector{2}(c) for c in combinations(@SVector[1,2,3,4],2)]
#     T = eltype(eltype(s.vertices))
#     P = Simplex{3,1,2,2,T}
#     Edges = Vector{P}(undef, length(C))
#     for c in C
#         q = relorientation(c, @SVector[1,2,3,4])
#         Q = abs(q)
#         Edges[Q] = q > 0 ?
#             simplex(s.vertices[c[1]], s.vertices[c[2]]) :
#             simplex(s.vertices[c[2]], s.vertices[c[1]])
#     end
#     return Edges
# end

function edges(s::AbstractSimplex{3,2})
    return [
        simplex(verticeslist(s)[2], verticeslist(s)[3]),
        simplex(verticeslist(s)[3], verticeslist(s)[1]),
        simplex(verticeslist(s)[1], verticeslist(s)[2])]
end

function faces(c)
    @SVector[
        simplex(c[2], c[3]),
        simplex(c[3], c[1]),
        simplex(c[1], c[2]),
    ]
end

function faces(c::AbstractSimplex{3,3,0,4,T}) where {T}
    @SVector[
        simplex(c[4],c[3],c[2]),
        simplex(c[1],c[3],c[4]),
        simplex(c[1],c[4],c[2]),
        simplex(c[1],c[2],c[3])
    ]
end


"""
    ReferenceSimplex{Dimension, CoordType, NumVertices}

This domain is defined to bootstrap the quadrature generation strategy. The generic definition
of numquads on a chart pulls back to the domain. For a limit set of reference domains, explicit
quadrature rules are defined. The weights and points are then pushed forward to the configuaration
space element over which integration is desired.

For more details see the implementation in quadpoints.jl
"""
struct ReferenceSimplex{D,T,N}
    # N: number of defining points
    # D: dimension of the simplex
    # T: type of the coordinates
    simplex::Simplex{D,D,0,N,T}
end

function ReferenceSimplex{D,T,N}() where {D,T,N}
    P = SVector{D,T}[]
    for i in 1:D
        a = zeros(T,D)
        a[i] = 1
        p = SVector{D,T}(a)
        push!(P,p)
    end

    o = zero(SVector{D,T})
    push!(P,o)
    ReferenceSimplex{D,T,N}(simplex(P...))
end


barytocart(ch::ReferenceSimplex, u) = barytocart(ch.simplex, u)
carttobary(ch::ReferenceSimplex, p) = carttobary(ch.simplex, p)

domain(ch::AbstractSimplex{U,D,C,N,T}) where {U,D,C,T,N} = ReferenceSimplex{D,T,N}()
neighborhood(ch::ReferenceSimplex, u) = u

"""
    tangents(splx)

Returns a matrix whose columns are the tangents of the simplex `splx`.
"""
tangents(splx::Simplex) = hcat((splx.tangents)...)
tangents(splx::MirroredSimplex) = tangents(splx.simplex)

vertices(splx::Simplex) = hcat((splx.vertices)...)
vertices(splx::MirroredSimplex) = vertices(splx.simplex)

"""
    verticeslist(simplex)

Returns the vertices as a list.
"""
verticeslist(splx::Simplex) = splx.vertices
verticeslist(splx::MirroredSimplex) = verticeslist(splx.simplex)

tangentlist(splx::Simplex) = splx.tangents
tangentlist(splx::AbstractSimplex) = tangents(splx.simplex)
export verticeslist
export tangentlist