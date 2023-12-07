"""
    intersect(chartA, chartB) -> [chart1, chart2, ...]

Compute the intersection of two charts of equal dimension.

Compute an array of charts such that the disjoint union of those charts produces the intersection of the two charts provided as inputs. In particular the sum of integrals over the returned charts will equal the integral over the intersection of the two given charts.
"""
function intersection(p1::Simplex{U,1,C,2,T}, p2::Simplex{U,1,C,2,T}) where {U,C,T}

    tol = sqrt(eps(T))
    PT = eltype(p1.vertices)

    W = [u for u in p2.vertices]
    clipconvex!(W, p1.vertices[1],  p1.tangents[1])
    clipconvex!(W, p1.vertices[2], -p1.tangents[1])

    # For consistency an array needs to be return. In higher
    # dimensions the intersection could be the union of
    # multiple simplices.
    [simplex(W[1], W[2])]
end

function clipconvex!(W, v, m)
    m2 = dot(m,m)
    for i in 1:length(W)
        s = dot(W[i]-v,m)
        s = min(0, s)
        W[i] = v + s / m2 * m
    end
end

"""
    intersection(triangleA, triangle B)

The output inherits the orientation from the first operand
"""
function intersection(p1::Simplex{U,2,C,3}, p2::Simplex{U,2,C,3}) where {U,C}
  pq = sutherlandhodgman(p1.vertices, p2.vertices)
  return [ simplex(pq[1], pq[i], pq[i+1]) for i in 2:length(pq)-1 ]
end

function intersection(p1::Simplex{3,2,1,3}, p2::Simplex{3,2,1,3};tol=eps()) 
    pq = sutherlandhodgman(p1.vertices, p2.vertices)
    nonoriented_simplexes = [ simplex(pq[1], pq[i], pq[i+1]) for i in 2:length(pq)-1 ]
    nonoriented_simplexes = nonoriented_simplexes[volume.(nonoriented_simplexes).>tol]
    signs = Int.(sign.(dot.(normal.(nonoriented_simplexes),Ref(normal(p1)))))
    flip_normal.(nonoriented_simplexes,signs)
  end


function intersection(p1::Simplex{U,3,C,4}, p2::Simplex{U,3,C,4}) where {U,C}
  @assert overlap(p1,p2)
  volume(p1) <= volume(p2) ? [p1] : [p2]
end

"""
    intersection2(triangleA, triangleB)

returns intersection in which the first operand inhirits orientation from first argument
and second argument inhirits orientation of second argument.
"""
function intersection2(p1::Simplex, p2::Simplex)
    a = intersection(p1,p2)
    return [(i,i) for i in a]
end

function intersection2(p1::Simplex{3,2,1,3}, p2::Simplex{3,2,1,3}; tol=eps())
    pq = sutherlandhodgman(p1.vertices, p2.vertices)
    nonoriented_simplexes = [ simplex(pq[1], pq[i], pq[i+1]) for i in 2:length(pq)-1 ]
    nonoriented_simplexes = nonoriented_simplexes[volume.(nonoriented_simplexes).>tol]
    signs1 = Int.(sign.(dot.(normal.(nonoriented_simplexes),Ref(normal(p1)))))
    
    signs2 = Int.(sign.(dot.(normal.(nonoriented_simplexes),Ref(normal(p2)))))
    return [(flip_normal(s,signs1[i]),flip_normal(s,signs2[i])) for (i,s) in enumerate(nonoriented_simplexes)]
end
export intersection2

"""
    intersectline(a,b,p,q)

Computes the intersection of the lines (in a 2D space) defined
by points [a,b] and [p,q]
"""
function intersectlines(a,b,p,q)

    P = typeof(a)

    d1 = det(@SMatrix [a[1] a[2]; b[1] b[2]])
    d2 = det(@SMatrix [p[1] p[2]; q[1] q[2]])

    id = one(eltype(a))
    d1x = det(@SMatrix [a[1] id; b[1] id])
    d1y = det(@SMatrix [a[2] id; b[2] id])

    d2x = det(@SMatrix [p[1] id; q[1] id])
    d2y = det(@SMatrix [p[2] id; q[2] id])

    den = det(@SMatrix [d1x d1y; d2x d2y])
    @assert !isinf(1/den)

    P(
        det(@SMatrix [d1 d1x; d2 d2x]) / den,
        det(@SMatrix [d1 d1y; d2 d2y]) / den,
    )

end



"""
    inside(p,a,b)

Determines is p is on the interior side of the segment b of the boundary,
assuming that the boundary is oriented counter-clockwise.
"""
function leftof(p,a,b)

    tol = sqrt(eps(eltype(a)))
    ap = @SVector [ p[1]-a[1], p[2]-a[2] ]
    ab = @SVector [ b[1]-a[1], b[2]-a[2] ]
    ap[1]*ab[2] - ap[2]*ab[1] <= tol ? true : false

end

function _leftof_extended(p,a,b)
    tol = sqrt(eps(eltype(a)))
    ap = @SVector [ p[1]-a[1], p[2]-a[2] ]
    ab = @SVector [ b[1]-a[1], b[2]-a[2] ]
    det = ap[1]*ab[2] - ap[2]*ab[1] 
    det <= -tol && return 1
    det >= tol && return -1
    return 0
end

#export sutherlandhodgman2d

"""
    sutherlandhodgman2d(subject,clipper)

Computes the intersection of the coplanar triangles
subject and clipper.
"""
function sutherlandhodgman2d(subject,clipper)

    PT = eltype(subject)

    clipped = Array(copy(subject))
    sizehint!(clipped, 8)

    input = Array(copy(clipped))
    sizehint!(input, 8)

    b = last(clipper)
    for a in clipper

        resize!(input, length(clipped))
        copyto!(input, clipped)
        resize!(clipped, 0)

        isempty(input) || (q = last(input))

        for p in input
            if leftof(p, b, a)
                if !leftof(q, b, a)
                    ist = intersectlines(b,a,q,p)
                    push!(clipped, ist)
                end

                push!(clipped, p)

            elseif leftof(q, b, a)

                ist = intersectlines(b,a,q,p)
                push!(clipped, ist)
            end

            q = p
        end
        b = a
    end

    return clipped
end

function sutherlandhodgman2d_full(subject,clipper)
    
    PT = eltype(subject)

    clipped = Array(copy(subject))
    sizehint!(clipped, 8)

    input = Array(copy(clipped))
    sizehint!(input, 8)

    b = last(clipper)
    other = Array(copy(subject))
    sizehint!(other,8)
    set_of_sets = []
    
    for a in clipper
        resize!(other,0)
        resize!(input, length(clipped))
        copyto!(input, clipped)
        resize!(clipped, 0)

        isempty(input) || (q = last(input))

        for p in input
            if _leftof_extended(p, b, a) == 1
                if _leftof_extended(q, b, a) == -1
                    ist = intersectlines(b,a,q,p)
                    push!(clipped, ist)
                    push!(other,ist)
                end

                push!(clipped, p)
            elseif _leftof_extended(p, b, a) == 0
                push!(clipped,p)
                push!(other,p)

            elseif _leftof_extended(q, b, a) == 1
                push!(other,p)
                ist = intersectlines(b,a,q,p)
                push!(clipped, ist)
                push!(other,ist)
            else
                push!(other,p)
            end

            q = p
        end

        
        t = [other[i]-other[end] for i in 1:length(other)-1]
        if rank(hcat((t...))) < 2 
            b = a
            continue
        end

        ind = findall(i->rank(hcat(([a-b,other[i]-b]...)))==1,1:length(other))
        if length(other[ind]) > 1
            p1,p2 = other[ind]
            sign(dot(p1-a,p2-a)) == -1 && (push!(other, a))
            sign(dot(p1-b,p2-b)) == -1 && (push!(other, b))
        end
        push!(set_of_sets,unique(order_counter_clockwise(other)))
        b=a
        
    end
    
    return set_of_sets,unique_tol(clipped)
end
function order_counter_clockwise(vertices) #should be convex
    middle = sum(vertices)/length(vertices)
    transl = vertices .- Ref(middle)
    angles =  [atan(t[2],t[1]) for t in transl]
    return [i[2] for i in sort(zip(angles,vertices))]
end
#export sutherlandhodgman

"""
    sutherlandhodgman(subject, clipper)

Compute the intersection of two coplanar triangles, potentially
embedded in a higher dimensional space.
"""
function sutherlandhodgman(subject, clipper)

    triangle = simplex(clipper, Val{2})
    subject2d = [carttobary(triangle,p) for p in subject]
    for p in subject2d for x in p @assert !isinf(x) end end
    clipper2d = [carttobary(triangle,q) for q in clipper]
    for p in clipper2d for x in p @assert !isinf(x) end end
    clipped2d = sutherlandhodgman2d(subject2d, clipper2d)
    clipped = [barytocart(triangle,q) for q in clipped2d ]

end
function sutherlandhodgman_full(subject, clipper)

    triangle = simplex(clipper, Val{2})
    subject2d = [carttobary(triangle,p) for p in subject]
    for p in subject2d for x in p @assert !isinf(x) end end
    clipper2d = [carttobary(triangle,q) for q in clipper]
    for p in clipper2d for x in p @assert !isinf(x) end end
    set_of_sets2d,clipped2d = sutherlandhodgman2d_full(subject2d, clipper2d)
    set_of_sets = [[barytocart(triangle,i) for i in q] for q in set_of_sets2d]
    clipped = [barytocart(triangle,q) for q in clipped2d ]
    return set_of_sets,clipped
end

function complementary_mesh(p1::Simplex{U,2,C,3}, p2::Simplex{U,2,C,3}; tol=eps()) where {U,C}
    set_of_sets1,pq = sutherlandhodgman_full(p1.vertices, p2.vertices)
    set_of_sets2,_ = sutherlandhodgman_full(p2.vertices, p1.vertices)

    nonoriented_pq = [ simplex(pq[1], pq[i], pq[i+1]) for i in 2:length(pq)-1 ]
    nonoriented_pq = nonoriented_pq[volume.(nonoriented_pq).>tol]

    nonoriented_set2 = vcat(create_triangles.(set_of_sets2)...)
    nonoriented_set1 = vcat(create_triangles.(set_of_sets1)...)

    signs1_pq = Int.(sign.(dot.(normal.(nonoriented_pq),Ref(normal(p1)))))
    signs2_pq = Int.(sign.(dot.(normal.(nonoriented_pq),Ref(normal(p2)))))
    signs_set1 = Int.(sign.(dot.(normal.(nonoriented_set1),Ref(normal(p1)))))
    signs_set2 = Int.(sign.(dot.(normal.(nonoriented_set2),Ref(normal(p2)))))

    out1 = [flip_normal.(nonoriented_pq,signs1_pq);flip_normal.(nonoriented_set1,signs_set1)]
    out2 = [flip_normal.(nonoriented_pq,signs2_pq);flip_normal.(nonoriented_set2,signs_set2)]
    return out1,out2

end

function create_triangles(p)
    p = unique_tol(p)
    set = []
    l = length(p)
    for j in 0:l
        f = true
        resize!(set,0)
        for i in 2:l-1
            f *= rank(hcat(p[mod1(1+j,l)]-p[mod1(i+j,l)],p[mod1(i+1+j,l)]-p[mod1(i+j,l)])) == 2
            f && (push!(set,simplex(p[mod1(1+j,l)],p[mod1(i+j,l)],p[mod1(i+1+j,l)])))
        end
        if f
            return set
        end
    end
    println(p)
    return nothing
end

# function create_triangles(p)
    
#     l = length(p)
#     l==3 && return [simplex(p[1],p[2],p[3])]
#     mid = sum(p)/length(p)
#     set = [simplex(mid,p[end],p[1])]
    
#     for i in 1:l-1
#         (push!(set,simplex(mid,p[i],p[i+1])))
#     end 
#     return set
# end


function complementary_mesh(p1::Simplex{U,2,0,3}, p2::Simplex{U,2,0,3}; tol=eps()) where {U}
    set_of_sets1,pq = sutherlandhodgman_full(p1.vertices, p2.vertices)
    set_of_sets2,_ = sutherlandhodgman_full(p2.vertices, p1.vertices)

    nonoriented_pq = [ simplex(pq[1], pq[i], pq[i+1]) for i in 2:length(pq)-1 ]
    nonoriented_pq = nonoriented_pq[volume.(nonoriented_pq).>tol]

    nonoriented_set2 = vcat(create_triangles.(set_of_sets2)...)
    nonoriented_set1 = vcat(create_triangles.(set_of_sets1)...)

    out1 = [nonoriented_pq; nonoriented_set1]
    out2 = [nonoriented_pq; nonoriented_set2]
    return out1,out2
    
end


u = simplex((@SVector [0.0,0.0,0.0]),(@SVector [1.0,0.0,0.0]),(@SVector [0.0,1.0,0.0]))
v = simplex((@SVector [0.0,0.0,0.0]),(@SVector [-1.0,0.0,0.0]),(@SVector [0.0,1.5,0.0]))