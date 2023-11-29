using CompScienceMeshes
using StaticArrays
using Combinatorics
using LinearAlgebra

function edges(v)
    return [simplex(v[3],v[1]);[simplex(v[i],v[i+1]) for i in 1:2]]
end
function onedge(p,e)
    v = verticeslist(e)
    p ≈ v[1] && return false
    p ≈ v[2] && return false
    sign(dot(v[1]-p,v[2]-p)) == -1 && return true
    return false
end
parallel(e1,e2) = abs(det(hcat(tangents(e1,1),tangents(e2,1)))) < eps()
function test_mesh(s1,s2,list)#s1 and s2 original similpices, list the list of simplices.
    ## volume test
    l = [list[1];list[2]]
    @assert sum(CompScienceMeshes.volume.([s1,s2])) ≈ sum(CompScienceMeshes.volume.(l))
    tol = sqrt(eps())
    ## find all edges
    edge = []
    for s in l
    push!(edge,edges(verticeslist(s))...)
    end

    for (e1,e2) in combinations(edge,2)
        if overlap(e1,e2)
            i = intersection(e1,e2)
            @assert length(i) == 1
            if volume(i[1]) > tol
                if !(sort(verticeslist(e1)) ≈ sort(verticeslist(e2)))
                    println((sort(verticeslist(e1)) , sort(verticeslist(e2))))
                end
                @assert sort(verticeslist(e1)) ≈ sort(verticeslist(e2))
            end
        elseif parallel(e1,e2)
            

        else
            pnew = CompScienceMeshes.intersectlines(verticeslist(e1)[1],verticeslist(e1)[2],verticeslist(e2)[1],verticeslist(e2)[2])
            @assert !(onedge(pnew,e1) && onedge(pnew,e2))
        end
    end

end

simp1 = simplex((@SVector [1.0,0.0]), (@SVector [0.0,1.0]),(@SVector [0.0,0.0]))

simps = [
simplex((@SVector [0.0,0.0]), (@SVector [1.0,0.0]), (@SVector [0.0,1.0])),
simplex((@SVector [0.0,0.5]), (@SVector [0.5,0.0]), (@SVector [0.5,0.5])),
simplex((@SVector [-0.1,0.5]), (@SVector [0.5,-0.1]), (@SVector [0.6,0.5])),
simplex((@SVector [1.0,1.0]), (@SVector [1.0,0.0]), (@SVector [0.0,1.0])),
simplex((@SVector [1.0,1.0]), (@SVector [2.0,-1.0]), (@SVector [0.0,1.0])),
simplex((@SVector [1.0,1.0]), (@SVector [1.0,0.0]), (@SVector [-1.0,2.0])),
simplex((@SVector [1.0,1.0]), (@SVector [2.0,-1.0]), (@SVector [-1.0,2.0]))]

@time out = CompScienceMeshes.complementary_mesh.(Ref(simp1),simps)

test_mesh.(Ref(simp1),simps,out)
