
using StaticArrays
using LinearAlgebra
using CompScienceMeshes
using Test

const x̂ = point(1,0,0)

m = meshrectangle(1.0,1.0,0.5,3)
f = m.faces
f[1] = @SVector [1,2,5]
f[2] = @SVector [1,5,4]
f[7] = @SVector [5,9,8]
f[8] = @SVector [5,6,9]

ℑ₁ = Mesh(m.vertices,f)
CompScienceMeshes.mirrormesh
ℑ₂ = CompScienceMeshes.rotate(ℑ₁, 0.5π * x̂)
ℑ₃ = CompScienceMeshes.rotate(ℑ₁, 1.0π * x̂)
Γ₁ = weld(ℑ₁,-ℑ₂)
Γ₂ = weld(ℑ₂,-ℑ₃)
T = CompScienceMeshes.remove_nonused_vertices(weld(Γ₁,ℑ₃))

new_mesh = CompScienceMeshes.double_mesh(T,1,1)

@test numcells(new_mesh) == 48

e = skeleton(new_mesh,1)
D = connectivity(e,new_mesh)
for (i,face) in enumerate(new_mesh.faces)
    n = CompScienceMeshes.find_neighboors(D,i)
    @assert length(n)==3
    for k in n
        @assert length(k[2])==1 
        n1 = normal(chart(new_mesh,i))
        n2 = normal(chart(new_mesh,k[2][1]))
        t = sum([1 0;0 -1]*e.vertices[e.faces[k[1]]])
        c = sum(e.vertices[e.faces[k[1]]])/2
        c1 = sum(new_mesh.vertices[face])/3
        c2 = sum(new_mesh.vertices[new_mesh.faces[k[2][1]]])/3
        s1 = sign(det(hcat(c-c1,t,n1)))
        s2 = sign(det(hcat(c-c2,t,n2)))
        @assert s1==-s2
    end
end

