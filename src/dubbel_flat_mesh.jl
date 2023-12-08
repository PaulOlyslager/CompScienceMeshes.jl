#stap 1: geef mesh + normaal op 1 cell
#stap 2: loop over structuur en maak charts aan met: punten, normaal, maak lijst met [1,-1] voor alle faces
# stap 2 --> recursief, wanneer al op bestaande chart --> stop.

#stap 3 : voeg voor elke chart de hoekpunten toe, 
#stap  : neem elke chart en vouw over elke edge en kijk welke chart je bereikt, stitch vertices van gemeenschappelijke edge
using CompScienceMeshes
using SparseArrays
struct Chart
    vertices # the real ones
    normal::Int
    oiginal_face::Int
end

mutable struct _VertexMapper
    list::Vector{Vector{Int}}
end
_VertexMapper(n::Int) = _VertexMapper([Int[] for i in 1:n])
function _add_new_vertex!(v::_VertexMapper,index)
    refs = v.list[index]
    push!(v.list[index],length(refs)+1)
    return length(refs)
end
function _combine!(v::_VertexMapper,index,tup)
    refs = v.list[index]
    ref1, ref2 = refs[tup]
    v.list[index][refs.==ref2] .= Ref(ref1)
    return nothing
end
#### IMPORTANT: CHANGE FIRST 2 VERTICES OF A FACE TO CHANGE NORMAL
#double_mesh(mesh::CompScienceMeshes.AbstractMesh, startcell::Int, normal::Int) = nothing
function double_mesh(mesh::CompScienceMeshes.Mesh, startcell::Int, normal::Int)
    @warn "Mesh cannot have edges with 2 sides at the boundary while edge is not on boundary"
    edges = skeleton(mesh,1)
    vertexmapper = _VertexMapper(length(mesh.vertices))
    original_normal = zeros(Int,numcells(mesh)) #if zero chart does not exist, otherwise it gives index of the chart.
    oposite_normal = zeros(Int,numcells(mesh))
    normal_info = (original_normal,oposite_normal)

    cell = []
    for i in mesh.faces[startcell]
        push!(cell,(i,_add_new_vertex!(vertexmapper,i)))
    end
    if normal == -1
        n = 2
        cell[[1,2]] .= cell[[2,1]]
    elseif normal == 1
        n = 1
    else
        @error "normal was not 1 or -1"
    end
    normal_info[n][startcell] = 1
    newcells_ind = Int[startcell] #given new cell gives old cell
    newcells = Vector{Tuple{Int,Int}}[cell] #list of newcells, (index vertex,index vertexmapper)
    


    D = connectivity(edges,mesh)
    display(D)
    _double_mesh!(1,mesh.faces,normal_info,D,vertexmapper,newcells,newcells_ind,edges.faces,mesh.vertices)
    return _generate_mesh(mesh,vertexmapper,newcells)
end
#cell: index in newcell array
#oldcells: mesh.faces van original mesh
#nomral_info: tuple with folowing info:
    #orignal normal: contains newcell indices of the faces with this normal
    #oposite_normal: contains newcell indices of the faces with oposite normal

#D : connectivitye matrix oldcells, oldedges
#vertexmapper: techniek to double vertices
#newcells: array with all new cells: [(v1_old,index_vertex_mapper),(v2,...),(v3,...)]
#newcells_ind: given new cell index, this yields the old cell index
#edges: might not be nessecary but list of edges.
function _permutation(cell1,cell2) #permutation in first 2 indices!!!
    cell1[1:2]==cell2[1:2] && return 1
    return -1
end
function angle(x,y)
    a = atan(y,x)
    a < 0 && (a+=2*pi)
    return a
end

function _double_mesh!(cell,oldcells,normal_info,
    D, vertexmapper::_VertexMapper,newcells,newcells_ind,edges,vertices)

    
    cell_old_ind = newcells_ind[cell]
    cell_old = oldcells[cell_old_ind]
    new_cell = newcells[cell]
    new_cell_old_notation = [i[1] for i in new_cell]
    neighboors = find_neighboors(D,cell_old_ind)

    for e in neighboors
        edge = edges[e[1]]
        if length(e[2])==1
            neighboorcell_old_ind = e[2][1] #index in old 
            neighboorcell_old = Vector(copy(oldcells[neighboorcell_old_ind]))
            
        else
            cell_normal = normal(simplex(vertices[cell_old]...))*_permutation(cell_old,new_cell_old_notation)
            edge_middle = sum(vertices[edge])/2
            cell_middle = sum(vertices[cell_old])/3
            edge_tangent = normalize(tangents(simplex(vertices[edge]...),1))
            a = cell_middle-edge_middle
            edge_cell = normalize(a-edge_tangent*dot(a,edge_tangent))
            
            other_vectors = []
            for ind in e[2]# voegt de vertices toe die niet op de edge liggen
                other_vectors = [other_vertices;[vertices[i]-cell_middle for i in oldcells[ind] if i ∉ edge]]
            end
            angles = []
            for vect in other_vectors
                x = dot(vect,edge_cell)
                y = dot(vect,cell_normal)
                push!(angles,angle(x,y))
            end
            ind = argmin(angles)
            @warn "not tested yet"
            neighboorcell_old_ind = e[2][ind] #index in old 
            neighboorcell_old = Vector(copy(oldcells[neighboorcell_old_ind]))
        end

        if D[cell_old_ind,e[1]]*_permutation(cell_old,new_cell_old_notation) == D[neighboorcell_old_ind,e[1]] 
            neighboorcell_old[[1,2]] .= neighboorcell_old[[2,1]]
            s = 2
        else
            s = 1
        end
        neighboorcell_new = Tuple{Int,Int}[]
        for i in 1:3
            if neighboorcell_old[i] ∈ edge
                ind = findfirst(==(neighboorcell_old[i]),new_cell_old_notation)
                push!(neighboorcell_new,new_cell[ind])
            else
                push!(neighboorcell_new,(neighboorcell_old[i],_add_new_vertex!(vertexmapper,neighboorcell_old[i])))
            end
        end

        existing_ind = normal_info[s][neighboorcell_old_ind]
        ### stap 1: kijken of deze al bestaat, zo ja paste vertices samen.
        if existing_ind !=0
            existing = newcells[existing_ind]
            for i in 1:3
                @assert existing[i][1] == neighboorcell_new[i][1]
                _combine!(vertexmapper,existing[i][1],[existing[i][2],neighboorcell_new[i][2]])
            end
        else
            ### toevoegen aan newcells
            ### toevoegen aan new_cell_ind
            ### toevoegen aan normaal array
            push!(newcells,neighboorcell_new)
            push!(newcells_ind, neighboorcell_old_ind)
            new_cell_ind = length(newcells_ind)
            normal_info[s][neighboorcell_old_ind] = new_cell_ind
            _double_mesh!(new_cell_ind,oldcells, normal_info, D,
            vertexmapper,newcells,newcells_ind, edges, vertices)
        end
        
    end
    return nothing
end


function _generate_mesh(oldmesh::CompScienceMeshes.Mesh{U,D1,T}, vertexmapper::_VertexMapper,cells) where {U, D1,T}
    new_vertices = SVector{U,T}[]
    for (i,v) in enumerate(oldmesh.vertices)
        numbers = unique(vertexmapper.list[i])
        n = length(numbers)
        new_numbers = [i for i in 1:n] .+ length(new_vertices)
        vertexmapper.list[i] = replace(vertexmapper.list[i], (numbers .=> new_numbers)...)
        for _ in 1:n ; push!(new_vertices,v) end
    end
    new_faces = SVector{D1,Int}[]
    for cell in cells
        inds = []
        for (v,vm) in cell
            push!(inds,vertexmapper.list[v][vm])
        end
        push!(new_faces,SVector{3,Int}(inds))
    end
    return CompScienceMeshes.Mesh(new_vertices,new_faces)
end

function find_neighboors(D,chart)
    neighboors = []
    rows = rowvals(D)
    m,n = size(D)
    for edge in 1:n
        r = rows[nzrange(D,edge)]
        if chart∈r
            if length(r)==1
                push!(neighboors,(edge,r))
            else
                push!(neighboors,(edge,[i for i in r if i!=chart]))
            end
        end
    end
    return neighboors
end



# geef mesh + normaal op 1 cell
# loop recursief over structuur, telken een face kantelen over de 3 edges, zien waar je uitkomt.
    # optie1: nieuw hoekpunt nog nergens gebruikt, gebruik dit dan voor de nieuwe face
    # optie2: nieuw hoekpunt al ergens gebruikt, maak een kopie aan
    # optie3: kantelen op face die al bestaat -> stop en stitch common edge.
#
#
function test!(out,i,l)
    print(length(out))
    println(l)
    push!(out,i)
    if i<1
        return 0
    end
    return [test!(out,i-1,"a"), test!(out,i-1,"b"), test!(out,i-1,"c")]
end

using StaticArrays

mesh = meshrectangle(1.0,1.0,0.5,3)
f = mesh.faces
f[1] = @SVector [1,2,5]
f[2] = @SVector [1,5,4]
f[7] = @SVector [5,9,8]
f[8] = @SVector [5,6,9]
mesh = Mesh(mesh.vertices,f)
new_mesh = double_mesh(mesh,1,1)
import PlotlyJS
display(PlotlyJS.plot(normals(new_mesh)))