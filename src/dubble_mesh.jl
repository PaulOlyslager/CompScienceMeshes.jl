using SparseArrays
using LinearAlgebra

mutable struct _VertexMapper
    list::Vector{Vector{Int}}
end
_VertexMapper(n::Integer) = _VertexMapper([Int[] for i in 1:n])
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
"""
    double_mesh(mesh::CompScienceMeshes.Mesh, startcell::Int, normal::Int)

doubles the mesh on the flat surfaces, keeps it the same on parts with volume. Can also
be used to generate continuous normal field. startcell is teh index of the cell from witch 
the algorithm starts, normal is a 1 or -1, depending on the normal on startcell.
"""
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
    a <= 0 && (a+=2*pi)
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
        #Select neighboorcell
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
                other_vectors = [other_vectors;[vertices[i]-cell_middle for i in oldcells[ind] if i ∉ edge]]
            end
            angles = []
            for vect in other_vectors
                x = dot(vect,edge_cell)
                y = dot(vect,cell_normal)
                push!(angles,angle(x,y))
            end
            ind = argmin(angles)
            neighboorcell_old_ind = e[2][ind] #index in old 
            neighboorcell_old = Vector(copy(oldcells[neighboorcell_old_ind]))
        end
        # Set orientation of neighboorcell correct
        if D[cell_old_ind,e[1]]*_permutation(cell_old,new_cell_old_notation) == D[neighboorcell_old_ind,e[1]] 
            neighboorcell_old[[1,2]] .= neighboorcell_old[[2,1]]
            s = 2
        else
            s = 1
        end
        # Select vertices for the neighboorcell
        neighboorcell_new = Tuple{Int,Int}[]
        for i in 1:3
            if neighboorcell_old[i] ∈ edge
                ind = findfirst(==(neighboorcell_old[i]),new_cell_old_notation)
                push!(neighboorcell_new,new_cell[ind])
            else
                push!(neighboorcell_new,(neighboorcell_old[i],_add_new_vertex!(vertexmapper,neighboorcell_old[i])))
            end
        end
        # Add neighboorcell or combine if it already exists
        existing_ind = normal_info[s][neighboorcell_old_ind]
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
"""
    find_neighboors(Connectivity matrix, face index)
yields a list with tuples: (edge index, [face indices])
"""
function find_neighboors(D,chart::Int)
    neighboors = Tuple{Int,Vector{Int}}[]
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

"""
    blowmesh(mesh,d)
moves for each face the vertices a distance d in the direction of the normal
"""
function blowmesh(mesh,distance)
    vert = copy(mesh.vertices)
    for (i,f) in enumerate(mesh.faces)
        vert[f] .= vert[f] .+ Ref(distance.*normal(chart(mesh,i)))
    end
    return Mesh(vert,mesh.faces)
end
"""
    remove_nonused_vertices(mesh)
removes the vertices that are not used in a face
"""
function remove_nonused_vertices(mesh)
    used = zeros(Bool,length(mesh.vertices))
    faces = typeof(mesh.faces[1])[]
    for face in mesh.faces
        used[face] .= true
    end
    for face in mesh.faces
        f = [sum(used[1:ind]) for ind in face]
        push!(faces,f)
    end
    Mesh(mesh.vertices[used],faces)

end
