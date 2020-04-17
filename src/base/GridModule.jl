module GridModule

#> "grid array" type to hold grid quantities
mutable struct Grid_Array
    number_grid_x::Int
    x_min::Float64
    x_max::Float64
    x::Array{Float64,1}
    number_grid_y::Int
    y_min::Float64
    y_max::Float64
    y::Array{Float64,1}
    number_grid_z::Int
    z_min::Float64
    z_max::Float64
    z::Array{Float64,1}
end

function GridArrayInit(
        number_grid_x,x_min,x_max,
        number_grid_y,y_min,y_max,
        number_grid_z,z_min,z_max
    )
    x = range(x_min,stop=x_max,length=number_grid_x)
    y = range(y_min,stop=y_max,length=number_grid_y)
    z = range(z_min,stop=z_max,length=number_grid_z)
    grid_array = Grid_Array(
        number_grid_x,x_min,x_max,x,
        number_grid_y,y_min,y_max,y,
        number_grid_z,z_min,z_max,z
    )
    return grid_array
end

end