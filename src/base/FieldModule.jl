module FieldModule

using FFTW
using GridModule
using LinearAlgebra

#> Field struct
mutable struct Field_Array
    phi::Array{Float64,1}
    ion_density::Array{Float64,1}
    electron_density::Array{Float64,1}
end

# outer constructor that zeroes fields
function Fields_Init(grid) 
    fields = Field_Array(
        fill(0.0,length(grid.x)),
        fill(0.0,length(grid.x)),
        fill(0.0,length(grid.x))
    )
    return fields
end

function SolvePoissonEquation_FFT_1D(
        RHS_array::Array{Float64,1},
        grid_array::GridModule.Grid_Array
    )
    #> Put the right-hand-side quantity into k-space
    RHS_fft=FFTW.rfft(RHS_array[1:end-1])
    
    #> Find the wavenumber array in k-space
    k_fft=range(0,stop=length(RHS_fft)-1,length=length(RHS_fft))
    k_fft*=2.0*pi/(grid_array.x[end]-grid_array.x[1])
    
    #> k^2 (e*phi/T) = ni - ne
    field_fft_out = RHS_fft./(k_fft.^2)
    
    #> For quasi-neutrality, k=0 is dropped.
    field_fft_out[1] = 0.0+0.0im
    #field_fft_out[3:end].*=0.0+0.0im
    
    #> Pull back into real-space
    field_out=FFTW.irfft(field_fft_out,length(RHS_array[1:end-1]))
    push!(field_out,field_out[1])
    
    return field_out
end

function CalculateElectricField_FFT_1D(
        phi::Array{Float64,1},
        grid_array::GridModule.Grid_Array
    )
    #> Put the right-hand-side quantity into k-space
    RHS_fft=FFTW.rfft(phi[1:end-1])
    
    #> Find the wavenumber array in k-space
    k_fft=range(0,stop=length(RHS_fft)-1,length=length(RHS_fft))
    k_fft*=2.0*pi/(grid_array.x[end]-grid_array.x[1])
    
    #> -i k (phi) = E
    field_fft_out = -1.0im*RHS_fft.*(k_fft)
    #field_fft_out[3:end].*=0.0+0.0im
    #field_fft_out[1]=0.0+0.0im
    
    #> Pull back into real-space
    field_out=FFTW.irfft(field_fft_out,length(phi[1:end-1]))
    push!(field_out,field_out[1])
    
    return field_out
end

function CalculateElectricField_FD_1D(
        phi::Array{Float64,1},
        grid_array::GridModule.Grid_Array
    )
    field_out = fill(0.0,length(phi))
    dx = grid_array.x[2]-grid_array.x[1]
    
    for i in 2:length(phi)-1
        field_out[i] = -(phi[i+1]-phi[i-1])/(2.0*dx)
    end
    field_out[1]=-(phi[2]-phi[1])/dx
    field_out[end]=-(phi[end]-phi[end-1])/dx
    return field_out
end

function CalculateElectricField_FD_1D_conducting_boundaries(
        phi::Array{Float64,1},
        grid_array::GridModule.Grid_Array
    )
    field_out = fill(0.0,length(phi))
    dx = grid_array.x[2]-grid_array.x[1]
    
    for i in 2:length(phi)-1
        field_out[i] = -(phi[i+1]-phi[i-1])/(2.0*dx)
    end
    field_out[1]=0.0
    field_out[end]=0.0
    return field_out
end


function ApplyPeriodicBoundary!(
        field_x::Array{Float64,1}
    )
    field_x[1] += field_x[end]
    field_x[end] = field_x[1]
end

#> "grid array" type to hold grid quantities
mutable struct Poisson_1D
    A::LinearAlgebra.SymTridiagonal
    phi_boundary_left::Float64
    phi_boundary_right::Float64
    RHS_boundary_left::Float64
    RHS_boundary_right::Float64    
end

function Poisson_1D_Init(
        grid_x,
        phi_bound_left,phi_bound_right
    )
    dx = grid_x[2]-grid_x[1]
    #> Note that the default length should NOT 
    #> be the same length. The boundaries won't
    #> be accounted for so subtract 2.
    b = fill(+2.0/(dx*dx),length(grid_x)-2)
    #> The upper and lower diagonals are shorted by 1.
    a = fill(-1.0/(dx*dx),length(b)-1)
    A_Mat = LinearAlgebra.SymTridiagonal(b,a)
    poisson_struct = Poisson_1D(
        A_Mat,
        phi_bound_left,
        phi_bound_right,
        +phi_bound_left/(dx*dx),
        +phi_bound_right/(dx*dx),
    )
    return poisson_struct
end

function Poisson_1D_Solve(
        Poisson_struct,RHS_density
    )
    RHS = copy(RHS_density)
    RHS[2] += Poisson_struct.RHS_boundary_left
    RHS[end-1] += Poisson_struct.RHS_boundary_right
    
    phi = fill(0.0,length(RHS))
    phi[1] = Poisson_struct.phi_boundary_left
    phi[end] = Poisson_struct.phi_boundary_right
    phi[2:end-1] = Poisson_struct.A\RHS[2:end-1];
    return phi
end

end