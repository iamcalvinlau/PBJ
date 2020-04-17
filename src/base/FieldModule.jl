module FieldModule

import FFTW
import GridModule

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

function ApplyPeriodicBoundary!(
        field_x::Array{Float64,1}
    )
    field_x[1] += field_x[end]
    field_x[end] = field_x[1]
end

end