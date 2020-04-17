module ParticleGridInteractionModule

import FFTW
import ParticlesModule
import GridModule

function ScatterFieldToParticle_1D_single(
        field_x::Array{Float64,1},
        particle_x::Float64,
        grid_x::Array{Float64,1}
    )
    ind = searchsortedlast(grid_x, particle_x)
    ind = min(ind,length(grid_x)-1)
    ind = max(ind,1)
    weight_upper = (particle_x-grid_x[ind])/(grid_x[ind+1]-grid_x[ind])
    weight_lower = 1.0-weight_upper
    output_quantity_tmp = (weight_lower*field_x[ind])+(weight_upper*field_x[ind+1])
    return output_quantity_tmp
end

function ScatterFieldToParticle_Ex!(
        electric_field_x::Array{Float64,1},
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    for ip in 1:particle_array.number_particles
        particle_array.particle[ip].Ex = ScatterFieldToParticle_1D_single(electric_field_x,particle_array.particle[ip].x,grid_array.x)
    end
end

function ScatterParticleToGrid_1D_single(
        particle_x::Float64,
        grid_x::Array{Float64,1}
    )
    output_quantity_tmp = zero(grid_x)
    ind = searchsortedlast(grid_x, particle_x)
    ind = min(ind,length(grid_x)-1)
    ind = max(ind,1)
    weight_upper = (particle_x-grid_x[ind])/(grid_x[ind+1]-grid_x[ind])
    weight_lower = 1.0-weight_upper
    output_quantity_tmp[ind] = weight_lower
    output_quantity_tmp[ind+1] = weight_upper
    return output_quantity_tmp
end

function ScatterParticleToGrid_2D_single(
        particle_x::Float64,
        grid_x::Array{Float64,1},
        particle_y::Float64,
        grid_y::Array{Float64,1},
    )
    output_quantity_tmp = zeros(Float64,length(grid_x),length(grid_y))
    
    ind_x = searchsortedlast(grid_x, particle_x)
    ind_x = min(ind_x,length(grid_x)-1)
    ind_x = max(ind_x,1)
    weight_upper_x = (particle_x-grid_x[ind_x])/(grid_x[ind_x+1]-grid_x[ind_x])
    weight_lower_x = 1.0-weight_upper_x
    
    ind_y = searchsortedlast(grid_y, particle_y)
    ind_y = min(ind_y,length(grid_y)-1)
    ind_y = max(ind_y,1)
    weight_upper_y = (particle_y-grid_y[ind_y])/(grid_y[ind_y+1]-grid_y[ind_y])
    weight_lower_y = 1.0-weight_upper_y
    
    output_quantity_tmp[ind_x,ind_y] = weight_lower_x*weight_lower_y
    output_quantity_tmp[ind_x+1,ind_y] = weight_upper_x*weight_lower_y
    output_quantity_tmp[ind_x+1,ind_y+1] = weight_upper_x*weight_upper_y
    output_quantity_tmp[ind_x,ind_y+1] = weight_lower_x*weight_upper_y
    return output_quantity_tmp
end

function ScatterParticleToGrid_x(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    output_quantity=zero(grid_array.x)
    for ip in 1:particle_array.number_particles
        output_quantity[:] += ScatterParticleToGrid_1D_single(
            particle_array.particle[ip].x,grid_array.x
        )*particle_array.particle[ip].f_over_g
    end
    return output_quantity
end

function ScatterParticleToGrid_x_vx(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array,
        vgrid_array::GridModule.Grid_Array
    )
    output_quantity = zeros(Float64,length(grid_array.x),length(vgrid_array.x))
    
    for ip in 1:particle_array.number_particles
        output_quantity[:,:] += ScatterParticleToGrid_2D_single(
            particle_array.particle[ip].x,grid_array.x,
            particle_array.particle[ip].vx,vgrid_array.x
        )*particle_array.particle[ip].f_over_g
    end
    output_quantity=transpose(output_quantity)
    return output_quantity
end

function ApplyFiniteParticleShape_x(
        field_x::Array{Float64,1},
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    field_fft=FFTW.rfft(field_x[1:end-1])
    
    k_fft=range(0,stop=length(field_fft)-1,length=length(field_fft))
    k_fft*=2.0*pi/(grid_array.x[end]-grid_array.x[1])
    
    field_fft_out = field_fft.*exp.(-(k_fft*particle_array.particle_shape_x).^2)
    
    field_out=FFTW.irfft(field_fft_out,length(field_x[1:end-1]))
    push!(field_out,field_out[1])

    return field_out
end

function ApplyPeriodicParticleBoundary!(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    for ip in 1:particle_array.number_particles
        particle_array.particle[ip].x = mod(particle_array.particle[ip].x-grid_array.x[1],grid_array.x[end]-grid_array.x[1])+grid_array.x[1]
        particle_array.particle[ip].y = mod(particle_array.particle[ip].y-grid_array.y[1],grid_array.y[end]-grid_array.y[1])+grid_array.y[1]
        particle_array.particle[ip].z = mod(particle_array.particle[ip].z-grid_array.z[1],grid_array.z[end]-grid_array.z[1])+grid_array.z[1]
    end
end

end