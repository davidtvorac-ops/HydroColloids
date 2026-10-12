using Statistics, Random, LinearAlgebra, StaticArrays, Roots
using DelimitedFiles, Dates
using Plots, LaTeXStrings

include("init_config.jl")                              # Initializes the system
include("force.jl")                             # Calculates the external force at position x
include("equilibrate.jl")                           
include("calculateCollisionTime.jl")
include("calculateCollisions.jl")
# include("step.jl")
include("run.jl")
include("angular_displacement.jl")

### PARAMETERS ###
struct Collision
    time::Float64
    i::Integer
    j::Integer
end
Base.isless(c1::Collision, c2::Collision) = c1.time < c2.time

struct Parameters
    N::Integer                      # number of particles
    L::Float64                      # Circumference of the potential
    Tr::Integer                     # Number of optical traps
    U0::Float64                     # amplitude of the cosine potential
    w::Float64                      # half-width of the gaussian potential
    μ₀::Float64                     # bare diffusion coefficient
    dt::Float64                     # time-step of the simulation
    teq::Float64                    # equilibration duration at the beginning of the simulation (recommended >200)
    tf::Float64                     # duration of the simulation after equilibration
    τ::Float64                      # period of rotation switching
    R::Float64                      # Radius of optic traps system
    a::Float64                      # radius of the particles
    a_eff::Float64                  # percentage of radius of the particles For HI
    λ::Float64                      # Wavelength in μm
    z::Float64                      # distance of the particles from the coverslip surface (z = 1.1λ)
    ω::Float64                      # angular velocity of rotation
    wall::Bool                   # whether wall effects are considered
    η::Float64                      # viscosity
    μbase::Matrix{Float64}          # self-part of the mobility tensor
    RW::Matrix{Float64}             # Auxiliary term for exact wall-single particle 
end

### Walls effects not-included #### 
include("master_run.jl");
wall = true
runs = 1
include("tensors.jl")

include("velocities.jl"); include("step.jl");
τval= [0.23, 0.255, 0.28, 0.305, 0.33, 0.355, 0.38, 0.405, 0.43, 0.455, 0.5, 0.55, 0.6, 0.65, 0.7, 0.85, 1]

for Δω in [-6,6]
Jcache = master(Δω,wall,τval,runs)
writedlm("Osnabruck/test,rot_$(sign(Δω)).csv",[τval Jcache],';')
end


J₊ = readdlm("Osnabruck/test,rot_-1.csv",';')
J₋ = readdlm("Osnabruck/test,rot_1.csv",';')

# N = 6
data = readdlm("Jonas/summed-flux_13_12_7_6_1.csv",',')

diam = Shape([-0.5,0.5,0.5,-0.5,-0.5],[-0.5,-0.5,0.5,0.5,-0.5])
plot(J₊[:,1],mean(J₊[:,2:end],dims=2),label="ω = 6",lw=3,markershape=:circle,c=:blue4)
scatter!(data[2:end,1],data[2:end,6],markersize=14,msw=2,markershape=diam,color=:transparent,markerstrokecolor=:black)

plot!(J₋[:,1],mean(J₋[:,2:end],dims=2),label="ω = -6",lw=3,markershape=:circle,c=:red4)
scatter!(data[2:end,1],data[2:end,5],markersize=14,msw=2,markershape=diam,color=:transparent,markerstrokecolor=:black)






#### ANIMATION #####
####################
using Plots, DelimitedFiles, Statistics, LaTeXStrings, Measures
j = [-11:0.25:11;]
#create 2d grid from j
grd = [[y,x,p.z] for x in j, y in j]
fgrd = pot(p, grd, 0)
# force!(f, p, grd[:], t, traps)
# h = contour(j,j,reshape(norm.(reinterpret(SVector{3,Float64},f)),201,201),levels=10, lw=1.5,
#         aspect_ratio=1, size=(404,400),widen=false,
#         left_margin=-10pt,bottom_margin=-10pt,
#         axis=false,legend=false,grid=false)


        
# # h = contourf(j,j,norm.(fgrd), linealpha=0.0, linewidth=0,
# #         aspect_ratio=1, size=(404,400),widen=false,
# #         left_margin=-10pt,bottom_margin=-10pt,c=cgrad(:thermal,rev=false),
# #         axis=false,legend=false,grid=false)
# # scatter!(h,x[1:3:end],x[2:3:end],markersize=20,msw=1.5,alpha=0.65,c=:green4)
# surface(j,j,-norm.(fgrd), size=(404,400),widen=false,alpha=0.75,
#         left_margin=-10pt,bottom_margin=-10pt,c=cgrad([:khaki2,:khaki3,:darkorange2,:darkorange2,:darkorange3,:darkorange3,:darkorange3,:darkorange3,:darkorange3],rev=true),
#         zlims=(-800,0),
#         clims=(-200,0),
#         camera=(10,35),
#         lw=10,
#         grid=false,
#         axis=false,legend=false)

h = contour(j,j,norm.(fgrd),levels=[50;0:100:600;], lw=3,
        aspect_ratio=1, size=(404,400),widen=false,
        left_margin=-10pt,bottom_margin=-10pt,c=:rebeccapurple,#cgrad(:thermal,rev=false),
        axis=false,legend=false,grid=false,zlims=(0,600))

min_k = 18001
len_k = size(x_traj,3)
max_k = 20001
anim = @animate for k in min_k:20:max_k
    t = k/len_k*p.tf#*(max_k-min_k)/len_k
    hh = deepcopy(h)
    ω₀ = Δω*2π/τ
    x_rot = [cos(ω₀*t) sin(ω₀*t); -sin(ω₀*t) cos(ω₀*t)] * x_traj[1:2,:,k]
    scatter!(hh,x_rot[1,:],x_rot[2,:],markersize=18,msw=3.5,alpha=0.55,c=:green4,fontfamily="Computer Modern",title="Co-rotating frame")
    annotate!(hh, -3, 0.75, text("$(sign(Δω) < 0 ? "↻"*L"\sim J_+" : "↺"*L"\sim J_-" )",:left, 14, :black))
    annotate!(hh, -3, -0.75, text(L"t = %$(round(t-min_k/len_k*p.tf,digits=2))\,\mathrm{s}",:left, 14, :black))
    annotate!(hh, -3, -2.25, text(L"ω = %$(round(abs(ω₀),digits=2))\,\mathrm{rad\,s^{-1}}",:left, 14, :black))

    fgrd = pot(p, grd, t)
    lab_hh = contour(j,j,norm.(fgrd),levels=[50;0:100:600;], lw=3,
        aspect_ratio=1, size=(404,400),widen=false,
        left_margin=-10pt,bottom_margin=-10pt,c=:rebeccapurple,#cgrad(:thermal,rev=false),
        axis=false,legend=false,grid=false,zlims=(0,600))
    scatter!(lab_hh,x_traj[1,:,k],x_traj[2,:,k],markersize=18,msw=3.5,alpha=0.55,c=:green4,familyfont="Computer Modern",title="Laboratory frame")
    plot(hh,lab_hh,layout=(1,2),size=(808,400))
end
# gif(anim, "anim.gif", fps=15)
# # gif(anim, "anim/both_noisy/omega=$(τval[k]).gif", fps=15)
mp4(anim, "anim.mp4", fps=15)








####### AZIMUTHAL FORCE PLOT ######
###################################
# pp = plot()
# hline!(pp,[0],c=:black,lw=1,label=false)
# vline!(pp,[0],c=:black,lw=1,label=false)
#     scatter!(pp,[-100],[100], markersize=3, msw=0,  alpha=0.55, c=:firebrick1,label=" No HI" )
#     scatter!(pp,[-100],[100], markersize=3, msw=0,  alpha=0.55, c=:blue4,label=" Far HI")
# mu0 = 0.061/2.1^2
# anim = @animate for d = -pi/12:0.005:pi/12
#     pol_rot = copy(pol)     
#     hh = deepcopy(h)
#     pol_rot[1:2,2] = [cos(d) -sin(d); sin(d) cos(d)] * pol[1:2,2]
#         xs = reinterpret(SVector{3,Float64}, pol_rot)
#     tensors(p, xs, μ)
#     force!(f, p, xs, 0, traps)
#     r = pol_rot[:,2]
#     eφ = [-r[2], r[1], 0.0]      # tangential direction (z × r)
#     eφ /= norm(r)               # normalize
#     f0 = 3.15τ
#     f[4:6] .-= f0 * eφ

#     scatter!(hh,pol_rot[1,:],pol_rot[2,:],markersize=20,msw=1.5,alpha=0.55,c=:green4)

#     quiver!(hh,[pol_rot[1,2]], [pol_rot[2,2]], quiver=([(-μ*f)[4]],[(-μ*f)[5]]), alpha=0.55, c=:blue4, lw=2, label=false)

#     pol_perp = cross( [0,0,1],pol_rot[:,2])
#     pol_perp ./= norm(pol_perp)
#     ff = dot((-μ*f)[4:6],pol_perp)
#     scatter!(pp,[d], [ff], xlims=(-pi/12,pi/12),ylims=(-2,6), fontfamily="Computer Modern", markersize=3, msw=0, alpha=0.55, c=:blue4,label=false,grid=false)
    
#     muu = diagm(mu0*ones(3N))
#     ff = dot((-muu*f)[4:6],pol_perp)
#     scatter!(pp,[d], [ff], xlims=(-pi/12,pi/12), markersize=3, msw=0, framestyle=:box, alpha=0.55, c=:firebrick1,label=false)
    
#        quiver!(hh,[pol_rot[1,2]], [pol_rot[2,2]], quiver=([(-muu*f)[4]],[(-muu*f)[5]]),
#         left_margin=10pt, c=:firebrick1, lw=2, label=false)
#     plot!(hh,pp,layout=(1,2),size=(800,400),widen=false,bottom_margin=10pt,ylabel=["" L"(\mathbf{\mu f})_2 \cdot \mathbf{e}_{\varphi_2}"], xlabel=["" L"\varphi - \varphi_0~(\mathrm{rad})"])
# end
# gif(anim, "anim.gif", fps=15)


# mu0 = 0.061/2.1^2
# dd = -pi/12:0.005:pi/12
# nohi = zeros(length(dd))
# farhi = zeros(length(dd))
# hi = zeros(length(dd))
# for (i, d) in enumerate(dd)
#     pol_rot = copy(pol)     
#     pol_rot[1:2,2] = [cos(d) -sin(d); sin(d) cos(d)] * pol[1:2,2]
#     xs = reinterpret(SVector{3,Float64}, pol_rot)

#     force!(f, p, xs, 0, traps)
#     r = pol_rot[:,2]
#     eφ = [-r[2], r[1], 0.0]      # tangential direction (z × r)
#     eφ /= norm(r)               # normalize
#     f0 = 3.15τ
#     f[4:6] .-= f0 * eφ

#     pol_perp = cross( [0,0,1],pol_rot[:,2])
#     pol_perp ./= norm(pol_perp)
#     ff = dot((-μ*f)[4:6],pol_perp)
#     include("tensors.jl")
#     tensors(p, xs, μ)
#     hi[i] = ff
#     include("tensors_far.jl")
#     tensors(p, xs, μ)
#     farhi[i] = ff
#     muu = diagm(mu0*ones(3N))
#     ff = dot((-muu*f)[4:6],pol_perp)
#     nohi[i] = ff
# end

# vline([0], c=:black, lw=1, label=false,ls=:dash)
# hline!([0], c=:black, lw=1, label=false,ls=:dash)
# plot!(dd,hi, xlims=(-pi/12,pi/12), ylims=(-0.35,0.35), c=:dodgerblue2,label=" 2D: far HI", grid=false)

# plot!(dd, nohi,msw=0, framestyle=:box, c=:black,label=" 2D: no HI")
# plot!(xlabel=L"\varphi - \varphi_0~(\mathrm{stable})", ylabel=L"\mathrm{azimuthal\,force}", fontfamily="Computer Modern", size=(500,400), legend=:bottomright)
