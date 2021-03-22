using JuMP, Ipopt, Plots
##Supplementary Files
include("CollMat.jl")
##Required Parameters
Nx = 2;
Nu = 1;
NFE = 20;
NCP = 3;
x0= [0, 1.0];
dx0=[0,0];
u0 = 0;
q0 = 0;
dq0 = 0;

##For Parameters----------------------------------------

##Multi Stage Data---------------------------------------
M  = 3;
Nr = 3;
Ns = M^Nr; 

p = [0.8, 1.0, 1.2]
w_p = [1/3, 1/3, 1/3]
#Building the parameter matrix
p_mat = NaN*ones(Ns, NFE)
w_mat = NaN*ones(Ns,Nr)

#Till Nr
for nfe in 1:Nr
    # nfe = 3
    global Gs, Bs
    Gs = M^(nfe);
    Bs = M^(Nr - nfe);
    
    for gs in 1:Gs
    # gs = 2
        if mod(gs, M) == 0
            p_mat[(gs-1)*Bs + 1 : (gs)*Bs, nfe]     .=  p[M]
            w_mat[(gs-1)*Bs + 1 : (gs)*Bs, nfe]     .=  w_p[M]
        else
            p_mat[(gs-1)*Bs + 1 : (gs)*Bs, nfe]     .=  p[mod(gs,M)]
            w_mat[(gs-1)*Bs + 1 : (gs)*Bs, nfe]     .=  w_p[mod(gs,M)]
            # @constraint(m1,  u[(gs-1)*Bs+1 : (gs)*Bs-1, 1:Nu, nfe] .== u[(gs-1)*Bs+2 : (gs)*Bs, 1:Nu, nfe])
        end
end
end

p_mat
w_mat

#Remaining (from Nr+1 to NFE)
p_mat[:, Nr + 1 : end] .= p_mat[:, Nr]
 
p_mat


#Probability of each scenario
w = NaN*ones(Ns)

for i in 1:Ns
    w[i] = prod( w_mat[i, 1:end] )
end

w

sum(w)





# #p = [0.8, 1.0, 1.2, 0.8, 1.0, 1.2, 0.8, 1.0, 1.2];
# p = [0.8, 1.0, 1.2, 0.8, 1.0, 1.2, 0.8, 1.0, 1.2, 0.8, 1.0, 1.2, 0.8, 1.0, 1.2, 0.8, 1.0, 1.2, 0.8, 1.0, 1.2, 0.8, 1.0, 1.2, 0.8, 1.0, 1.2];


#p = hcat(p, p, p)
MColl = Collocation_Matrix(NCP)

##Probabilites----------------------------------------
# pr = [1/3 1/3 1/3];
# w = NaN*ones(size(pr, 2)^Nr);
# a = 1;
# for i in 1:size(pr, 2)
    # global a
    # for j in 1:size(pr, 2)
        # w[a] = pr[i] * pr[j]
        # a = a+1;
    # end
# end
##Model Defining
# w

m1 = Model(Ipopt.Optimizer)

##------------------------------------------------------------
@variable(m1,  x[1:Ns, 1:Nx, 1:NFE, 1:NCP+1])
@variable(m1, dx[1:Ns, 1:Nx, 1:NFE, 1:NCP])
@variable(m1,  u[1:Ns, 1:Nu, 1:NFE       ])
@variable(m1,  q[1:Ns, 1,    1:NFE, 1:NCP])
@variable(m1, dq[1:Ns, 1,    1:NFE, 1:NCP])
##-----------------------------------------------------------------
for ns in 1:Ns, nx in 1:Nx, nfe in 1:NFE, ncp in 1:NCP, nu in 1:Nu
    #Bounds
    set_lower_bound(x[ns, nx, nfe, ncp], -10)
    set_upper_bound(x[ns, nx , nfe, ncp], 10)
    set_lower_bound(u[ns, nu, nfe], -2)
    set_upper_bound(u[ns, nu, nfe], 2)
    #Start Values
    set_start_value(x[ns, nx, nfe, ncp],    x0[nx])
    set_start_value(dx[ns, nx, nfe, ncp], dx0[nx])
    set_start_value(u[ns, nu, nfe],         u0[nu])
    set_start_value(q[ns, 1, nfe, ncp],     q0)
    set_start_value(dq[ns, 1, nfe, ncp],    dq0)
end
##Set the ODEs-----------------------------------------------
@NLconstraints(m1, begin
    ODE1[ns in 1:Ns, nfe in 1:NFE, ncp in 1:NCP, nu in 1:Nu],    dx[ns, 1, nfe, ncp]      == (p_mat[ns,nfe] - x[ns, 2, nfe, ncp] ^2) * x[ns, 1, nfe, ncp] - x[ns, 2, nfe, ncp] + u[ns, nu, nfe];
    ODE2[ns in 1:Ns, nfe in 1:NFE, ncp in 1:NCP],    dx[ns, 2, nfe, ncp]      == x[ns, 1, nfe, ncp];
    #ODEn[nfe in 1:NFE, ncp in 1:NCP],   ...
    
end)
##Quadrature---------------------------------------------------------
@NLconstraints(m1, begin
      Quad_Diff[ns in 1:Ns, nfe in 1:NFE, ncp in 1:NCP, nu in 1:Nu], dq[ns, 1, nfe, ncp]  == (x[ns, 1, nfe, ncp])^2 + (x[ns, 2, nfe, ncp])^2 + (u[ns, nu, nfe])^2
end)
##Collcation!---------------------------------------------------------
@NLconstraints(m1, begin
        Coll_Eq_Diff[    ns in 1:Ns,     nx in 1:Nx,  nfe in 1:NFE,    ncp in 1:NCP],     x[ns, nx, nfe, ncp+1] == x[ns, nx, nfe, 1] + sum(MColl[ncp, i] * dx[ns, nx, nfe, i] for i in 1:NCP)
        Cont_Eq_First[   ns in 1:Ns,     nx in 1:Nx],                                     x[ns, nx, 1, 1]       == x0[nx]
        Cont_Eq_rest[    ns in 1:Ns,     nx in 1:Nx,  nfe in 2:NFE],                      x[ns, nx, nfe, 1]     == x[ns, nx, nfe-1, end]
        Coll_Eq_Quad0[   ns in 1:Ns,                                   ncp in 1:NCP],     q[ns, 1, 1, ncp]      == q0  + sum(MColl[ncp, i] * dq[ns, 1, 1, i] for i in 1:NCP)
        Coll_Eq_Quad[    ns in 1:Ns,                  nfe in 2:NFE,    ncp in 1:NCP],     q[ns, 1, nfe, ncp]    == q[ns, 1, nfe-1, NCP]  + sum(MColl[ncp, i] * dq[ns, 1, nfe, i] for i in 1:NCP)
end)


##Non-Anticipativity constraints---------------------------------------------------------
#=
for i in 1:Nr
    if i != Nr
        for j in 1 : (M^i)
            @constraint(m1,  u[(j-1)*(M^(Nr-i))+1 : j*(M^(Nr-i))-1, 1:Nu, i] .== u[(j-1)*(M^(Nr-i))+2 : j*(M^(Nr-i)), 1:Nu, i])
        end
    else
        for j in 1: (M^(Nr-1))
            @constraint(m1,  u[(j-1)*M+1 : j*M-1, 1:Nu, i] .== u[(j-1)*M+2 : j*M, 1:Nu, i])
        end
    end
end
=#

for nfe in 1:Nr
    global Gs, Bs
    Gs = M^(nfe - 1);
    Bs = M^(Nr - nfe + 1);
    for gs in 1:Gs
        @constraint(m1,  u[(gs-1)*Bs+1 : (gs)*Bs-1, 1:Nu, nfe] .== u[(gs-1)*Bs+2 : (gs)*Bs, 1:Nu, nfe])
    end
end


##Objective Function---------------------------------------------------------
@NLobjective(m1, Min,  sum(w[ns] * q[ns, 1, end, end] for ns in 1:Ns))
##---------------------------------------------------------
optimize!(m1)
JuMP.termination_status(m1)
JuMP.solve_time(m1)
##---------------------------------------------------------
Solution = JuMP.value.(u)[:,:,:]
x       = JuMP.value.(x)[:,:,:,:]
u_plot  = cat(x0[1], Solution[1:NFE], dims = 1)
x1      = cat(x0[1], x[1,:,:,:], dims = 1)
x2      = cat(x0[2], x[2,:,:,:], dims = 1)
t_plot = collect(1:NFE) 
##---------------------------------------------------------
#choose backend for plots
plotlyjs()




px1 = plot(t_plot, x[1, 1, :, end],         label = "x1_1")
for i in 2:Ns

    px1 = plot!(t_plot, x[i, 1, :, end],         label = "x1_$i")

end
px2 = plot( t_plot, x[1, 2, :, end],         label = "x2_1")
for i in 2:Ns

    px2 = plot!(t_plot, x[i, 2, :, end],         label = "X2_$i")

end

pu = plot(t_plot, Solution[1,1,:], label = "U1", linetype = :steppost)
for i in 2:Ns

    pu = plot!(t_plot, Solution[i,1,:], label = "U_2", linetype = :steppost)

end
fig1 = plot(px1, px2, pu, layout = (3, 1), xaxis = "Time")
##---------------------------------------------------------
#choose backend for plots
plotlyjs()
px = plot(t_plot, x[1, 1, :, end],         label = "x1_1")
for i in 2:Ns

    px = plot!(t_plot, x[i, 1, :, end],         label = "x1_$i")

end
px = plot!( t_plot, x[1, 2, :, end],         label = "x2_1")
for i in 2:Ns

    px = plot!(t_plot, x[i, 2, :, end],         label = "X2_$i")

end

pu = plot(t_plot, Solution[1,1,:], label = "U1", linetype = :steppost)
for i in 2:Ns

    pu = plot!(t_plot, Solution[i,1,:], label = "U_2", linetype = :steppost)

end
fig1 = plot(px, pu, layout = (2, 1), xaxis = "Time")
##---------------------------------------------------------
scatter(Solution[1:end,1,1])
##---------------------------------------------------------
scatter(Solution[1:end,1,2])
##---------------------------------------------------------
scatter(Solution[1:end,1,3])
##---------------------------------------------------------
scatter(Solution[1:end,1,1:3])
##---------------------------------------------------------
scatter(Solution[19:27,1,1])
##---------------------------------------------------------
scatter(Solution[1:end,1,1])
##---------------------------------------------------------
scatter(Solution[4:6,1,2])
##---------------------------------------------------------
scatter(Solution[10:18,1,2])
##---------------------------------------------------------
scatter(Solution[19:27,1,2])
##---------------------------------------------------------
scatter(Solution[4:6,1,2])
##---------------------------------------------------------
scatter(Solution[1:end,1,1:4])