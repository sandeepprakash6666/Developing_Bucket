using JuMP
using Ipopt

##Supplementary Files

include("CollMat.jl")

##Required Parameters

Nx = 2;
Nu = 1;
NFE = 20;
NCP = 3;
x0= [0,0];
dx0=[0,0];
u0 = 0;
q0 = 0;
dq0 = 0;
##For Parameters----------------------------------------
#pm = [0.8, 1.0, 1.2];
#p = zeros(Ns, NFE)
#p[1:M:end, :] .= pm[1]








p = [0.8, 1.0, 1.2, 0.8, 1.0, 1.2, 0.8, 1.0, 1.2];


##Probabilites----------------------------------------
pr = [1/3 1/3 1/3];
w = zeros(size(pr, 2)^2);
a = 1;
for i in 1:size(pr, 2)
    global a
    for j in 1:size(pr, 2)
        w[a] = pr[i] * pr[j]
        a = a+1;
    end
end
##--------------------------------------------



MColl = Collocation_Matrix(NCP)

#Multi Stage Data
M  = 3;
Nr = 2;
Ns = M^Nr; 
##Model Defining

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
    set_lower_bound(x[ns, nx, nfe, ncp], 0)
    set_upper_bound(x[ns, nx , nfe, ncp], 10)
    set_lower_bound(u[ns, nu, nfe], -1)
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
    ODE1[ns in 1:Ns, nfe in 1:NFE, ncp in 1:NCP, nu in 1:Nu],    dx[ns, 1, nfe, ncp]      == (p[ns] - x[ns, 2, nfe, ncp] ^2) * x[ns, 1, nfe, ncp] - x[ns, 2, nfe, ncp] + u[ns, nu, nfe];
    ODE2[ns in 1:Ns, nfe in 1:NFE, ncp in 1:NCP],    dx[ns, 2, nfe, ncp]      == x[ns, 1, nfe, ncp];
    #ODEn[nfe in 1:NFE, ncp in 1:NCP],   ...
end)
##Quadrature---------------------------------------------------------
@NLconstraints(m1, begin
      Quad_Diff[ns in 1:Ns, nfe in 1:NFE, ncp in 1:NCP, nu in 1:Nu], dq[ns, 1, nfe, ncp]  == (x[ns, 1, nfe, ncp])^2 + (x[ns, 2, nfe, ncp])^2 + (u[ns, nu, nfe])^2
end)
##Input Penalty---------------------------------------------------------
@NLconstraints(m1, begin

      Input_Penalty_1[ns in 1:Ns, nfe in 1:1, nu in 1:Nu],     -0.05 <=   u[ns, nu, nfe] - u0[1]       <=  0.05
      Input_Penalty_2[ns in 1:Ns, nfe in 2:NFE, nu in 1:Nu],   -0.05 <=   u[ns, nu, nfe] - u[ns, nu, nfe-1]    <=  0.05
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
    for j in 1:M^(i-1)
        for k in 1:M-1
            @constraint(m1, u[k + M * (j-1), 1:Nu, i] .== u[k+1 + M * (j-1), 1:Nu, i])
        end
    end
end
=#

@constraint(m1, u[1 : Ns - 1, 1:Nu, 1] .== u[2 : Ns, 1:Nu, 1]);

if Nr >= 2
    for i in 2:Nr
        for j in 1:M^(Nr-1)
                @constraint(m1,  u[(j-1) * M + 1 : j * M - 1, 1:Nu, i] .== u[(j-1) * M + 2 : j * M, 1:Nu, i])
                #@constraint(m1, constr_non[]  u[(j-1) * M + 1 : j * M - 1, 1:Nu, i] .== u[(j-1) * M + 2 : j * M, 1:Nu, i])
        end
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
##---------------------------------------------------------
##---------------------------------------------------------
##---------------------------------------------------------
##---------------------------------------------------------
##---------------------------------------------------------
##---------------------------------------------------------

#=

temp_pred_constraints = [ @constraint(m, dot(TM[i,:], u - u_cur_v) + T_cur[i] == T_pred[i]) for i = 1:128]


for i = 1:128
    set_name(temp_pred_constraints[i], "temp_pred_constraints[$(i)]")
end
=#