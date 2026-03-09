using DifferentialEquations
using Plots

# Parameters
m = 1.0        # mass of car
M = 0.2        # mass of pendulum
l = 0.5
g = 9.81
I = (1/3)*M*l^2

A = 2        # amplitude of force
ω = 1.5       # frequency

function car_pendulum!(du,u,p,t)

    x = u[1]
    xdot = u[2]
    θ = u[3]
    θdot = u[4]

    F = A*sin(ω*t)

    A11 = m + M
    A12 = M*l*cos(θ)
    A21 = M*l*cos(θ)
    A22 = M*l^2 + I

    b1 = F + M*l*(θdot^2)*sin(θ)
    b2 = M*g*l*sin(θ)

    A_mat = [A11 A12; A21 A22]
    b_vec = [b1; b2]

    acc = A_mat\b_vec

    xdd = acc[1]
    θdd = acc[2]

    du[1] = xdot
    du[2] = xdd
    du[3] = θdot
    du[4] = θdd
end

tspan = (0.0,20.0)

u0 = [0.0,0.0,0.3,0.0]

prob = ODEProblem(car_pendulum!,u0,tspan)

sol = solve(prob)

anim = @animate for i in 1:5:length(sol.t)

    x = sol[1,i]
    θ = sol[3,i]

    px = x + l*sin(θ)
    py = l*cos(θ)

    plot([-20,20],[0,0], lw=3, color=:black, legend=false)
    scatter!([x],[0], markersize=8)        # car

    plot!([x,px],[0,py], lw=3)              # pendulum rod
    scatter!([px],[py], markersize=6)       # pendulum mass

    xlims!(-1,20)
    ylims!(-0.9,0.9)

end

gif(anim,"pendulum_cart.gif",fps=2)