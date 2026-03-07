using DifferentialEquations
using Plots

# Parameters
m = 1.0      # mass of the cart
M = 0.2      # mass of the pendulum
l = 0.5      # pendulum length
g = 9.81
I = (1/3)*M*l^2   # inertia of the pendulum

function cart_pendulum!(du,u,p,t)

    x = u[1]
    xdot = u[2]
    θ = u[3]
    θdot = u[4]

    # Matrix coefficients
    A11 = m + M
    A12 = M*l*cos(θ)
    A21 = M*l*cos(θ)
    A22 = M*l^2 + I

    # Right side terms
    b1 = M*l*(θdot^2)*sin(θ)
    b2 = M*g*l*sin(θ)

    A = [A11 A12; A21 A22]
    b = [b1; b2]

    acc = A\b

    xdd = acc[1]
    θdd = acc[2]

    du[1] = xdot
    du[2] = xdd
    du[3] = θdot
    du[4] = θdd
end

# Simulation time
tspan = (0.0,10.0)

# Initial conditions
u0_1 = [0.0,0.0,0.2,0.0]
u0_2 = [0.0,0.0,0.5,0.0]
u0_3 = [0.0,0.0,1.0,0.0]

# Define problems
prob1 = ODEProblem(cart_pendulum!,u0_1,tspan)
prob2 = ODEProblem(cart_pendulum!,u0_2,tspan)
prob3 = ODEProblem(cart_pendulum!,u0_3,tspan)

# Solve
sol1 = solve(prob1)
sol2 = solve(prob2)
sol3 = solve(prob3)

# Plot cart position
plot(sol1.t,sol1[1,:],label="IC1",xlabel="t",ylabel="x",title="Cart position")
plot!(sol2.t,sol2[1,:],label="IC2")
plot!(sol3.t,sol3[1,:],label="IC3")

# Plot pendulum angle
plot(sol1.t,sol1[3,:],label="IC1",xlabel="t",ylabel="θ",title="Pendulum angle")
plot!(sol2.t,sol2[3,:],label="IC2")
plot!(sol3.t,sol3[3,:],label="IC3")