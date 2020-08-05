using PID
using Plots

t = 0:0.001:1
y = sin.(t)
dy = cos.(t)
iy = -cos.(t) .+ 1 #integral of y from 0 to t
numderiv = zeros(length(t)) 
numint = zeros(length(t))
dt = t[2] - t[1]
for i in eachindex(t)
    if i > 2
        numderiv[i] = PID.updatederivative(y[i], y[i-1], y[i-2], dt)
    end 
end

for i in eachindex(t)
    if i > 1
        numint[i] = PID.updateintegral(numint[i-1], y[i], y[i-1], dt)
    end
end

@assert sum(a -> a^2, dy[3:end] .- numderiv[3:end])/(length(t)-3) < dt^2 
@assert sum(a -> a^2, iy .- numint)/(length(t)-3) < dt^2 


function simulateheating(;setpoint=90.0, maxpower=120.0, kp=0.1, ki=0.0, kd=0.0, steps=10000)
    M = 0.2
    Cp = 800
    kt = 1
    dt = 0.01

    #basic model for the change over time of the temperature of a resistor as a function of its current temperature and power applied
    dT(p, T, dt) = dt*(p - kt*T) / (M*Cp)

    #dT/dt = 0 => p = kt*T => T = p/kt which is the steady state temperature for a given power

    datachannel = Channel{Float64}(1)
    controlchannel = Channel{Float64}(1)
    stopchannel = Channel{Bool}(1)

    # PID.runcontrolloop(kp, ki, kd, dt, setpoint, datachannel, controlchannel, stopchannel)
    @async PID.runcontrolloop(kp, ki, kd, dt, setpoint, datachannel, controlchannel, stopchannel)

    T = 0 #T is the offset from ambient temperature
    temps = Vector{Float64}()
    powers = Vector{Float64}()
    push!(temps, T)
    push!(datachannel, T)
    println("starting steps")
    for i in 1:steps
        if i > 2
            p = take!(controlchannel) * maxpower
        else
            p = 0.0
        end
        T = T + dT(p, T, dt)
        push!(datachannel, T)
        push!(powers, p)
        push!(temps, T)
        # println("Done with step $i of $steps")
    end
    println("Setting trigger to end control loop")
    push!(stopchannel, true)
    println("Adding data to push control loop")
    push!(datachannel, T)
    println("Taking control value to finish control loop")
    take!(controlchannel)
    close(datachannel)
    close(controlchannel)
    close(stopchannel)
    return temps, powers
end

results = [simulateheating(setpoint = 80, maxpower=120, kp = 0.2, ki = ki, steps = 100000) for ki in (0.00005, 0.0001, 0.000115, 0.0002, 0.0004, 0.0008, 0.0016)]
p1 = plot(mapreduce(a -> a[1], hcat, results))
p2 = plot(mapreduce(a -> a[2], hcat, results))
plot(p1, p2, layout = 2)

# steps = 100000
# plot((1:steps).*0.01, simulateheating(p = 100, steps = steps)[2:end])