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
println("Derivative test passed")
@assert sum(a -> a^2, iy .- numint)/(length(t)-3) < dt^2 
println("Integration test passed")

function simulateheating(;setpoint=90.0, pb = 10.0, maxpower=120.0, kp=0.1, ki=0.0, kd=0.0, steps=10000)
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
    @async PID.runcontrolloop(kp, ki, kd, pb, dt, setpoint, datachannel, controlchannel, stopchannel)

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

numsteps = 25000
maket(numsteps) = dt .*(0:numsteps) 
zoomind = 17550:17600

#example with infinitely narrow proportional band, corresponds to a system that is either on or off
result1 = simulateheating(setpoint=80, maxpower=120, kp=1.0, pb=0.0, steps = numsteps)
p1 = plot(maket(numsteps), result1[1], legend=false, xaxis = "time", yaxis = "temperature")
p2 = plot(maket(numsteps)[2:end], result1[2], legend=false, xaxis="time", yaxis="power")
p3 = plot(maket(numsteps)[zoomind], result1[1][zoomind], legend=false, xaxis="time", yaxis="temperature")
p4 = plot(maket(numsteps)[zoomind], result1[2][zoomind], legend=false, xaxis="time", yaxis="power")
p = plot(p1, p2, p3, p4, layout=4)
savefig(p, "0degreeband.png")

#example with a proportional band of 10 degrees so power varies smoothly from 
#maximum at 10 degrees too low until 0 at the setpoint
result2 = simulateheating(setpoint=80, maxpower=120, kp=0.1, pb=10.0, steps = numsteps)
p1 = plot(maket(numsteps), result2[1], legend=false, xaxis = "time", yaxis = "temperature")
p2 = plot(maket(numsteps)[2:end], result2[2], legend=false, xaxis="time", yaxis="power")
p = plot(p1, p2, layout=2)
savefig(p, "10degreeband.png")

#10 degree band but with an integral constant as well to correct for a system
#remaining stable below the setpoint
numsteps = 100000
for ki = (0.0, 0.001, 0.002, 0.004, 0.01, 0.1, 1.0)
    result3 = simulateheating(setpoint=80, maxpower=120, kp=0.1, pb=10.0, ki=ki, steps = numsteps)
    p1 = plot(maket(numsteps), result3[1], legend=false, xaxis = "time", yaxis = "temperature")
    p2 = plot(maket(numsteps)[2:end], result3[2], legend=false, xaxis="time", yaxis="power")
    p = plot(p1, p2, layout=2)
    savefig(p, "10degreeband_ki=$ki.png")
end