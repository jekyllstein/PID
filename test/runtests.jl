using PID

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