module PID

"""
Updates an existing integral of errors for the previous timestep, i1, with the current time step.  e0 and e1 are the current and previous 
values respectively for the error and dt is the time between each step.  The integral is updated by adding
a term a to the original integral that is the trapezoid rule value given by the two recent values.

All input arguments should be a subtype of AbstractFloat

`inew = updateintegral(i1 e0, e1, dt)` 
"""
function updateintegral(i1::T, e0::T, e1::T, dt::T) where T <: AbstractFloat 
    a = dt * (e0 + e1) / 2
    i1 + a
end

"""
Calculates a derivative value for the errors using the last 3 error values, e0, e1, and e2 where e0 is the current value.
The 2nd order accuracy version of the backward finite difference approximation is used where dt is the time between
each error measurement. See https://www.wikiwand.com/en/Finite_difference_coefficient

All input arguments should be a subtype of AbstractFloat

`d0 = updatederivative(e0::T, e1::T, e2::T, dt::T) where T <: AbstractFloat`
""" 
updatederivative(e0::T, e1::T, e2::T, dt::T) where T <: AbstractFloat = (T(1.5)*e0 - 2*e1 +  T(0.5)*e2) / dt
    

"""
Calculates pid controller output for the current error value e0 based on the user defined coefficients kp, ki, kd.
Requires the most recent 3 error values as well as the previous integral value and the time spacing between measurements.
The controller with clamp the output value between the min and max values for the controller which is 0 to 1 by default.  
In most cases this will correspond to the minimum and maximum power available to apply to the control system.

All input arguments should be a subtype of AbstractFloat

`output = pidcontrol(e0, e1, e2, i1, kp, ki, kd, dt)`
"""
pidcontrol(e::T, i::T, d::T, kp::T, ki::T, kd::T; maxoutput=1.0, minoutput=0.0) where T <: AbstractFloat = clamp(kp*e + ki*i + kd*d, minoutput, maxoutput)

"""
This function is a mockup of a control loop where datachannel represents some place where the controller can get the most recent measured value.
controlchannel represents a place where the controller outputs its desired value.  kp, ki, kd are the user defined constants that regulate the behavior
of the control system.  dt is the fixed time between measurements that appear in datachannel.  setpoint is the desired value for the system to reach.
stopchannel represents a way to halt the control loop if something is placed in the channel.

E.g. 
setpoint = 100 degrees C
datachannel contains temperature measurements at a fixed interval
controlchannel receives the desired power to apply to the heating element.  In practice this could be a value can be normalized to 0 to 1 where 0 is 
no applied power and 1 is the maximum applied power.

The purpose of the control loop is to reach and maintain a given setpoint.  If a new setpoint is desired, the loop can be stopped and restarted with a new setpoint.
This could also be designed so the setpoint is changeable during the loop, I just chose to make it a fixed value. 

Also note that this control loop only measures one input and controls one output so it isn't designed to use a combination of control elements to maintain a setpoint such as 
power to a heating element and a fan.  If the setpoint is lowered one approach would be to stop the pid loop, run the fans at full power until the temperature hits the new setpoint,
and then start a new pid control loop at the new setpoint.
"""
function runcontrolloop(kp, ki, kd, dt, setpoint, datachannel, controlchannel, stopchannel)
    tsteps = 0
    #calculates the difference between the setpoint and the process variable (pv)
    calcerror(pv) = setpoint - pv
    while empty(stopchannel)
        #If during the first 3 iterations fill in e0, e1, e2 and start accumulating integral.  After that point shift the values backwards replacing e0 with the current measurement
        if tsteps == 0
            e2 = calcerror(take!(datachannel))
        elseif tsteps == 1
            e1 = calcerror(take!(datachannel))
            i = updateintegral(0.0, e1, e2, dt)
        elseif tsteps >= 2
        #after 3 steps we have all the information to start outputing control values but past the 3rd step we have to shift back the error values to make room to add a new e0
            if tsteps > 2
                e2 = e1
                e1 = e0
            end
            e0 = calcerror(take!(datachannel))
            i = updateintegral(i, e0, e1, dt) 
            d = updatederivative(e0, e1, d2, dt)
            output = pidcontrol(e0, i, d, kp, ki, kd)
            push!(controlchannel, output)
        end
    end
end

end # module
