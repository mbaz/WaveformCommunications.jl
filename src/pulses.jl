struct Pulse{T,F}
    Tp :: T
    f  :: F
end

(p::Pulse)(t) = p.f(t)

"""cosinepulse(T ; normalize = true)

A cosine pulse of duration `T`. This pulse is also known as a "time-domain
raised cosine" pulse."""

function cosinepulse(T ; normalize = true)
    p1(t) = begin
        if 0 ≤ t < T
            return 1.0-cospi(2*t/T)
        else
            return 0.0
        end
    end
    N = 1.0
    if normalize
        area, _ = quadgk(t -> p1(t)^2, 0.0, T)
        N = 1.0/sqrt(area)
    end
    return Pulse(T, t -> N*p1(t))
end

"""halfinepulse(T ; normalize = true)

A half-sine pulse of duration `T`."""
function  halfsinepulse(T ; normalize = true)
    p1(t) = begin
        if 0 ≤ t < T
            return sinpi(t/T)
        else
            return 0.0
        end
    end
    N = 1.0
    if normalize
        area, _ = quadgk(t -> p1(t)^2, 0.0, T)
        N = 1.0/sqrt(area)
    end
    return Pulse(T, t -> N*p1(t))
end

"""
    rcpulse(T, β, duration = 10T ; normalize = true)

A causal raised-cosine pulse with main lobe duration `T`, roll-off
factor (excess bandwidth) `β`, and duration `d`.
"""
function rcpulse(T, β, duration = 10T ; normalize = true)
    if β > 1 || β < 0
        @error "Excess bandwidth β must be between 0 and 1, inclusive." β
    end
    p1(t) = begin
        t ≥ duration/2 && return 0.0
        t < -duration/2 && return 0.0
        β == 0 && return sinc(t)
        if t ≈ T/(2β) || t ≈ -T/(2β)
            return (π/(4T))*sinc(1/(2β))
        end
        return (1/T)*sinc(t/T)*cospi(β*t/T)/(1-(2β*t/T)^2)
    end
    N = 1.0
    if normalize
        area, _ = quadgk(t -> p1(t)^2, -duration/2, duration/2)
        N = 1.0 / sqrt(area)
    end
    return Pulse(T, t -> N * p1(t - duration/2))
end

"""
    srrcpulse(T, β, duration = 10T ; normalize = true)

A causal square-root raised cosine pulse with main lobe duration `T`, roll-off
factor (excess bandwidth) `β``, and duration `d`.
"""
function srrcpulse(T, β, duration = 10T ; normalize = true)
    if β > 1 || β < 0
        @error "Excess bandwidth β must be between 0 and 1, inclusive." β
    end
    p1(t) = begin
        t ≥ duration/2 && return 0.0
        t < -duration/2 && return 0.0
        β == 0 && return sinc(t)
        t ≈ 0 && return (1 - β + 4β / π) / sqrt(T)
        if t ≈ T / (4β) || t ≈ -T / (4β)
            x = (β * (π + 2)) / (π * sqrt(2T)) * sinpi(1 / (4β))
            y = (β * (π - 2)) / (π * sqrt(2T)) * cospi(1 / (4β))
            return x + y
        else
            x = sinpi((1 - β) * t / T)
            y = (4β * t / T) * cospi((1 + β) * t / T)
            z = sqrt(T) * (π * t / T) * (1 - (4β * t / T)^2)
            return (x + y) / z
        end
    end
    N = 1.0
    if normalize
        area, _ = quadgk(t -> p1(t)^2, -duration/2, duration/2)
        N = 1.0 / sqrt(area)
    end
    return Pulse(T, t -> N * p1(t - duration/2))
end

"""gaussianpulse(B, T, duration=10 ; normalize = true)

A causal gaussian pulse with 3-dB bandwidth `B`, symbol time
`T`, and the specified `duration`.

See Molisch, "Wireless Communications", 11.2.2 and 11.3.9.
"""
function gaussianpulse(B, T, duration = 4T; normalize = true)
    p1(t) = begin
        t ≥ duration/2 && return 0.0
        t < -duration/2 && return 0.0
        σ = sqrt(log(2))/(2π*B*T)
        return 1/(sqrt(2π)*σ*T)*exp(-t^2/(2σ^2*T^2))
    end
    N = 1.0
    if normalize
        area, _ = quadgk(t -> p1(t)^2, -duration / 2, duration / 2)
        N = 1.0 / sqrt(area)
    end
    return Pulse(T, t -> N * p1(t - duration / 2))
end