module WaveformCommunications

using QuadGK, Statistics

export cosinepulse, halfsinepulse, rcpulse, srrcpulse, gaussianpulse,
       Constellation, pam, squareqam, psk,
       Pulse, pulseshaper

include("pulses.jl")
include("constellations.jl")

"""
    pulseshaper(c, pulse, nsyms = 20)

Return a function of time that generates a pulse-shaped waveform.

Symbols are drawn at random from constellation `c`. The pulse shape is specified by
`pulse`. The waveform starts at time 0 and ends after `nsyms` symbols have been generated.
"""
function pulseshaper(c::Constellation, pulse::Pulse, nsyms = 20)
    syms = rand(c, nsyms)
    return t -> begin
        s = zero(ComplexF64)
        for k = 1:nsyms
            s += syms[k].*pulse.(t.-(k-1)*pulse.Tp)
        end
        return s
    end
end

end