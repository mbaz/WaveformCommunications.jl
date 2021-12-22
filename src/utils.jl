"""
eyediag(wf, Tp, sps::Integer=10 ; traces = 100, n = 3, tini = 0.0)

Return the eye diagram of the waveform `wf`.

This function returns a matrix with `n*sps` rows and `traces` columns. Each column
corresponds to one trace of the waveform `wf`; each trace includes `n` symbols,
each of duration `Tp`. Each symbol is sampled `sps` times. The eye diagram is
evaluated starting at time `tini`.

To obtain a plot of the eye diagram, plot the columns of the returned matrix.
"""
function eyediag(wf, Tp, sps::Integer=10 ; traces = 100, n = 3, tini = 0.0)
    ed = zeros(ComplexF64, sps*n, traces)
    tend = tini + n*Tp
    for tr in 1:traces
        t = range(tini, tend, length=sps*n)
        ed[:,tr] = wf.(t)
        tini = tend
        tend = tini + n*Tp
    end
    return ed
end