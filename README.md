# WaveformCommunications

Generate baseband waveforms using the most common pulses and constellations.

Supported pulses: half-sine, cosine, raised cosine, square-root raised cosine,
gaussian.

Supported constellations: any PAM, square QAM and PSK constellation, with arbitrary
energy per symbol and rotation.

## Example

Generate a waveform to transmit 1000 symbols from a 16-QAM modulation, at 100 baud, with
symbol energy 2 J, and using half-sine pulses.

1. Instantiate the constellation

```
using WaveformCommunications

c = qam(M = 16, Es = 2)
```

2. Generate the pulse function.

```
Rp = 10  # transmit at 10 baud
p = halfsinepulse(1/Rp)
```

3. Instantiate the pulse shaper

```
nsyms = 100   # number of symbols to generate
ps = pulseshaper(c, p, nsyms)
```

4. Generate and plot the waveform

The waveform can be sampled at any desired rate, and evaluated at any desired time `t`. Note that the waveform is complex.

```
fs = 500 # sampling frequency
t = range(0, nsyms/Rp, step=1/fs)
waveform = ps.(t)
plot(t, real(waveform))
```

5. Plot the eye diagram

The waveform's eye diagram may be plotted by running

```
ed = eyediag(wf, 1/Rp)
```

and plotting the columns of `ed`.