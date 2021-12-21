
"""
    Constellation

Constellation type for digital communications simulation."""
struct Constellation
    label   :: String
    symbols :: Vector{ComplexF64}  # symbol alphabet
end

avgenergy(v) = mean(abs2.(v))
avgenergy(c::Constellation) = avgenergy(c.symbols)
Base.length(c::Constellation) = length(c.symbols)
Base.getindex(c::Constellation, index) = c.symbols[index]
Base.rand(c::Constellation, n::Integer) = rand(c.symbols, n)
Base.rand(c::Constellation) = rand(c.symbols, 1)
#Base.show(c::Constellation) = plot(c.symbols, w="p")

""" pam(; M = 2, Es = 1.0, ϕ = 0.0) -> Constellation

Return a PAM constellation with `M` elements and average symbol energy `Es`.
The constellation is rotated ϕ radians."""
function pam(; M = 2, Es = 1.0, ϕ = 0.0)
    isinteger(log2(M)) || error("M must be a power of 2 (M = $M)")
    v = M/2
    base = range(-(2v-1), (2v-1), step = 2)
    E = avgenergy(base)
    symbols = exp(1im*ϕ)*sqrt(Es/E).*base
    return Constellation("$M-PAM", collect(symbols))
end

""" squareqam(; M = 4, Es = 1.0, ϕ = 0.0)

Return a square QAM constellation with `M` elements and average symbol energy `Es`.
The constellation is rotated ϕ radians."""
function squareqam(; M = 4, Es = 1.0, ϕ = 0.0)
    iseven(log2(M)) || error("M must be an even power of 2 (M = $M)")
    M2 = Int(log2(M))
    c = real(pam(M=M2).symbols)
    symbols = [ComplexF64(i, q) for i in c for q in c]
    E = avgenergy(symbols)
    symbols .*= sqrt(Es/E)*exp(1im*ϕ)
    return Constellation("$M-QAM", symbols)
end

""" psk(; M = 4, Es = 1.0, ϕ = 0.0)

Return a PSK constellation with `M` elements and average symbol energy `Es`.
The constellation is rotated ϕ radians."""
function psk(; M = 4, Es = 1.0, ϕ = 0.0)
    isinteger(log2(M)) || error("M must be a power of 2 (M = $M)")
    θ = range(0, step = 2π/M, length = M)
    symbols = exp.(1im*θ)
    E = avgenergy(symbols)
    symbols .*= sqrt(Es/E)*exp(1im*ϕ)
    return Constellation("$M-PSK", symbols)
end