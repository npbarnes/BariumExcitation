module BariumExcitation

using CSV
using Unitful
using Missings
using Interpolations
using LinearAlgebra
using PhysicalConstants.CODATA2018: m_e, c_0, h

# Transition matrix for the 4 state model (Carlsten 1974)
const k1i = 0.0014u"s^-1"
const k2i = 0.033u"s^-1"
const kti = 0.071u"s^-1"
const k12 = 0.58u"s^-1"
const k21 = 4.05u"s^-1"
const kt1 = 0.4u"s^-1"
const k1t = 0.55u"s^-1"
const k2t = 12.19u"s^-1"
const kt2 = 2.48u"s^-1"
const QQ = [ # Transition Matrix (transposed)
    0u"s^-1"       k1i            k2i            kti;
    0u"s^-1" -(k1i+k12+k1t)       k21            kt1;
    0u"s^-1"       k12      -(k2i+k21+k2t)       kt2;
    0u"s^-1"       k1t            k2t      -(kti+kt1+kt2)
]

# Elementary charge in Franklins
const e = 4.803206799E-10u"g^(1/2)*cm^(3/2)*s^-1"
const C = π*e^2/(m_e*c_0^3*h)
const quadrupole_rate = 2.0u"s^-1"

const solarspectrum = CSV.File("solarspectrum.csv")
const interpolator = linear_interpolation(solarspectrum.Wavelength*u"μm", solarspectrum.var"Energy Flux Density"*u"W*m^-2*μm^-1")
interpolated_flux(λ) = passmissing(interpolator)(λ)

function transition_wavelengths()
    # They're listed in this order because this matches the order of Table 2 in StenBaek-Nielsen (1989).
    λ = Matrix{Union{typeof(1.0u"Å"), Missing}}(undef, 21, 21)
    λ .= missing
    λ[1,4] = 5535.48u"Å"
    λ[6,4] = 15000u"Å"
    λ[15,4] = 11305u"Å"
    λ[1,3] = 3501.11u"Å"
    λ[6,3] = 5826.3u"Å"
    λ[1,2] = 3071.58u"Å"
    λ[6,2] = 4726.44u"Å"
    λ[6,5] = 8559u"Å"
    λ[13,5] = 7120.33u"Å"
    λ[18,5] = 7417.53u"Å"
    λ[6,7] = 6482.91u"Å"
    λ[18,7] = 5805.69u"Å"
    λ[1,9] = 7911.38u"Å"
    λ[13,9] = 27751u"Å"
    λ[15,9] = 29224u"Å"
    λ[1,8] = 3889.33u"Å"
    λ[13,8] = 5997.09u"Å"
    λ[15,8] = 6063.12u"Å"
    λ[6,11] = 6865.69u"Å"
    λ[13,11] = 5907.64u"Å"
    λ[15,11] = 5971.70u"Å"
    λ[18,11] = 6110.78u"Å"
    λ[15,10] = 4591.82u"Å"
    λ[18,10] = 4673.62u"Å"
    λ[1,12] = 4132.43u"Å"
    λ[13,12] = 6595.33u"Å"
    λ[15,12] = 6675.27u"Å"
    λ[13,14] = 6450.85u"Å"
    λ[15,14] = 6527.31u"Å"
    λ[18,14] = 6693.84u"Å"
    λ[15,16] = 6341.68u"Å"
    λ[18,16] = 6498.76u"Å"
    λ[6,19] = 9370u"Å"
    λ[13,19] = 7672.09u"Å"
    λ[15,19] = 7780.48u"Å"
    λ[15,21] = 7280.30u"Å"
    λ[18,21] = 7488.08u"Å"
    λ[13,17] = 3909.91u"Å"
    λ[15,17] = 3937.87u"Å"
    λ[15,20] = 3935.72u"Å"
    λ[18,20] = 3995.66u"Å"
    #λ[1,6] = 8773u"Å" # Quadrupole Transition
    λ
end

function oscillator_strengths()
    # They're listed in this order because this matches the order of Table 2 in StenBaek-Nielsen (1989).
    f = Matrix{Union{Float64, Missing}}(undef, 21, 21)
    f .= missing
    f[1,4] = 1.64
    f[6,4] = .005
    f[15,4] = .02
    f[1,3] = .16
    f[6,3] = .13
    f[1,2] = .18
    f[6,2] = .062
    f[6,5] = .5
    f[13,5] = .15
    f[18,5] = .008
    f[6,7] = .17
    f[18,7] = .011
    f[1,9] = .0102
    f[13,9] = .0086
    f[15,9] = .0147
    f[1,8] = .01
    f[13,8] = .14
    f[15,8] = .18
    f[6,11] = .024
    f[13,11] = .017
    f[15,11] = .084
    f[18,11] = .22
    f[15,10] = .0029
    f[18,10] = .0008
    f[1,12] = .01
    f[13,12] = .26
    f[15,12] = .078
    f[13,14] = .062
    f[15,14] = .21
    f[18,14] = .078
    f[15,16] = .09
    f[18,16] = .31
    f[6,19] = .5
    f[13,19] = .25
    f[15,19] = .067
    f[15,21] = .33
    f[18,21] = .047
    f[13,17] = .106
    f[15,17] = .014
    f[15,20] = .084
    f[18,20] = .012
    f
end

function einstein_coefficients()
    # They're listed in this order because this matches the order of Table 2 in StenBaek-Nielsen (1989).
    # The transpose is because A is for transitions from the higher energy state to the lower
    A_transpose = Matrix{Union{typeof(1.0u"s^-1"), Missing}}(undef, 21, 21)
    A_transpose .= missing
    A_transpose[1,4] = 1.18u"1e8*s^-1"
    A_transpose[6,4] = .002u"1e8*s^-1"
    A_transpose[15,4] = .022u"1e8*s^-1"
    A_transpose[1,3] = .3u"1e8*s^-1"
    A_transpose[6,3] = .43u"1e8*s^-1"
    A_transpose[1,2] = .42u"1e8*s^-1"
    A_transpose[6,2] = .3u"1e8*s^-1"
    A_transpose[6,5] = .46u"1e8*s^-1"
    A_transpose[13,5] = .12u"1e8*s^-1"
    A_transpose[18,5] = .014u"1e8*s^-1"
    A_transpose[6,7] = .19u"1e8*s^-1"
    A_transpose[18,7] = .022u"1e8*s^-1"
    A_transpose[1,9] = .0036u"1e8*s^-1"
    A_transpose[13,9] = .0007u"1e8*s^-1"
    A_transpose[15,9] = .0019u"1e8*s^-1"
    A_transpose[1,8] = .015u"1e8*s^-1"
    A_transpose[13,8] = .26u"1e8*s^-1"
    A_transpose[15,8] = .54u"1e8*s^-1"
    A_transpose[6,11] = .034u"1e8*s^-1"
    A_transpose[13,11] = .02u"1e8*s^-1"
    A_transpose[15,11] = .16u"1e8*s^-1"
    A_transpose[18,11] = .56u"1e8*s^-1"
    A_transpose[15,10] = .009u"1e8*s^-1"
    A_transpose[18,10] = .036u"1e8*s^-1"
    A_transpose[1,12] = .013u"1e8*s^-1"
    A_transpose[13,12] = .40u"1e8*s^-1"
    A_transpose[15,12] = .20u"1e8*s^-1"
    A_transpose[13,14] = .062u"1e8*s^-1"
    A_transpose[15,14] = .33u"1e8*s^-1"
    A_transpose[18,14] = .16u"1e8*s^-1"
    A_transpose[15,16] = .106u"1e8*s^-1"
    A_transpose[18,16] = .48u"1e8*s^-1"
    A_transpose[6,19] =.38u"1e8*s^-1"
    A_transpose[13,19] = .17u"1e8*s^-1"
    A_transpose[15,19] = .073u"1e8*s^-1"
    A_transpose[15,21] = .3u"1e8*s^-1"
    A_transpose[18,21] = .056u"1e8*s^-1"
    A_transpose[13,17] = .27u"1e8*s^-1"
    A_transpose[15,17] = .062u"1e8*s^-1"
    A_transpose[15,20] = .26u"1e8*s^-1"
    A_transpose[18,20] = .049u"1e8*s^-1"

    A = transpose(A_transpose)
    A
end

function branching_ratios()
    R = Matrix{Union{Float64, Missing}}(undef, 21, 21)
    for I in CartesianIndices(R)
        γ, β = Tuple(I)
        R[γ, β] = A[γ,β]/sum(skipmissing(@view(A[γ,:])))
    end
    R
end

function solar_fluxes()
    E = Matrix{Union{typeof(1.0u"W*m^-2*nm^-1"), Missing}}(undef, 22, 22)
    E .= missing
    E[1,4] = 1.26u"W*m^-2*nm^-1"
    E[6,4] = 0.3u"W*m^-2*nm^-1"
    E[15,4] = 0.54u"W*m^-2*nm^-1"
    E[1,3] = 0.92u"W*m^-2*nm^-1"
    E[6,3] = 1.86u"W*m^-2*nm^-1"
    E[1,2] = .56u"W*m^-2*nm^-1"
    E[6,2] = 2.21u"W*m^-2*nm^-1"
    E[6,5] = .99u"W*m^-2*nm^-1"
    E[13,5] = 1.39u"W*m^-2*nm^-1"
    E[18,5] = 1.3u"W*m^-2*nm^-1"
    E[6,7] = 1.55u"W*m^-2*nm^-1"
    E[18,7] = 1.81u"W*m^-2*nm^-1"
    E[1,9] = 1.17u"W*m^-2*nm^-1"
    E[13,9] = 0.04u"W*m^-2*nm^-1"
    E[15,9] = 0.03u"W*m^-2*nm^-1"
    E[1,8] = 0.62u"W*m^-2*nm^-1"
    E[13,8] = 1.78u"W*m^-2*nm^-1"
    E[15,8] = 1.77u"W*m^-2*nm^-1"
    E[6,11] = 1.47u"W*m^-2*nm^-1"
    E[13,11] = 1.85u"W*m^-2*nm^-1"
    E[15,11] = 1.78u"W*m^-2*nm^-1"
    E[18,11] = 1.76u"W*m^-2*nm^-1"
    E[15,10] = 2.19u"W*m^-2*nm^-1"
    E[18,10] = 2.14u"W*m^-2*nm^-1"
    E[1,12] = 1.55u"W*m^-2*nm^-1"
    E[13,12] = 1.60u"W*m^-2*nm^-1"
    E[15,12] = 1.57u"W*m^-2*nm^-1"
    E[13,14] = 1.65u"W*m^-2*nm^-1"
    E[15,14] = 1.47u"W*m^-2*nm^-1"
    E[18,14] = 1.56u"W*m^-2*nm^-1"
    E[15,16] = 1.69u"W*m^-2*nm^-1"
    E[18,16] = 1.60u"W*m^-2*nm^-1"
    E[6,19] = 0.84u"W*m^-2*nm^-1"
    E[13,19] = 1.23u"W*m^-2*nm^-1"
    E[15,19] = 1.20u"W*m^-2*nm^-1"
    E[15,21] = 1.35u"W*m^-2*nm^-1"
    E[18,21] = 1.27u"W*m^-2*nm^-1"
    E[13,17] = 0.53u"W*m^-2*nm^-1"
    E[15,17] = 0.74u"W*m^-2*nm^-1"
    E[15,20] = 0.40u"W*m^-2*nm^-1"
    E[18,20] = 1.89u"W*m^-2*nm^-1"
    E
end

const λ = transition_wavelengths()
const f = oscillator_strengths()
const A = einstein_coefficients()
const R = branching_ratios()
const E = solar_fluxes()

# The k functions are the effective transition rates used by Carlsten (1974) and Stenbaek-Nielsen (1989)
#k(α, γ, β) = C * E(λ[α, γ]) * λ[α,γ]^3*f[α,γ]*R[γ,β]
#k(α, β) = sum(skipmissing( k(α, γ, β) for γ in 1:21 ); init=0.0u"s^-1")

function upward_transition_rate(α, β; fluxes)
    if fluxes=="modern"
        flux = interpolated_flux(λ[α,β])
    elseif fluxes=="classic"
        flux = E[α,β]
    else
        error("`fluxes` must be either \"modern\" or \"classic\"")
    end
    C * flux * λ[α, β]^3 * f[α, β]
end

function downward_transition_rate(α, β)
    if α == 6 && β == 1
        return quadrupole_rate
    end
    A[α, β]
end

function transition_rate(α, β; fluxes)
    up = upward_transition_rate(α, β; fluxes)
    down = downward_transition_rate(α, β)
    if ismissing(up) && !ismissing(down)
        return down
    elseif !ismissing(up) && ismissing(down)
        return up
    elseif ismissing(up) && ismissing(down)
        return 0.0u"s^-1"
    else
        error("A transition can't be both upward and downward")
    end
end

function ionization_rate(α)
    if α == 1
        return 1/700u"s"
    elseif α == 6
        return 1/30u"s"
    elseif α ∈ (13,15,18)
        return 1/14u"s"
    else
        return 0.0u"s^-1"
    end
end

function transition_matrix(; states="22", fluxes="modern")
    if states ∉ ("4","22")
        error("`states` must be either \"4\" or \"22\"")
    end
    if fluxes ∉ ("modern", "classic")
        error("`fluxes` must be either \"modern\" or \"classic\"")
    end
    if states == "4"
        if fluxes=="classic"
            return QQ
        else
            error("4 state model with modern fluxes not implemented")
        end
    end

    T = typeof(1.0u"s^-1")
    Q = Matrix{T}(undef, 22, 22)
    S = Matrix{T}(undef, 21,21)
    S .= Inf*u"s^-1"
    for j in 1:21
        for i in 1:21
            if i == j
                continue
            end
            S[i,j] = transition_rate(j,i; fluxes)
        end
    end
    for i in 1:21
        @views S[i,i] = -(sum(S[1:i-1,i]) + sum(S[i+1:end, i]) + ionization_rate(i))
    end
    Q = [
        zero(T) ionization_rate.(1:21)...
        zeros(T, (21,1)) S
    ]
    Q
end

const Q = transition_matrix()

end # module BariumExcitation
