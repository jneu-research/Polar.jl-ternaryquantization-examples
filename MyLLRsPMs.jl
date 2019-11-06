__precompile__()

module MyLLRsPMs

    # LLRs

    using StaticArrays

    import Base: zero, -
    using Polar.PolarDecoding: LLRType
    using Polar.GF2n: GF2Element, GF2_0, GF2_1
    import Polar.PolarDecoding: cn_update, vn_update
    
    export GenuineLLR, MinApproxLLR, Q3lLLR
    export val, label
    

    
    struct GenuineLLR <: LLRType
        v::Float64
    end

    GenuineLLR() = GenuineLLR(0.00)
    zero(::Type{GenuineLLR}) = GenuineLLR()

    function cn_update(a::GenuineLLR, b::GenuineLLR)
      return GenuineLLR( 2*atanh( tanh(a.v / 2) * tanh(b.v / 2) ) )
    end

    function vn_update(a::GenuineLLR, b::GenuineLLR, c::GF2Element)
        if c == GF2_0
            v =  a.v + b.v
        else # c == GF2_1
            v = -a.v + b.v
        end

        if isnan(v)
            @assert isinf(a.v) && isinf(b.v)
            v = 0.0
        end

        @assert !isnan(v)
        return GenuineLLR(v)
    end

    function -(a::GenuineLLR)
        return GenuineLLR(-a.v)
    end
    
    function val(a::GenuineLLR)
        return a.v
    end
    
    label(::Type{GenuineLLR}) = "GenuineLLR"

    

    struct MinApproxLLR <: LLRType
        v::Float64
    end

    MinApproxLLR() = MinApproxLLR(0.00)
    zero(::Type{MinApproxLLR}) = MinApproxLLR()

    function cn_update(a::MinApproxLLR, b::MinApproxLLR)
        return MinApproxLLR(sign(a.v)*sign(b.v)*min(abs(a.v), abs(b.v)))
    end

    function vn_update(a::MinApproxLLR, b::MinApproxLLR, c::GF2Element)
        if c == GF2_0
            v =  a.v + b.v
        else # c == GF2_1
            v = -a.v + b.v
        end

        @assert !isnan(v)
        return MinApproxLLR(v)
    end

    function -(a::MinApproxLLR)
        return MinApproxLLR(-a.v)
    end
    
    function val(a::MinApproxLLR)
        return a.v
    end
    
    label(::Type{MinApproxLLR}) = "MinApproxLLR"
    


    struct Q3lLLR <: LLRType
        v::UInt8

        function Q3lLLR(v::UInt8)
            # @assert v in (1, 2, 3)
            return new(v)
        end
    end

    const Q3lΔ = 1.0
    # const Q3lpmincm∞ = ln(1.0 + exp( Q3lΔ))
    # const Q3lpminc00 = ln(2.0)
    # const Q3lpmincp∞ = ln(1.0 + exp(-Q3lΔ))
    const Q3lm∞ = Q3lLLR(UInt8(1))
    const Q3l00 = Q3lLLR(UInt8(2))
    const Q3lp∞ = Q3lLLR(UInt8(3))

    # export Q3lm∞, Q3l00, Q3lp∞

    Q3lLLR() = Q3l00
    zero(::Type{Q3lLLR}) = Q3lLLR()

    Q3lLLR(v::AbstractFloat) = abs(v) < 1e-10 ? Q3l00 : v > 0.0 ? Q3lp∞ : v < 0.0 ? Q3lm∞ : Q3l00

    const cn_update_lut = SMatrix{3,3}([
        # Q3lm∞ Q3l00 Q3lp∞
          Q3lp∞ Q3l00 Q3lm∞ ;  # Q3lm∞
          Q3l00 Q3l00 Q3l00 ;  # Q3l00
          Q3lm∞ Q3l00 Q3lp∞ ;  # Q3lp∞
    ])
    const vn_update_lut = SMatrix{3,3}([
        # Q3lm∞ Q3l00 Q3lp∞
          Q3lm∞ Q3lm∞ Q3l00 ;  # Q3lm∞
          Q3lm∞ Q3l00 Q3lp∞ ;  # Q3l00
          Q3l00 Q3lp∞ Q3lp∞ ;  # Q3lp∞
    ])
    const vn_update_lut_0 = SMatrix{3,3}([ vn_update_lut[a.v, b.v] for a in (Q3lm∞, Q3l00, Q3lp∞), b in (Q3lm∞, Q3l00, Q3lp∞) ])
    const vn_update_lut_1 = SMatrix{3,3}([ vn_update_lut[4-a.v, b.v] for a in (Q3lm∞, Q3l00, Q3lp∞), b in (Q3lm∞, Q3l00, Q3lp∞) ])

    function cn_update(a::Q3lLLR, b::Q3lLLR)
        @inbounds begin
            return cn_update_lut[a.v, b.v]
        end
    end

    function vn_update(a::Q3lLLR, b::Q3lLLR, c::GF2Element)
        @inbounds begin
            if c == GF2_0
                return vn_update_lut_0[a.v, b.v]
            else # c == GF2_1
                return vn_update_lut_1[a.v, b.v]
            end
        end
    end

    function -(a::Q3lLLR)
        if a == Q3lm∞
            return Q3lp∞
        elseif a == Q3lp∞
            return Q3lm∞
        else
            return Q3l00
        end
    end

    function val(a::Q3lLLR)
        if a == Q3lm∞
            return -Q3lΔ
        elseif a == Q3lp∞
            return +Q3lΔ
        else
            return 0.0
        end
    end
    
    label(::Type{Q3lLLR}) = "Q3lLLR"



    # PMs

    import Base: isless, zero
    using Polar.PolarDecoding: PMType, LLRType
    using Polar.GF2n: GF2Element, GF2_0, GF2_1
    import Polar.PolarDecoding: update
    using Polar.CommunicationsUtils

    export GenuinePM, TwoPWLPM, ThreePWLPM

    
    struct GenuinePM <: PMType
        v::Float64
    end

    GenuinePM() = GenuinePM(0.00)
    
    
    struct TwoPWLPM <: PMType
        v::Float64
    end

    TwoPWLPM() = TwoPWLPM(0.00)
    
    
    struct ThreePWLPM <: PMType
        v::Float64
    end

    ThreePWLPM() = ThreePWLPM(0.00)
    
    
    
    zero(::Type{Tpm}) where {Tpm<:Union{GenuinePM,TwoPWLPM,ThreePWLPM}} = Tpm()

    function update(a::Tpm, b::Tllr, c::GF2Element, i::Int) where {Tllr<:LLRType,Tpm<:Union{GenuinePM,TwoPWLPM,ThreePWLPM}}
        if c == GF2_0
            return Tpm(a.v + pm_increment(Tpm,  b, i))
        else # c == GF2_1
            return Tpm(a.v + pm_increment(Tpm, -b, i))
        end
    end

    function isless(x::Tpm, y::Tpm) where {Tpm<:Union{GenuinePM,TwoPWLPM,ThreePWLPM}}
        return isless(x.v, y.v)
    end

    
    
    function pm_increment(::Type{GenuinePM}, b::Tllr, i::Int) where {Tllr<:LLRType}
        return ln(1 + exp(-val(b)))
    end
    
    function pm_increment(::Type{TwoPWLPM}, b::Tllr, i::Int) where {Tllr<:LLRType}
        v = val(b)
        if v >= 0.0
            return 0.0
        else
            return -v
        end
    end
    
    function pm_increment(::Type{ThreePWLPM}, b::Tllr, i::Int) where {Tllr<:LLRType}
        Δ = 2*ln(2)
        λ = val(b)
        if λ <= -Δ
            return -λ
        elseif λ <= +Δ
            return -(λ-Δ)/2
        else
            return 0.0
        end
    end
    
    
    
    label(::Type{GenuinePM}) = "GenuinePM"
    label(::Type{TwoPWLPM}) = "TwoPWLPM"
    label(::Type{ThreePWLPM}) = "ThreePWLPM"

end
