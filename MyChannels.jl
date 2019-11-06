__precompile__()

module MyChannels

	using SpecialFunctions
	using Polar.GF2n: GF2Element, GF2_0, GF2_1, toint
	using Polar.CommunicationsUtils
	using Polar.Channels
	import Polar.Channels: intype, outtype, get_w, get_f, get_llr


	export SeptenaryM3, SeptenaryM2, SeptenaryM1, Septenary0, SeptenaryP1, SeptenaryP2, SeptenaryP3
	export Q3LBiAWGNChannel, OptQ3LBiAWGNChannel
	export Q7LBiAWGNChannel, OptQ7LBiAWGNChannel
	

	@enum SeptenarySymbol SeptenaryM3=1 SeptenaryM2=2 SeptenaryM1=3 Septenary0=4 SeptenaryP1=5 SeptenaryP2=6 SeptenaryP3=7
	const SeptenaryAlphabet = SeptenarySymbol


	struct Q3LBiAWGNChannel <: AbstractCommunicationsChannel
		σ²::Float64
		δ::Float64
		_p::Float64
		_e::Float64
		_notpe::Float64

		function Q3LBiAWGNChannel(σ²::Float64, δ::Float64)
			F = x′ -> 1/2*(1+erf( (x′-mod_bpsk(GF2_0))/sqrt(2*σ²) ))
			_p = F(-δ)
			_e = F(+δ) - _p
			_notpe = 1-_p-_e
			return new(σ², δ, _p, _e, _notpe)
		end
	end
	intype(::Type{Q3LBiAWGNChannel}) = BinaryAlphabet
	outtype(::Type{Q3LBiAWGNChannel}) = TernaryAlphabet
	function get_w(ch::Q3LBiAWGNChannel, x::BinaryAlphabet)
		function w(y::TernaryAlphabet)::Float64
			if y == TernaryPΔ
				return x == GF2_0 ? ch._notpe : ch._p
			elseif y == Ternary0
				return ch._e
			else # y == TernaryMΔ
				return x == GF2_0 ? ch._p : ch._notpe
			end
		end
		return w
	end
	function get_f(ch::Q3LBiAWGNChannel)
	    _p_0 = ch._p
	    _pe_0 = ch._p + ch._e
	    _p_1 = ch._notpe
	    _pe_1 = ch._notpe + ch._e
		function f(x::BinaryAlphabet)::TernaryAlphabet
			p = x == GF2_0 ? _p_0 : _p_1
			pe = x == GF2_0 ? _pe_0 : _pe_1
			r = rand()
			if r <= p
				return TernaryMΔ
			elseif r <= pe
				return Ternary0
			else
				return TernaryPΔ
			end
		end
		return f
	end
	function get_llr(ch::Q3LBiAWGNChannel)
		Δ = ln(ch._notpe/ch._p)
		function llr(y::TernaryAlphabet)
			if y == TernaryMΔ
				return -Δ
			elseif y == TernaryPΔ
				return +Δ
			else
				return 0.0
			end
		end
		return llr
	end


	struct OptQ3LBiAWGNChannel <: AbstractCommunicationsChannel
		ch::Q3LBiAWGNChannel
		_δ::Float64

		function OptQ3LBiAWGNChannel(σ²::Float64)
			Plog2Q(p, q) = p > 1e-300 ? p*log2(q) : 0.0
			# http://www.lntwww.de/cgi-bin/extern/uni.pl?uno=hyperlink&due=block&b_id=3264&hyperlink_typ=block_verweis&hyperlink_fenstergroesse=blockverweis_gross&session_id=
			C_BSEC(p, e) = Plog2Q(1-p-e, (1-p-e)/(1/2*(1-e))) + Plog2Q(p, p/(1/2*(1-e)))
			function Q3LBiAWGNCapacity_getδopt_helper(σ²; Δ=0.01, lims=(0, 5))
			    @assert length(lims) == 2
			    Rbest = 0.0
			    δbest = 0.0
			    F = x -> 1/2*(1+erf( (x-1)/sqrt(2*σ²) ))
			    
			    for δ in lims[1]:Δ:lims[2]
			        p = F(-δ)
			        e = F(+δ) - p
			        R = C_BSEC(p, e)
			        if R >= Rbest
			            Rbest = R
			            δbest = δ
			        end
			    end
			    
			    return (Rbest, δbest)
			end

			_δ = Q3LBiAWGNCapacity_getδopt_helper(σ²)[2]

			ch = Q3LBiAWGNChannel(σ², _δ)
			return new(ch, _δ)
		end
	end
	intype(::Type{OptQ3LBiAWGNChannel}) = BinaryAlphabet
	outtype(::Type{OptQ3LBiAWGNChannel}) = TernaryAlphabet
	get_w(ch::OptQ3LBiAWGNChannel, x::BinaryAlphabet) = get_w(ch.ch, x)
	get_f(ch::OptQ3LBiAWGNChannel) = get_f(ch.ch)
	get_llr(ch::OptQ3LBiAWGNChannel) = get_llr(ch.ch)


	struct Q7LBiAWGNChannel <: AbstractCommunicationsChannel
		σ²::Float64
		δ::Float64
		_pM3::Float64
		_pM2::Float64
		_pM1::Float64
		_p0::Float64
		_pP1::Float64
		_pP2::Float64
		_pP3::Float64

		function Q7LBiAWGNChannel(σ²::Float64, δ::Float64)
			F = x′ -> 1/2*(1+erf( (x′-mod_bpsk(GF2_0))/sqrt(2*σ²) ))
			_pM3 = F(-5*δ)
			_pM2 = F(-3*δ) - F(-5*δ)
			_pM1 = F(-1*δ) - F(-3*δ)
			_p0 = F(+1*δ) - F(-1*δ)
			_pP1 = F(+3*δ) - F(+1*δ)
			_pP2 = F(+5*δ) - F(+3*δ)
			_pP3 = 1 - F(+5*δ)
			@assert _pM3+_pM2+_pM1+_p0+_pP1+_pP2+_pP3 ≈ 1.0

			return new(σ², δ, _pM3, _pM2, _pM1, _p0, _pP1, _pP2, _pP3)
		end
	end
	intype(::Type{Q7LBiAWGNChannel}) = BinaryAlphabet
	outtype(::Type{Q7LBiAWGNChannel}) = SeptenaryAlphabet
	function get_w(ch::Q7LBiAWGNChannel, x::BinaryAlphabet)
		function w(y::SeptenaryAlphabet)::Float64
			if y == SeptenaryM3
				return x == GF2_0 ? ch._pM3 : ch._pP3
			elseif y == SeptenaryM2
				return x == GF2_0 ? ch._pM2 : ch._pP2
			elseif y == SeptenaryM1
				return x == GF2_0 ? ch._pM1 : ch._pP1
			elseif y == Septenary0
				return ch._p0
			elseif y == SeptenaryP1
				return x == GF2_0 ? ch._pP1 : ch._pM1
			elseif y == SeptenaryP2
				return x == GF2_0 ? ch._pP2 : ch._pM2
			else #if y == SeptenaryP3
				return x == GF2_0 ? ch._pP3 : ch._pM3
			end
		end
		return w
	end
	function get_f(ch::Q7LBiAWGNChannel)
		_p_cum_M3 = ch._pM3
		_p_cum_M2 = _p_cum_M3 + ch._pM2
		_p_cum_M1 = _p_cum_M2 + ch._pM1
		_p_cum_0 = _p_cum_M1 + ch._p0
		_p_cum_P1 = _p_cum_0 + ch._pP1
		_p_cum_P2 = _p_cum_P1 + ch._pP2
		_p_cum_P3 = _p_cum_P2 + ch._pP3
		@assert _p_cum_P3 ≈ 1.0

		function f(x::BinaryAlphabet)::SeptenaryAlphabet
			r = rand()
			if r <= _p_cum_M3
				return x == GF2_0 ? SeptenaryM3 : SeptenaryP3
			elseif r <= _p_cum_M2
				return x == GF2_0 ? SeptenaryM2 : SeptenaryP2
			elseif r <= _p_cum_M1
				return x == GF2_0 ? SeptenaryM1 : SeptenaryP1
			elseif r <= _p_cum_0
				return Septenary0
			elseif r <= _p_cum_P1
				return x == GF2_0 ? SeptenaryP1 : SeptenaryM1
			elseif r <= _p_cum_P2
				return x == GF2_0 ? SeptenaryP2 : SeptenaryM2
			else #if r <= _p_cum_P3
				return x == GF2_0 ? SeptenaryP3 : SeptenaryM3
			end
		end
		return f
	end
	function get_llr(ch::Q7LBiAWGNChannel)
		_llr_M3 = ln(ch._pM3/ch._pP3)
		_llr_M2 = ln(ch._pM2/ch._pP2)
		_llr_M1 = ln(ch._pM1/ch._pP1)
		_llr_0 = ln(1.0) #ln(ch._p0/ch._p0)
		_llr_P1 = ln(ch._pP1/ch._pM1)
		_llr_P2 = ln(ch._pP2/ch._pM2)
		_llr_P3 = ln(ch._pP3/ch._pM3)

		function llr(y::SeptenaryAlphabet)
			if y == SeptenaryM3
				return _llr_M3
			elseif y == SeptenaryM2
				return _llr_M2
			elseif y == SeptenaryM1
				return _llr_M1
			elseif y == Septenary0
				return _llr_0
			elseif y == SeptenaryP1
				return _llr_P1
			elseif y == SeptenaryP2
				return _llr_P2
			else #if y == SeptenaryP3
				return _llr_P3
			end
		end
		return llr
	end


	struct OptQ7LBiAWGNChannel <: AbstractCommunicationsChannel
		ch::Q7LBiAWGNChannel
		_δ::Float64

		function OptQ7LBiAWGNChannel(σ²::Float64)
			function I_Q7LBiAWGNC(ch::Q7LBiAWGNChannel)
				W_0 = get_w(ch, GF2_0)
				W_1 = get_w(ch, GF2_1)

				function p_X(x::BinaryAlphabet)
					return 0.5
				end

				function p_XY(x::BinaryAlphabet, y::SeptenaryAlphabet)
					if x == GF2_0
						return p_X(x) * W_0(y)
					else #if x == GF2_1
						return p_X(x) * W_1(y)
					end
				end

				function p_Y(y::SeptenaryAlphabet)
					return p_XY(GF2_0, y) + p_XY(GF2_1, y)
				end

				alphabet_X = [GF2_0, GF2_1]
				alphabet_Y = [SeptenaryM3, SeptenaryM2, SeptenaryM1, Septenary0, SeptenaryP1, SeptenaryP2, SeptenaryP3]

				return sum([ Plog2Q(p_XY(x, y), p_XY(x, y)/(p_X(x)*p_Y(y))) for x in alphabet_X, y in alphabet_Y ])
			end

			function Q7LBiAWGNCapacity_getδopt_helper(σ²; Δ=0.01, lims=(0, 5))
			    @assert length(lims) == 2
			    Rbest = -1.0
			    δbest = 0.0
			    
			    for δ in lims[1]:Δ:lims[2]
			    	_ch = Q7LBiAWGNChannel(σ², δ)
			    	R = I_Q7LBiAWGNC(_ch)
			        if R >= Rbest
			            Rbest = R
			            δbest = δ
			        end
			    end
			    
			    @assert δbest != lims[2]
			    @assert Rbest >= 0.0

			    return (Rbest, δbest)
			end

			_δ = Q7LBiAWGNCapacity_getδopt_helper(σ²)[2]

			ch = Q7LBiAWGNChannel(σ², _δ)
			return new(ch, _δ)
		end
	end
	intype(::Type{OptQ7LBiAWGNChannel}) = BinaryAlphabet
	outtype(::Type{OptQ7LBiAWGNChannel}) = SeptenaryAlphabet
	get_w(ch::OptQ7LBiAWGNChannel, x::BinaryAlphabet) = get_w(ch.ch, x)
	get_f(ch::OptQ7LBiAWGNChannel) = get_f(ch.ch)
	get_llr(ch::OptQ7LBiAWGNChannel) = get_llr(ch.ch)

end
