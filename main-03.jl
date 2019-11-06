using Distributed

@everywhere begin
    using Pkg
    Pkg.activate(".")
end

@everywhere begin
	EXID = "main-03"

	using Random
	Random.seed!(2342+myid())
end

@everywhere begin
	include("MyChannels.jl")
	include("MyLLRsPMs.jl")
	include("MySimulationUtils.jl")

	include("MyEnhancedQ3lTypes.jl")
end

@everywhere begin
	using Polar.PolarDecoding
	using Polar.PolarCodes
	using Polar.Channels
	using Polar.GF2n
	using Polar.CommunicationsUtils
	using Polar.SimulationUtils
	using Polar.Utils

	using .MyChannels
	using .MyLLRsPMs
	using .MySimulationUtils

	using .MyEnhancedQ3lTypes
end


# codes

# n=256
# DE: 3L-Q-BiAWGN + 3L-Q-LLRs (from 001/31-...)
const C_own_Q3LBiAWGN_Q3LLLR_EbN045_001_31 = PolarCode([true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, false, true, true, true, true, true, true, true, false, true, true, true, false, true, false, false, false, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, false, true, true, true, true, true, true, true, false, true, true, true, false, true, false, false, false, true, true, true, true, true, true, true, false, true, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, false, true, true, true, true, true, true, false, false, true, false, false, false, false, false, false, false, true, true, true, false, true, false, false, false, true, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, true, false, true, false, false, false, true, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false])


# simulation scenarios

const scenarios = [
	# ASILOMAR PAPER
	("q3biawgn-q3llrs-C_own_Q3LBiAWGN_Q3LLLR_EbN045_001_31", C_own_Q3LBiAWGN_Q3LLLR_EbN045_001_31, [32,], 1.00:0.25:5.25, OptQ3LBiAWGNChannel, Q3lLLR, ThreePWLPM),
	("q3biawgn-epmuq3llr01s-C_own_Q3LBiAWGN_Q3LLLR_EbN045_001_31", C_own_Q3LBiAWGN_Q3LLLR_EbN045_001_31, [32,], 1.00:0.25:5.00, OptQ3LBiAWGNChannel, EPMUQ3lLLR01, CustomPMUpdatePM),
]


# simulations

for scenario in scenarios

	scenario_name = scenario[1]
	scenario_code = scenario[2]
	scenario_Ls = scenario[3]
	scenario_EbN0dBs = scenario[4]
	scenario_channel = scenario[5]
	scenario_llrtype = scenario[6]
	scenario_pmtype = scenario[7]

	println("SCENARIO: $(scenario_name)")

	for L in scenario_Ls

		println("\t# L = $(L)")
		EXID_SUBID = "$(scenario_name)-L$(L)"
		SIM_DATA = Dict{Float64,SimulationResult{(:MLInListError,)}}()

		sim_code = scenario_code
		sim_decoder = SCLDecoder{Int(log2(PolarCodes.n(sim_code))),L,scenario_llrtype,scenario_pmtype,sim_code}()

		for EbN0dB in scenario_EbN0dBs
			σ² = esn0db2σ²(ebn0db2esn0db(EbN0dB, PolarCodes.R(sim_code)))
			sim_channel = scenario_channel(σ²)

			SIM_DATA[EbN0dB] = sim_parallel_experiments_fer_mlinlist01_refined(sim_code, sim_decoder, sim_channel; desired_ci_rel=0.05, desired_ci_event=:MLInListError)

			println("### EbN0=$(EbN0dB)dB: $(SIM_DATA[EbN0dB])")
			gnuplotexport("$(EXID)-$(EXID_SUBID).dat", "EbN0dB", SIM_DATA)
		end

	end

end
