__precompile__()

module MySimulationUtils

	using Random
	using Distributed

	using Polar.PolarCodes
	using Polar.PolarDecoding
	using Polar.Channels
	using Polar.GF2n
	using Polar.CommunicationsUtils
	using Polar.SimulationUtils


	# this helper is necessary for the likelihood computations required for the ML-LB

	function loglikelihood(ch_w_0, ch_w_1, ys, cs)
		n = length(ys)
		@assert length(ys) == n && length(cs) == n

		v = 0.0
		for i in 1:n
			if cs[i] == GF2_0
				v += ln(ch_w_0(ys[i]))
			else # if cs[i] == GF2_1
				v += ln(ch_w_1(ys[i]))
			end
		end
		return v
	end


	# some helper functions, not meant to be used outside this module ...

	# function sim_helper_prepare_channel_functions_02(sim_channel::T_channel, ::Type{sim_T_llr}) where {sim_T_llr<:LLRType,T_channel<:Union{BiAWGNChannel,Q3LBiAWGNChannel,OptQ3LBiAWGNChannel}}
	function sim_helper_prepare_channel_functions_02(sim_channel::T_channel, ::Type{sim_T_llr}) where {sim_T_llr<:LLRType,T_channel}
		@assert intype(T_channel) == BinaryAlphabet

	    ch_f = get_f(sim_channel)

	    ch_llr = get_llr(sim_channel)
	    ch_y2llr = y -> sim_T_llr(ch_llr(y))

	    ch_w_0 = get_w(sim_channel, GF2_0)
	    ch_w_1 = get_w(sim_channel, GF2_1)

	    return (ch_f, ch_y2llr, ch_w_0, ch_w_1)
	end


	# simulation experiment: simple experiment with PM error, list error and ML-LB

	export sim_parallel_experiments_fer_simple01_refined, sim_local_experiments_fer_simple01_refined

	function sim_local_experiments_fer_simple01_refined(sim_code, sim_decoder, sim_channel; N_iterations=100)
		sim_T_llr = llrtype(sim_decoder)
		sim_T_pm = pmtype(sim_decoder)

		(tmp_u, tmp_u_full, tmp_c, tmp_y, tmp_llrs, tmp_ĉ, tmp_û) = SimulationUtils.sim_helper_prepare_tmp_variables_01(sim_code, sim_channel, sim_T_llr)
		(ch_f, ch_y2llr, ch_w_0, ch_w_1) = sim_helper_prepare_channel_functions_02(sim_channel, sim_T_llr)

		result = SimulationResult(:PMError, :ListError, :MLError)

		for i_iteration in 1:N_iterations

			# count a new :ANY
			tic!(result)

			# generate dummy data and polar-encode
			Random.rand!(sim_code, tmp_u)
			polarencode!(sim_code, tmp_c, tmp_u)

			# simulate channel
			map!(ch_f, tmp_y, tmp_c)
			map!(ch_y2llr, tmp_llrs, tmp_y)

			# decode (SCL)
			decode!(sim_decoder, tmp_llrs)
			pm_opt = extract_best_pm_C!(sim_decoder, tmp_ĉ)


			# ERROR EVALUATION

			# PMError

			if tmp_c != tmp_ĉ
				tic!(result, :PMError)
			end


			# ListError

			(list_cs, list_pms) = extract_list_C(sim_decoder)
			tmp_list_error_has_happened = false
			if !any(l -> (list_cs[:,l] == tmp_c), 1:length(list_pms))
				tic!(result, :ListError)
				tmp_list_error_has_happened = true
			end


			# ML LB

			ll_true = loglikelihood(ch_w_0, ch_w_1, tmp_y, tmp_c)
			list_lls = [ loglikelihood(ch_w_0, ch_w_1, tmp_y, list_cs[:,l]) for l in 1:length(list_pms) ]

			if ll_true > maximum(list_lls)
				# we can't say for sure whether the ML would have made a mistake or not, the true codeword has better likelihood than all of the list at least (this is why ML lower BOUND)

			elseif ll_true == maximum(list_lls)
				# some list entries are as good as the true codeword, picking randomly among them, the ML might have made mistakes in some cases
				competing_paths = list_lls .== ll_true
				if tmp_list_error_has_happened
					n_competing_paths = count(competing_paths) + 1
				else
					n_competing_paths = count(competing_paths) + 0
				end
				if rand() > 1. / n_competing_paths
					tic!(result, :MLError)
				end

			else # if ll_true < maximum(list_lls)
				# the ML would have made a mistake
				tic!(result, :MLError)

			end

		end

		return result
	end

	function sim_parallel_experiments_fer_simple01_refined(sim_code, sim_decoder, sim_channel; desired_ci_rel=0.1, desired_ci_event=:PMError, N_jobsperbatch=96, N_framesperjob=200)
		result = SimulationResult(:PMError, :ListError, :MLError)

		@timev while ci_wilson_rel(result, desired_ci_event) > desired_ci_rel
			result = result + reduce(+, pmap(i -> sim_local_experiments_fer_simple01_refined(sim_code, sim_decoder, sim_channel; N_iterations=N_framesperjob), 1:N_jobsperbatch))
		end

		return result
	end



	# simulation experiment: experiment with ML-in-list

	export sim_parallel_experiments_fer_mlinlist01_refined, sim_local_experiments_fer_mlinlist01_refined

	function sim_local_experiments_fer_mlinlist01_refined(sim_code, sim_decoder, sim_channel; N_iterations=100)
		sim_T_llr = llrtype(sim_decoder)
		sim_T_pm = pmtype(sim_decoder)

		(tmp_u, tmp_u_full, tmp_c, tmp_y, tmp_llrs, tmp_ĉ, tmp_û) = SimulationUtils.sim_helper_prepare_tmp_variables_01(sim_code, sim_channel, sim_T_llr)
		(ch_f, ch_y2llr, ch_w_0, ch_w_1) = sim_helper_prepare_channel_functions_02(sim_channel, sim_T_llr)

		# channel function custom to ML-in-list decoding
		ch_llr = get_llr(sim_channel)

		# variables to set up ML-in-list decoding variables
		n = PolarDecoding.n(sim_decoder)
		L = PolarDecoding.L(sim_decoder)

		# custom variables for ML-in-list decoding
		tmp_cs = Matrix{GF2Element}(undef, n, L)
		tmp_pms = Vector{sim_T_pm}(undef, L)
		tmp_float_llrs = Vector{Float64}(undef, n)
		tmp_likelihoods = Vector{Float64}(undef, L)

		result = SimulationResult(:MLInListError)

		for i_iteration in 1:N_iterations

			# count a new :ANY
			tic!(result)

			# generate dummy data and polar-encode
			rand!(sim_code, tmp_u)
			polarencode!(sim_code, tmp_c, tmp_u)

			# simulate channel
			map!(ch_f, tmp_y, tmp_c)
			map!(ch_y2llr, tmp_llrs, tmp_y)

			# decode (SCL)
			L_used = decode!(sim_decoder, tmp_llrs)
			extract_list_C!(sim_decoder, tmp_cs, tmp_pms)
			@assert L_used == L


			# ERROR EVALUATION

			# calculate likelihoods of list entries
			for l in 1:L
				tmp_likelihoods[l] = loglikelihood(ch_w_0, ch_w_1, tmp_y, tmp_cs[:,l])
			end
			max_val, max_idx = findmax(tmp_likelihoods)
			copy!(tmp_ĉ, tmp_cs[:,max_idx])

			if tmp_c != tmp_ĉ
				tic!(result, :MLInListError)
			end

		end

		return result
	end

	function sim_parallel_experiments_fer_mlinlist01_refined(sim_code, sim_decoder, sim_channel; desired_ci_rel=0.1, desired_ci_event=:MLInListError, N_jobsperbatch=96, N_framesperjob=200)
		result = SimulationResult(:MLInListError)

		@timev while ci_wilson_rel(result, desired_ci_event) > desired_ci_rel
			result = result + reduce(+, pmap(i -> sim_local_experiments_fer_mlinlist01_refined(sim_code, sim_decoder, sim_channel; N_iterations=N_framesperjob), 1:N_jobsperbatch))
		end

		return result
	end

end
