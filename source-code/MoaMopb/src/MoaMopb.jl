module MoaMopb

using Printf

import HiGHS
import SCIP
import Gurobi
import CPLEX

import MathOptInterface as MOI
import MultiObjectiveAlgorithms as MOA

function remove_empty(toks)
    r = findall(x -> x == "", toks)
    deleteat!(toks, r)
end

function ensure_var(
    name::String,
    vars::Vector{MOI.VariableIndex},
    model::MOA.Optimizer,
)::MOI.VariableIndex
    @assert(startswith(name, 'x'))
    idx = parse(UInt, name[2:length(name)])
    if length(vars) < idx
        old_len = length(vars)
        append!(vars, MOI.add_variables(model, Int(idx - length(vars))))
        for i::UInt = old_len+1:idx
            MOI.add_constraint(model, vars[i], MOI.ZeroOne())
        end
    end
    @assert(length(vars) >= idx)
    return vars[idx]
end

function parse_func(
    toks::Vector{SubString{String}},
    vars::Vector{MOI.VariableIndex},
    model::MOA.Optimizer,
)::MOI.ScalarAffineFunction{Float64}
    func::Vector{MOI.ScalarAffineTerm{Float64}} = []
    offset::Float64 = 0
    for i::UInt = 1:length(toks)/2
        tok = toks[i*2]
        name = tok
        coef = parse(Float64, toks[i*2-1])
        if startswith(tok, "~x")
            name = tok[2:length(tok)]
            offset += coef
            coef = -coef
        else
            @assert(startswith(toks[i*2], 'x'))
        end
        push!(func, MOI.ScalarAffineTerm(coef, ensure_var(string(name), vars, model)))
    end
    return MOI.ScalarAffineFunction(func, offset)
end

function run(
    path::AbstractString,
    opt::Type,
    alg::MOA.AbstractAlgorithm,
    tlim::Union{Int,Nothing},
    print_sol::Bool,
    print_model::Bool,
    silent::Bool,
    threads::Union{Int,Nothing},
    tuned::Bool,
)
    model = MOA.Optimizer(opt)
    MOI.set(model, MOI.Silent(), silent)
    MOI.set(model, MOA.Algorithm(), alg)
    MOI.set(model, MOI.TimeLimitSec(), tlim)
    if tuned
        MOI.set(model, MOI.AbsoluteGapTolerance(), 1e-6)
        MOI.set(model, MOI.RelativeGapTolerance(), 0.0)
        if opt == CPLEX.Optimizer
            MOI.set(model, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Tolerances_Integrality"), 1e-9)
        elseif opt == HiGHS.Optimizer
            MOI.set(model, MOI.RawOptimizerAttribute("mip_feasibility_tolerance"), 1e-8)
        end
    end
    if !isnothing(threads)
        if opt == SCIP.Optimizer
            MOI.set(model, MOI.RawOptimizerAttribute("lp/threads"), threads)
            MOI.set(model, MOI.RawOptimizerAttribute("parallel/maxnthreads"), threads)
        else
            MOI.set(model, MOI.NumberOfThreads(), threads)
        end
    end

    println(
        "===[ Solver Info ]==============================================================",
    )
    println(
        MOI.get(model, MOI.SolverName()),
        " (",
        MOI.get(model, MOI.SolverVersion()),
        ")",
    )
    println("Time limit: ", tlim, "s")
    println("Threads: ", threads)
    println("Tuned: ", tuned)

    file = open(path, "r")

    objterms::Vector{MOI.VectorAffineTerm{Float64}} = []
    objoffsets::Vector{Float64} = []
    vars::Vector{MOI.VariableIndex} = []

    while !eof(file)
        raw_line = readline(file)
        if startswith(raw_line, '*')
            # comment
            continue
        end
        line = string(strip(raw_line))
        if isempty(line)
            # empty line
            continue
        end
        @assert(line[length(line)] == ';')
        toks = split(line[1:length(line)-1], " ")
        toks = remove_empty(toks)
        if toks[1] == "min:"
            # objective
            func = parse_func(toks[2:length(toks)], vars, model)
            push!(objoffsets, func.constant)
            idx = length(objoffsets)
            append!(objterms, [MOI.VectorAffineTerm(idx, term) for term in func.terms])
            continue
        end
        # constraint
        bound = parse(Float64, toks[length(toks)])
        func = parse_func(toks[1:length(toks)-2], vars, model)
        if toks[length(toks)-1] == ">="

            MOI.Utilities.normalize_and_add_constraint(model, func, MOI.GreaterThan(bound))
        else
            @assert(toks[length(toks)-1] == "=")
            MOI.Utilities.normalize_and_add_constraint(model, func, MOI.EqualTo(bound))
        end
    end

    obj = MOI.VectorAffineFunction(objterms, objoffsets)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.VectorAffineFunction{Float64}}(), obj)

    if print_model
        println(
            "===[ Model ]====================================================================",
        )
        print(model)
    end

    Base.exit_on_sigint(false)
    MOI.optimize!(model)

    println(
        "===[ Pareto Front ]=============================================================",
    )
    for i = 1:MOI.get(model, MOI.ResultCount())
        println(
            "---[ Non-dominated point ]------------------------------------------------------",
        )
        print("o")
        objvals = MOI.get(model, MOI.ObjectiveValue(i))
        for j in eachindex(objvals)
            print(" ", round(Int, objvals[j]))
        end
        println()
        print("f")
        objvals = MOI.get(model, MOI.ObjectiveValue(i))
        for j in eachindex(objvals)
            print(" ", objvals[j])
        end
        println()
        if print_sol
            print("v")
            for j = 1:MOI.get(model, MOI.NumberOfVariables())
                var = "-x$(@sprintf("%d", j))"
                if MOI.get(model, MOI.VariablePrimal(i), MOI.VariableIndex(j)) > 0.5
                    var = var[2:length(var)]
                end
                print(" ", var)
            end
            println()
        end
    end

    println(
        "===[ Summary ]==================================================================",
    )
    println("Termination status: ", MOI.get(model, MOI.TerminationStatus()))
    println("# Results: ", MOI.get(model, MOI.ResultCount()))
    println("Solve time: ", MOI.get(model, MOI.SolveTimeSec()), "s")
end

function main()
    print_model = false
    print_sol = false
    silent = true
    tlim = nothing
    threads = nothing
    tuned = false

    alg = MOA.Lexicographic()
    opt = HiGHS.Optimizer

    # parse args
    for x in ARGS[1:length(ARGS)-1]
        # possible arguments go here
        if x == "--lexicographic"
            # this is the default as well
            alg = MOA.Lexicographic()
        elseif x == "--epsilon"
            alg = MOA.EpsilonConstraint()
        elseif x == "--hierarchical"
            alg = MOA.Hierarchical()
        elseif x == "--tamby-vanderpooten"
            alg = MOA.TambyVanderpooten()
        elseif x == "--kirlik-sayin"
            alg = MOA.KirlikSayin()
        elseif x == "--dominguez-rios"
            alg = MOA.DominguezRios()
        elseif x == "--highs"
            opt = HiGHS.Optimizer
        elseif x == "--scip"
            opt = SCIP.Optimizer
        elseif x == "--gurobi"
            opt = Gurobi.Optimizer
        elseif x == "--cplex"
            opt = CPLEX.Optimizer
        elseif startswith(x, "--timelimit=")
            tlim = parse(Int, x[13:length(x)])
        elseif x == "--print-model"
            print_model = true
        elseif x == "--print-solutions"
            print_sol = true
        elseif x == "--non-silent"
            silent = false
        elseif startswith(x, "--threads=")
            threads = parse(Int, split(x, "=")[2])
        elseif x == "--tuned"
            tuned = true
        else
            @warn "Unknown option" x
        end
    end

    if isempty(ARGS)
        @warn "No instance file provided, exiting"
    else
        run(ARGS[length(ARGS)], opt, alg, tlim, print_sol, print_model, silent, threads, tuned)
    end
end

end # module MoaMopb
