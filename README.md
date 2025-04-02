# MOA MOPB

A tool to solve multi-objective pseudo-Boolean optimization problems (0-1 ILP)
in OPB format with [MultiObjectiveAlgorithms.jl](https://github.com/jump-dev/MultiObjectiveAlgorithms.jl).

## Usage

The tool is inteded to be used either from the command line (by calling the
`main` function, see CLI flags below) or via the `run` function.

For better performance of the CLI script, you can precompile it with the
collowing commands in the Julia REPL:
```julia-repl
julia> using PackageCompiler
pkg> activate .
julia> create_sysimage(["MoaMopb"], sysimage_path="MoaMopb.so", precompile_execution_file="precompile.jl")
```
The `precompile.jl` script will call the tool with certain configurations so
that all of the code is precompiled. Feel free to adjust the precompile script
to your usecase, for example comment out the Gurobi configuration if you do not
have a Gurobi license.

Afterwards, call the tool as either of the following two lines:
```bash
julia --project -JMoaMopb.so -e "using MoaMopb; MoaMopb.main()" <ARGS>
julia --project -JMoaMopb.so -- run.jl <ARGS>
```

## Arguments

| CLI Flag            | Description                                                          |
| ------------------- | -------------------------------------------------------------------- |
| `--lecicographic`   | Run the `MOA.Lexicographic` algorithm (default)                      |
| `--epsilon`         | Run the `MOA.EpsilonConstraint` algorithm                            |
| `--hierarchical`    | Run the `MOA.Hierarchical` algorithm                                 |
| `--dominguez-rios`  | Run the `MOA.DominguezRios` algorithm                                |
| `--highs`           | Use the HiGHS optimizer as backaend (default)                        |
| `--scip`            | Use the SCIP optimizer as backend                                    |
| `--gurobi`          | Use the Gurobi optimizer as backend                                  |
| `--cplex`           | Use the CPLEX optimizer as backend                                   |
| `--print-model`     | Print the read model before solving                                  |
| `--print-solutions` | Print the found solutions and not just the objective values          |
| `--non-silent`      | Activate output from the solving algorithm                           |
| `--threads=n`       | Limit the underlying MIP solver to `n` threads                       |
| `--tuned`           | Set some options for the underlying MIP solver to tune for exactness |
