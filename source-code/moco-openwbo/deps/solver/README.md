# HOW TO BUILD 

Run

```
$ mvn -DskipTests=true package
```

to generate the
./target/org.sat4j.threeAlgorithms-0.0.1-SNAPSHOT-jar-with-dependencies.jar
jar file.

# HOW TO RUN

```
java -jar  ./target/org.sat4j.moco.threeAlgorithms-0.0.1-SNAPSHOT-jar-with-dependencies.jar -alg 2 examples/example1.opb
```

(Never used, but should work)
# HOW TO RUN THE ANALYZER

The org.sat4j.moco.jar includes an analysis tool for evaluating the
quality of Pareto front approximations through indicators such as
hypervolume and inverted generational distance.

Run it using the following command:

```
$ java -cp org.sat4j.moco.jar org.sat4j.moco.analysis.Analyzer <instance file> [<label>:<output file>]+
```

`<output file>` is expected to be in the output format produced by the
MOCO solver.  If there exist multiple files for different runs of the
same algorithm, these should have the same `<label>`.

# HOW TO RUN THE TRANSLATOR

To print the encoding of the logical circuit of the selection
delimeter into ``out.mocnf``,


```
java -jar  ./target/org.sat4j.moco.threeAlgorithms-0.0.1-SNAPSHOT-jar-with-dependencies.jar  in.opb -ib 1 -o out.mocnf
```

If you would like to enforce upper limits for the iObj'th objective function,
use the ``ul`` option with the syntax ``<iObj>:<upperLimit>``,
where ``upperLimit`` is the exclusive upper limit, and ``iObj`` is the
index of the objective function, starting from 0. Uppper limits may be
set for any subset of the objective functions, separating the upper limits with a comma. 
To impose an upper limit of 1 to the objective function 0 and an upper limit of 2 to the objective function 4,
use ``-ul 0:1,4:2``

To set the ratios used in the digital basis, use the option ``rb <arg>``, with ``<arg>=<iObj>:[r1,...rN],...``
For exmple, to set the ratios of the first objective function `1` to ``[2,2,3]`` and the ratios of the objective function ``3`` to ``[4,4]``,
use ``-rb 0:[2,3,3],3:[4,4]``.

The last ratio will be inserted into the end of the list as many times
as needed to fully represent the attainable values.


NOTE: Sequential encoder is not working.
