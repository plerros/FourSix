# FourSix
Acute triangulation with CGAL

## 1 How to

### 1.1 Compile
```
cd build
cmake ..
make
```

### 1.2 Run
```
./FourSix -i input.json
```

Or to save the results to a file:

```
./FourSix -i input.json -o output.json
```

In the folder *test/* there are *input.json* files from CGshop, the Kapodestrian Universiy of Athens, and from me. The CGshop `.json` have been modified with additional parameters to pick which algorithm and parameters will be used.

Note: If the UI is enabled, it will display edges outside the shell. Those are an implementation artifact, and aren't part of the JSON result.

### 1.3 Benchmark

The .md files produced by the scripts can be diffed with tools like meld.

## 2 Implementation

### 2.1 File structure

- `src/configuration.hpp`
	* Project's Configuration
	* Output / Debugging options

- `src/config_cgal.hpp`
    * CGAL Preprocessor defs

- `src/main.cpp`

- `src/json/json.cpp` / `.hpp`
    * Read / Write json

- `src/data/data_in.cpp` / `.hpp`
    * Convert json to c++ types
    * Input checking
    * Reduce infomation duplicates

- `src/data/data.cpp` / `.hpp`
    * Convert c++ types to CGAL types
    * Adds mesh shell to constraints
    * Removes elements outside of shell
    * Additional useful datastructs

- `src/data/data_out.cpp` / `.hpp`
    * Combines data and triangulation to produce `json::value` output
    * Filters lines outside of shell

- `src/obtuse_polygon/obtuse_polygon.cpp` / `.hpp`
    * Polygon construction from neighboring obtuse triangles
    * Methods to add triangles to polygon
    * Implementation limits us to 2 non-obtuse neighbors per triangle	

- `src/triangulation/triangulation.cpp` / `.hpp`
    * Simple intrerface for CDT
    * Obtuse triangle counting
    * Methods to add steiner points

- `src/optimization_methods/mixed_recursive.cpp` / `.hpp`
    * Optimization methods used by the problem space search algorithms

- `src/optimization_methods/local_search.cpp` / `.hpp`
    * Local search algorithm

- `src/optimization_methods/simulated_annealing.cpp` / `.hpp`
    * Simulated annealing algorithm

- `src/optimization_methods/ant_colony.cpp` / `.hpp`
    * Anthill colony algorithm

## 3 Research Conclusions

### 3.1 Adding random steiner points

There's two random methods:

1. `steiner_neighbor_random()`: The initial implementation of it added random steiner points. I observed that random points in triangles with a constrained edge didn't benefit from this method, and updated it to skip those.
2. `steiner_constraint_random()`: Adds steiner points on constrained edges

These produce a list of random steiner points that reduce the obtuse triangle count. Three basic implementations to pick points from these lists:

1. 1st
2. Random
3. All of them

To interpret the following conclusions, here's a table that maps random methods and pick methods to the shorthand names used:

| | steiner_neighbor_random() | steiner_constraint_random() |
| --- | --- | --- |
| None | ∅ | ∅ |
| 1st | neighborFirst | constFirst |
| Random | neighbor1 | const1 |
| All of them | neighbor | const |

**Effects of {neighbor, const} in simulated annealing:**

1. Measurement data: N/A

2. Comparisons:

	- {`neighbor`, ∅} vs {∅, `const`} :
	    + Similar results & execution time
	- {`neighbor`, ∅} vs {`neighbor`, `const`}:
	    + Similar results & execution time
	- {`neighbor1`, ∅} vs {`neighborFirst`, ∅}:
	    + ΠSimilar results & execution time
	- {`neighborFirst`, ∅} vs {`neighbor`, ∅}:
	    + Similar results & execution time
	- {`neighborFirst`, ∅} vs {∅, ∅}
	    + Similar results & execution time

3. Conclusions:

	- Doesn't help

**Effects of {neighbor, const} for anthill colony:**

1. Measurement data : N/A

2. Comparisons:

	- {`neighbor`, ∅} vs {∅, `const`}:
		+ {∅, `const`} runs slower and produces worse results
	- {`neighbor`, ∅} vs {`neighbor`, `const`}:
	    + Similar results
	    + {`neighbor`, `const`} runs slower
	- {`neighbor1`, ∅} vs {`neighborFirst`, ∅}:
	    + Similar results & runtime
	- {`neighborFirst`, ∅} vs {`neighbor`, ∅}
	    + Similar results
		+ {`neighbor`, ∅} runs faster
	- {`neighborFirst`, ∅} vs {∅, ∅}
		+ {`neighborFirst`, ∅} results are much better (reduces added steiner points to half)
		+ {`neighborFirst`, ∅} is orders of magnitude slower

3. Conclusions:

	- Any choice other than {∅,∅} improves the results. We should use a random method.
	- `const` is always a bad choice
	- {`neighbor`, ∅} is the best choice


### 3.2 Method Parametes

simulated annealing (L,α,β):

- L: Should be kept less than 100. Higher values don't improve the triangulation.
- {α,β}: As long as α > β, the code performs equally well. In the case where α < β, the code prefers adding steiner points to obtuse triangle reduction and makes bad choices.

ant (L,α,β,χ,ψ,λ,κ):

- L: Should be kept less than 200. Higher values don't improve the triangulation.
- {α,β}: As long as α > β, the code performs equally well. In the case where α < β, the code prefers adding steiner points to obtuse triangle reduction and makes bad choices.
- {χ,ψ,λ,κ}: They have no effect

