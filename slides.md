# PyOP2: A High-Level Framework for Performance-Portable Simulations on Unstructured Meshes

!SLIDE left title

# PyOP2: A High-Level Framework for Performance-Portable Simulations on Unstructured Meshes

## Florian Rathgeber<sup>1</sup>, Graham Markall<sup>1</sup>, Lawrence Mitchell<sup>3</sup>, Nicolas Loriant<sup>1</sup>, David Ham<sup>1,2</sup>, Carlo Bertolli<sup>1</sup>, Paul Kelly<sup>1</sup>

### <sup>1</sup> Department of Computing
### <sup>2</sup> Grantham Institute for Climate Change
### Imperial College London
### <sup>3</sup> EPCC, University of Edinburgh

!SLIDE left

# Computational science is hard

## Radical paradigm shifts in CSE due to many-core computing
* Scientists have to rewrite and hand-tune for each emerging platform
* Barrier: inhibits efficient utilisation of state-of-the-art computational resources

## Generative (meta)programming offers a solution
* Generate platform-specific implementations instead of hand-coding them
* Tailor to characteristics of the hardware and problem

!NOTES

# FEM is a versatile tool for science and engineering

## Tsunami simulation of the Hokkaido-Nansei-Oki tsunami of 1993

<iframe width="640" height="360" src="http://www.youtube.com/embed/Y6mM_PCNhq0?rel=0" frameborder="0" allowfullscreen></iframe>

The simulation was carried out with the [Fluidity multi-phase CFD code](http://amcg.ese.ic.ac.uk/index.php?title=FLUIDITY) solving the non-hydrostatic Navier-Stokes equations, using a [free surface and wetting and drying algorithm](http://amcg.ese.ic.ac.uk/index.php?title=Wetting_and_Drying) (courtesy [Simon Funke](http://www3.imperial.ac.uk/people/s.funke09/)).

!SLIDE huge

# The challenge

> How do we get performance portability for the finite element method without sacrificing generality?

!SLIDE left

# The strategy

## Get the abstractions right
... to isolate numerical methods from their mapping to hardware

## Start at the top, work your way down
... and make decisions at the highest abstraction level possible

## Harness the power of DSLs
... for generative, instead of transformative optimisations

!SLIDE left

# The tools

## Embedded domain-specific languages

... capture and *efficiently express characteristics* of the application/problem domain

## Runtime code generation

... encapsulates *specialist expertise* to deliver *problem- and platform-specific optimisations*

## In combination, they

* raise the level of abstraction and incorporate domain-specific knowledge
* decouple problem domains from their efficient implementation on different hardware
* capture design spaces and open optimisation spaces
* enable reuse of code generation and optimisation expertise and tool chains

!SLIDE huge

# The big picture

!SLIDE

}}} images/mapdes_abstraction_layers_overview.svg

!SLIDE huge

# Higher level abstraction

## From the equation to the finite element implementation

!SLIDE left

# FFC<sup>1</sup> takes equations in UFL

## Helmholtz equation
@@@ python
f = state.scalar_fields["Tracer"]

v = TestFunction(f)
u = TrialFunction(f)

lmbda = 1
a = (dot(grad(v), grad(u)) - lmbda * v * u) * dx

L = v*f*dx

solve(a == L, f)
@@@

<sub><sup>1</sup> [FFC](https://launchpad.net/ffc) is the FEniCS Form Compiler developed by the [FEniCS project](http://fenicsproject.org/)</sub>

!NOTES

## Fluidity extensions

* **`state.scalar_fields`** interfaces to Fluidity: read/write field of given name
* **`solve`** records equation to be solved and returns `Coefficient` for solution field

!SLIDE left

# ... and generates local assembly kernels

## Helmholtz OP2 kernel
@@@ clike
void kernel(double A[1][1], double *x[2],
            int j, int k) {
  // Kij - Jacobian determinant
  // FE0 - Shape functions
  // Dij - Shape function derivatives
  // W3  - Quadrature weights
  for (unsigned int ip = 0; ip < 3; ip++) {
    A[0][0] += (FE0[ip][j] * FE0[ip][k] * (-1.0)
      + (((K00 * D10[ip][j] + K10 * D01[ip][j]))
        *((K00 * D10[ip][k] + K10 * D01[ip][k]))
      + ((K01 * D10[ip][j] + K11 * D01[ip][j]))
        *((K01 * D10[ip][k] + K11 * D01[ip][k]))
      )) * W3[ip] * det;
  }
}
@@@

!SLIDE huge

# Lower level abstraction

## From the finite element implementation to its efficient parallel execution

!SLIDE left

# PyOP2 – a high-level framework for unstructured mesh computations

## Abstractions for unstructured meshes

* **Sets** of entities (e.g. nodes, edges, faces)
* **Mappings** between sets (e.g. from edges to nodes)
* **Datasets** holding data on a set (e.g. fields in finite-element terms)

## Mesh computations as parallel loops

* execute a *kernel* for all members of an iteration set in arbitrary order
* datasets accessed through at most one level of indirection via a mapping
* *access descriptors* specify which data is passed to the kernel and how it is addressed

## Multiple hardware backends via *runtime code generation*

* partioning/colouring for efficient scheduling and execution on different hardware
* currently supports CUDA/OpenCL, MPI in development, AVX support planned

!SLIDE left

# PyOP2 for finite-element computations

## Finite element local assembly
... means computing the *same kernel* for *every* mesh entity (cell, facet): a perfect match for the PyOP2 abstraction

## PyOP2 abstracts away data marshaling and parallel execution
* controls whether/how/when a matrix is assembled
* local assembly kernel is *translated* for and *efficiently executed* on the target architecture

## Global asssembly and linear algebra operations
... implemented as a thin wrapper on top of backend-specific linear algebra packages:  
*PETSc4py* on the CPU, *Cusp* on the GPU

!SLIDE left

# Finite element assembly and solve in PyOP2

@@@ python
from pyop2 import op2, ffc_interface

def solve(A, x, b):
    # Generate kernels for matrix and rhs assembly
    mass = ffc_interface.compile_form(A, "mass")[0]
    rhs  = ffc_interface.compile_form(b, "rhs")[0]

    # Extract coordinates (coords) and forcing function (f_vals)

    # Construct OP2 matrix to assemble into
    sparsity = op2.Sparsity((elem_node, elem_node), sparsity_dim) 
    mat = op2.Mat(sparsity, numpy.float64)
    f = op2.Dat(nodes, 1, f_vals, numpy.float64)

    # Assemble and solve
    op2.par_loop(mass, elements(3,3),
             mat((elem_node[op2.i[0]], elem_node[op2.i[1]]), op2.INC),
             coords(elem_node, op2.READ))

    op2.par_loop(rhs, elements(3),
             b(elem_node[op2.i[0]], op2.INC),
             coords(elem_node, op2.READ),
             f(elem_node, op2.READ))

    op2.solve(mat, b, x)
@@@

!SLIDE left

# UFL equations in Fluidity

## For each UFL equation in each time step:

![Fluidity-UFL-PyOP2-toolchain](images/fluidity_pyop2_pipeline.svg)

Instead of calling Fluidity's built-in advection-diffusion solver:

* Shell out to Python, execute the user's UFL equation
* FFC generates C++ code for local assembly of FE forms
* Instant JIT-compiles kernels and the parallel loops invoking them

!NOTES
* FFC + Instant invocations are cached

!SLIDE huge

# Preliminary performance results

!SLIDE left

# Experimental setup

##Solver
CG with Jacobi preconditioning using PETSc 3.3 (PyOP2), 3.2 (DOLFIN)

##CPU
2x Intel Xeon E5650 Westmere 6-core (HT off), 48GB RAM

##GPU
NVIDIA GeForce GTX680 (Kepler)

##Mesh
2D unit square meshed with triangles (200 - 204800 elements)

##Dolfin
Revision 6906, Tensor representation, CPP optimisations on, form compiler optimisations off

!SLIDE huge

# Single core

!SLIDE

}}} images/seq_runtime_linear.svg

!SLIDE

}}} images/seq_speedup_linear.svg

!SLIDE huge

# Parallel

!SLIDE

}}} images/cuda_runtime_linear.svg

!SLIDE

}}} images/cuda_speedup_linear.svg

!SLIDE left

# Resources

All the code mentioned is open source and available on *GitHub* and *Launchpad*. Try it!

## PyOP2
<https://github.com/OP2/PyOP2>

## FFC
<https://code.launchpad.net/~mapdes/ffc/pyop2>

## Fluidity
<https://code.launchpad.net/~fluidity-core/fluidity/pyop2>

## Benchmarks
<https://github.com/OP2/PyOP2_benchmarks>

## This talk
<http://kynan.github.com/wolfhpc2012>
