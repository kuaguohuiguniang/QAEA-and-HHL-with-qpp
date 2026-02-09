## Project Overview

This project implements two example quantum algorithms using the **qpp** (Quantum++) library:

1. **HHL Algorithm**
   An implementation of the Harrow–Hassidim–Lloyd (HHL) algorithm for solving the **Quantum Linear Systems Problem (QLSP)**.

2. **QAEA Algorithm**
   An implementation of the Quantum Amplitude Estimation Algorithm (QAEA) applied to a quasi-one-dimensional Navier–Stokes equation, modeling fluid flow through a **convergent–divergent (de Laval) nozzle**.

---

## Build Instructions

### Compile the HHL module

```bash
make HHL
```

### Run the HHL executable

```bash
./hhl_exec
```

### Compile the nozzle module

```bash
make nozzle
```

### Run the nozzle executable

```bash
./nozzle_exec
```

---

## References

1. **HHL Algorithm**
   Harrow, Hassidim, and Lloyd.
   *Quantum algorithm for solving linear systems of equations.*
   [https://arxiv.org/abs/0811.3171](https://arxiv.org/abs/0811.3171)

2. **Quantum Fluid Simulation**
   *Finding flows of a Navier–Stokes fluid through quantum computing.*
   Nature Quantum Information.
   [https://www.nature.com/articles/s41534-020-00291-0](https://www.nature.com/articles/s41534-020-00291-0)
