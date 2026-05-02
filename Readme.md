# Mandelbrot Hypergroups

Companion code for the paper:

**"Canonical Hypergroups from External Angles of the Mandelbrot Set: Computation and Classification by Period"**  
by Mohammad F. Marashdeh (Mutah University)

## Overview

This repository contains the Python code that verifies every computational claim in the paper. The script:

1. **Enumerates root ray pairs** for periods $n = 2, 3, 4, 5$ using the shortest-chord criterion on orbit portraits
2. **Computes hyperoperation tables** via Minkowski sums in $\mathbb{Z}/(2^n-1)\mathbb{Z}$
3. **Checks all five hypergroup axioms** (commutativity, associativity, identity, involution, reversibility)
4. **Verifies the transposition axiom** (join space property) by exhaustive enumeration
5. **Computes the fundamental relation** $\beta^*$ (1-hypergroup property)
6. **Checks simplicity** (no proper non-trivial sub-hypergroups)

## Requirements

- Python 3.6 or later
- No external libraries (uses only `fractions` and `itertools` from the standard library)

## Usage

```bash
python3 mandelbrot_hypergroups.py
```

The script runs all verifications and prints results to stdout. Runtime is under 1 second on a standard laptop.

## Key Results Reproduced

| Period | \|H_n\| | Closed | Associative | Violations | Canonical Hypergroup |
|--------|---------|--------|-------------|------------|---------------------|
| 2      | 2       | Yes    | Yes         | 0/8        | **Yes** (Krein)     |
| 3      | 4       | Yes    | Yes         | 0/64       | **Yes**             |
| 4      | 8*      | Yes*   | No          | 230/512    | No                  |
| 5      | 16      | Yes    | No          | 1636/4096  | No                  |
| 7      | 64      | Yes    | No          | 197698/262144 | No               |

*Period 4: H_4 itself is not closed; the table shows the 8-element closure {e} ∪ H_2 ∪ H_4.

## Algorithm Details

### Ray Pair Enumeration

For each pair {a, b} of period-n angles in Z/(2^n - 1)Z:

1. Compute the orbit portrait: the n chords {2^k·a, 2^k·b} for k = 0, ..., n-1
2. Check the **non-crossing condition**: no two chords interleave on the circle
3. Check the **shortest-chord condition**: {a, b} has the smallest arc length among all chords in its orbit portrait

This criterion reproduces the known ray pairs for n = 2, 3, 4 (verified against Douady-Hubbard and Schleicher's enumerations).

### Associativity Check

For each ordered triple (x, y, z) ∈ H_n³:
- Compute LHS = (x ⊕ y) ⊕ z = ∪_{w ∈ x⊕y} (w ⊕ z)
- Compute RHS = x ⊕ (y ⊕ z) = ∪_{w ∈ y⊕z} (x ⊕ w)
- Compare the two sets

## Citation

```bibtex
@article{marashdeh2025mandelbrot,
  title={Canonical Hypergroups from External Angles of the {M}andelbrot Set: Computation and Classification by Period},
  author={Marashdeh, Mohammad F.},
  journal={Experimental Mathematics},
  year={2025},
  note={Submitted}
}
```

## License

MIT License
