#!/usr/bin/env python3
"""
mandelbrot_hypergroups.py

Computational verification for:
"Canonical Hypergroups from External Angles of the Mandelbrot Set:
 Computation and Classification by Period"
by Mohammad F. Marashdeh

This script reproduces every computational claim in the paper:
  1. Ray pair enumeration for periods 2, 3, 4, 5
  2. Hyperoperation (multiplication) tables
  3. Associativity verification (all triples)
  4. Transposition axiom (join space) verification
  5. Fundamental relation beta* computation
  6. Sub-hypergroup enumeration (simplicity check)

Usage:
    python3 mandelbrot_hypergroups.py

No external libraries required. Python 3.6+ needed (f-strings).
All arithmetic is exact (integer mod 2^n - 1).
"""

from fractions import Fraction
from itertools import combinations


# ============================================================
# Part 1: Ray Pair Enumeration
# ============================================================

def doubling_orbit(a, mod):
    """Compute the orbit of a under k -> 2k mod (mod)."""
    orbit = []
    x = a % mod
    for _ in range(mod + 1):
        if x in orbit:
            break
        orbit.append(x)
        x = (2 * x) % mod
    return orbit


def exact_period(a, mod):
    """Return the exact period of a under doubling mod (mod)."""
    return len(doubling_orbit(a, mod))


def chord_length(a, b, mod):
    """Arc length of the shorter chord between a/mod and b/mod on the circle."""
    fa, fb = Fraction(a, mod), Fraction(b, mod)
    lo, hi = min(fa, fb), max(fa, fb)
    short = hi - lo
    return min(short, 1 - short)


def check_non_crossing(a, b, mod, n):
    """Check if the orbit portrait of the pair {a, b} is non-crossing."""
    chords = []
    for k in range(n):
        c1 = Fraction((a * pow(2, k, mod)) % mod, mod)
        c2 = Fraction((b * pow(2, k, mod)) % mod, mod)
        chords.append((min(c1, c2), max(c1, c2)))
    for i in range(len(chords)):
        for j in range(i + 1, len(chords)):
            a1, b1 = chords[i]
            a2, b2 = chords[j]
            if a1 < a2 < b1 < b2:
                return False
            if a2 < a1 < b2 < b1:
                return False
    return True


def find_root_pairs(n):
    """
    Find all root ray pairs for period n using the shortest-chord criterion.

    A pair {a, b} of period-n angles is a root ray pair iff:
      (i)   its orbit portrait is non-crossing, AND
      (ii)  it is the shortest chord in its orbit portrait.

    This reproduces the known pairs for n = 2, 3, 4 (verified against
    Douady-Hubbard and Schleicher's enumerations) and computes the
    pairs for n = 5.

    Returns: list of (a, b) pairs, modulus
    """
    mod = 2**n - 1
    period_n = [k for k in range(1, mod) if exact_period(k, mod) == n]

    pairs = []
    for a, b in combinations(period_n, 2):
        if not check_non_crossing(a, b, mod, n):
            continue
        # Check that {a, b} is the shortest chord in its orbit portrait
        my_len = chord_length(a, b, mod)
        is_shortest = True
        for k in range(1, n):
            c1 = (a * pow(2, k, mod)) % mod
            c2 = (b * pow(2, k, mod)) % mod
            if chord_length(c1, c2, mod) < my_len:
                is_shortest = False
                break
        if is_shortest:
            pairs.append((a, b))

    return pairs, mod


# ============================================================
# Part 2: Hypergroup Construction
# ============================================================

class MandelbrotHypergroup:
    """
    Represents H_n = {e} ∪ {roots of period-n hyperbolic components}
    with the external-angle hyperoperation.

    For composite 2^n - 1, also supports the extended set
    {e} ∪ H_d ∪ H_n where d | n (for the period-4 closure).
    """

    def __init__(self, n, ray_pairs=None, extra_pairs=None):
        """
        Construct H_n.

        Args:
            n: period
            ray_pairs: list of (a, b) pairs; computed if None
            extra_pairs: additional ray pairs from lower periods
                         (e.g., the period-2 pair for the H_4 closure)
        """
        self.n = n
        self.mod = 2**n - 1

        if ray_pairs is None:
            ray_pairs, _ = find_root_pairs(n)
        self.ray_pairs = ray_pairs

        # Build equivalence map: element -> class representative
        self.equiv = {}
        for a in range(self.mod):
            self.equiv[a] = a
        all_pairs = list(ray_pairs)
        if extra_pairs:
            all_pairs += extra_pairs
        for a, b in all_pairs:
            r = min(a, b)
            self.equiv[a] = r
            self.equiv[b] = r

        # Build classes: representative -> list of members
        self.classes = {}
        for a in range(self.mod):
            r = self.equiv[a]
            if r not in self.classes:
                self.classes[r] = []
            self.classes[r].append(a)

        # Representatives (the elements of H_n)
        self.reps = sorted(self.classes.keys())

        # Involution: inv(p) = γ(-α) for α ∈ A(p)
        self.inv_map = {}
        for r in self.reps:
            neg_members = [(-a) % self.mod for a in self.classes[r]]
            self.inv_map[r] = self.equiv[neg_members[0]]

    def hyper(self, p, q):
        """Compute [p] ⊕ [q] = {eq(a+b) : a ∈ [p], b ∈ [q]}."""
        result = set()
        for a in self.classes[self.equiv[p]]:
            for b in self.classes[self.equiv[q]]:
                result.add(self.equiv[(a + b) % self.mod])
        return frozenset(result)

    def name(self, r):
        """Human-readable name for a class representative."""
        if r == 0:
            return "e"
        return f"[{r}]"

    def fmt(self, s):
        """Format a set of representatives."""
        elems = sorted(s, key=lambda x: self.reps.index(x))
        return "{" + ", ".join(self.name(x) for x in elems) + "}"


# ============================================================
# Part 3: Axiom Verification
# ============================================================

def check_associativity(H):
    """
    Check associativity for all |H|^3 ordered triples.
    Returns (is_associative, violation_count, first_violation).
    """
    violations = 0
    first = None
    for x in H.reps:
        for y in H.reps:
            for z in H.reps:
                xy = H.hyper(x, y)
                left = set()
                for w in xy:
                    left |= set(H.hyper(w, z))

                yz = H.hyper(y, z)
                right = set()
                for w in yz:
                    right |= set(H.hyper(x, w))

                if left != right:
                    violations += 1
                    if first is None:
                        first = (x, y, z, frozenset(left), frozenset(right))

    return violations == 0, violations, first


def check_transposition(H):
    """
    Check the transposition axiom: z ∈ x⊕y ⟺ x ∈ z⊕inv(y).
    Returns (holds, violation_count).
    """
    violations = 0
    for x in H.reps:
        for y in H.reps:
            for z in H.reps:
                z_in_xy = z in H.hyper(x, y)
                x_in_z_invy = x in H.hyper(z, H.inv_map[y])
                if z_in_xy != x_in_z_invy:
                    violations += 1
    return violations == 0, violations


def compute_beta_star(H):
    """
    Compute the fundamental relation β* (transitive closure of β).
    Returns the number of β*-classes.
    """
    # Find direct β-relations from products
    related = set()
    for i in H.reps:
        for j in H.reps:
            prod = list(H.hyper(i, j))
            for a in prod:
                for b in prod:
                    if a != b:
                        related.add((min(a, b), max(a, b)))

    # Transitive closure via union-find
    parent = {i: i for i in H.reps}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for a, b in related:
        union(a, b)

    classes = set(find(r) for r in H.reps)
    return len(classes)


def check_simplicity(H):
    """
    Check if H is simple (no proper non-trivial sub-hypergroups).
    Returns (is_simple, proper_sub) where proper_sub is None or a set.
    """
    non_id = [r for r in H.reps if r != 0]
    for size in range(1, len(non_id)):
        for combo in combinations(non_id, size):
            sub = {0} | set(combo)
            # Check involution closure
            if not all(H.inv_map[x] in sub for x in sub):
                continue
            # Check product closure
            closed = True
            for a in sub:
                for b in sub:
                    if not H.hyper(a, b).issubset(sub):
                        closed = False
                        break
                if not closed:
                    break
            if closed:
                return False, sub
    return True, None


# ============================================================
# Part 4: Multiplication Table Printing
# ============================================================

def print_multiplication_table(H, label=""):
    """Print the full multiplication table."""
    if label:
        print(f"\n{'='*60}")
        print(f"  Multiplication table for {label}")
        print(f"{'='*60}")

    # Header
    header = f"{'⊕':>8s} |"
    for y in H.reps:
        header += f" {H.name(y):>12s}"
    print(header)
    print("-" * len(header))

    for x in H.reps:
        row = f"{H.name(x):>8s} |"
        for y in H.reps:
            result = H.hyper(x, y)
            row += f" {H.fmt(result):>12s}"
        print(row)


# ============================================================
# Part 5: Main Verification
# ============================================================

def verify_period(n, label, ray_pairs=None, extra_pairs=None):
    """Run all verifications for a given period."""
    print(f"\n{'#'*70}")
    print(f"#  Period {n}: {label}")
    print(f"{'#'*70}")

    # Construct hypergroup
    H = MandelbrotHypergroup(n, ray_pairs=ray_pairs, extra_pairs=extra_pairs)

    # Ray pairs
    print(f"\nModulus: 2^{n} - 1 = {H.mod}")
    print(f"Number of classes: {len(H.reps)}")
    print(f"Ray pairs ({len(H.ray_pairs)}):")
    for a, b in sorted(H.ray_pairs):
        print(f"  {{{a}, {b}}}  =  {{{a}/{H.mod}, {b}/{H.mod}}}")
    if extra_pairs:
        print(f"Extra pairs (from lower periods):")
        for a, b in extra_pairs:
            print(f"  {{{a}, {b}}}")

    print(f"\nEquivalence classes:")
    for r in H.reps:
        print(f"  {H.name(r):>6s} = {H.classes[r]}")

    print(f"\nInvolution:")
    for r in H.reps:
        print(f"  inv({H.name(r)}) = {H.name(H.inv_map[r])}")

    # Multiplication table (skip for large H)
    if len(H.reps) <= 16:
        print_multiplication_table(H, label)

    # Closure check
    reps_set = set(H.reps)
    closure_ok = True
    for x in H.reps:
        for y in H.reps:
            result = H.hyper(x, y)
            if not result.issubset(reps_set):
                closure_ok = False
                escapes = result - reps_set
                print(f"\n  CLOSURE FAILURE: {H.name(x)} ⊕ {H.name(y)} "
                      f"produces {sorted(escapes)} outside H")
                break
        if not closure_ok:
            break
    print(f"\nClosure: {'PASS' if closure_ok else 'FAIL'}")

    # Associativity
    total_triples = len(H.reps) ** 3
    assoc_ok, assoc_violations, first_v = check_associativity(H)
    print(f"\nAssociativity ({total_triples} triples): "
          f"{'PASS' if assoc_ok else 'FAIL'}")
    if not assoc_ok:
        print(f"  Violations: {assoc_violations} / {total_triples}")
        if first_v:
            x, y, z, left, right = first_v
            print(f"\n  First counterexample: "
                  f"({H.name(x)}, {H.name(y)}, {H.name(z)})")
            xy = H.hyper(x, y)
            print(f"    {H.name(x)} ⊕ {H.name(y)} = {H.fmt(xy)}")
            print(f"    LHS = ({H.name(x)} ⊕ {H.name(y)}) ⊕ {H.name(z)}:")
            for w in sorted(xy, key=lambda r: H.reps.index(r)):
                print(f"      {H.name(w)} ⊕ {H.name(z)} = "
                      f"{H.fmt(H.hyper(w, z))}")
            print(f"    LHS = {H.fmt(left)}")
            yz = H.hyper(y, z)
            print(f"    {H.name(y)} ⊕ {H.name(z)} = {H.fmt(yz)}")
            print(f"    RHS = {H.name(x)} ⊕ ({H.name(y)} ⊕ {H.name(z)}):")
            for w in sorted(yz, key=lambda r: H.reps.index(r)):
                print(f"      {H.name(x)} ⊕ {H.name(w)} = "
                      f"{H.fmt(H.hyper(x, w))}")
            print(f"    RHS = {H.fmt(right)}")

    # Transposition axiom (join space)
    if closure_ok:
        trans_ok, trans_violations = check_transposition(H)
        print(f"\nTransposition axiom ({total_triples} triples): "
              f"{'PASS' if trans_ok else 'FAIL'}")
        if not trans_ok:
            print(f"  Violations: {trans_violations}")

    # Beta*
    if closure_ok:
        beta_classes = compute_beta_star(H)
        print(f"\nβ* classes: {beta_classes} "
              f"({'1-hypergroup' if beta_classes == 1 else 'NOT 1-hypergroup'})")

    # Simplicity
    if closure_ok and assoc_ok:
        is_simple, proper_sub = check_simplicity(H)
        print(f"\nSimplicity: {'SIMPLE' if is_simple else 'NOT SIMPLE'}")
        if not is_simple:
            print(f"  Proper sub-hypergroup: {sorted(proper_sub)}")

    # Summary
    is_hypergroup = closure_ok and assoc_ok
    print(f"\n  => {label} is "
          f"{'a canonical hypergroup' if is_hypergroup else 'NOT a canonical hypergroup'}.")

    return H


def main():
    print("=" * 70)
    print("  Mandelbrot Hypergroup Verification")
    print("  Companion code for the paper by M. F. Marashdeh")
    print("=" * 70)

    # ---- Period 2 ----
    verify_period(2, "H_2")

    # ---- Period 3 ----
    verify_period(3, "H_3")

    # ---- Period 4 (H_4 alone) ----
    # First show H_4 is not closed
    pairs4, mod4 = find_root_pairs(4)
    H4 = MandelbrotHypergroup(4, ray_pairs=pairs4)
    print(f"\n{'#'*70}")
    print(f"#  Period 4: H_4 (without closure)")
    print(f"{'#'*70}")
    print(f"\nRay pairs: {pairs4}")
    reps4 = set(H4.reps)
    for x in H4.reps:
        for y in H4.reps:
            result = H4.hyper(x, y)
            escapes = result - reps4
            if escapes:
                print(f"  H_4 NOT CLOSED: {H4.name(x)} ⊕ {H4.name(y)} "
                      f"= {H4.fmt(result)}")
                print(f"    Elements {sorted(escapes)} are outside H_4")
                break
        else:
            continue
        break

    # ---- Period 4 closure: {e} ∪ H_2 ∪ H_4 ----
    # The period-2 pair in Z/15Z is {5, 10}
    period2_pair_in_15 = [(5, 10)]
    verify_period(4, "{e} ∪ H_2 ∪ H_4 (closure)",
                  ray_pairs=pairs4, extra_pairs=period2_pair_in_15)

    # ---- Period 5 ----
    verify_period(5, "H_5")

    # ---- Period 7 ----
    verify_period(7, "H_7")

    # ---- Summary ----
    print(f"\n\n{'='*70}")
    print("  SUMMARY")
    print(f"{'='*70}")
    print("  n=1: H_1 = {e}, trivially a canonical hypergroup.")
    print("  n=2: H_2 is a 2-element canonical hypergroup (Krein). PASS.")
    print("  n=3: H_3 is a 4-element canonical hypergroup. PASS.")
    print("  n=4: H_4 is NOT closed (sums escape to period-2 roots).")
    print("        Closure {e}∪H_2∪H_4 fails associativity. FAIL.")
    print("  n=5: H_5 is closed (31 prime) but fails associativity. FAIL.")
    print("  n=7: H_7 is closed (127 prime) but fails associativity. FAIL.")
    print()


if __name__ == "__main__":
    main()
