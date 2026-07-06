# Prime-Leaf Tree Theory: Recap So Far (v0.1 — original informal draft)

> This is the original informal recap the paper grew from, kept for provenance.
> The formalized and extended version is `../prime-leaf-tree-theory.tex`.

## 1. Starting intuition

The central idea is that ordinary multiplication hides structure.

In standard arithmetic, we write 3 × 2 = 6 and treat the result as the same kind of
object as 1+1+1+1+1+1 = 6. But in the theory we are developing, these are **not the
same kind of construction**.

The key principle is:

> You never get "six liters" from multiplication.
> You get a structured object whose ordinary numerical shadow is 6.

So 6L obtained by addition is a flat additive quantity, while 3 × 2 represents a
structured grouping: three groups of two, or two groups of three, depending on the
order. Therefore, multiplication should not immediately collapse back into ordinary
numbers.

## 2. Addition and multiplication live in different worlds

**Additive world.** The additive world contains flat quantities. Addition behaves like
accumulation: + : 𝒜 × 𝒜 → 𝒜. For example L+L+L+L+L+L = 6L. This is a single-unit
quantity.

**Multiplicative world.** The multiplicative world contains structured objects.
Multiplication does not return an ordinary number. It returns a **multiplication
tree**: ⊠ : 𝒯 × 𝒯 → 𝒯. So instead of saying 2 × 3 = 6, we say L₂ ⊠ L₃ is a structured
multiplicative object whose shadow is 6.

## 3. Prime leaves

We use the term **prime leaves** instead of atoms. Prime leaves are the primitive
building blocks of multiplication trees: L₂, L₃, L₅, L₇, … Each L_p is a prime leaf
associated with the prime shadow p. In the first version of the theory, we allow
ourselves to label leaves by known primes. Later, a deeper version may start from
abstract leaves and derive their numerical shadows.

## 4. Multiplication trees

A multiplication tree is defined recursively: T ::= L_p or T ::= T₁ ⊠ T₂.
Examples: L₂; L₂ ⊠ L₃; L₃ ⊠ L₂; (L₂ ⊠ L₃) ⊠ L₅; L₂ ⊠ (L₃ ⊠ L₅).
At the fundamental level, these are structured objects, not ordinary numbers.

## 5. Non-commutativity

In ordinary arithmetic 2 × 3 = 3 × 2. In prime-leaf tree theory this is not assumed:
L₂ ⊠ L₃ ≠ L₃ ⊠ L₂. The first means "a structure built from L₂, then L₃"; the second
"a structure built from L₃, then L₂". They may have the same numerical shadow, but they
are not the same tree.

## 6. Non-associativity

In ordinary arithmetic (2×3)×5 = 2×(3×5). In the tree theory this is also not assumed:
(L₂ ⊠ L₃) ⊠ L₅ ≠ L₂ ⊠ (L₃ ⊠ L₅). The parenthesization stores construction history.
Multiplication is not just about the final value; it records the way the object was
built.

## 7. The shadow map / norm

We introduce a projection map N : 𝒯 → ℕ which forgets the tree structure and returns
the ordinary numerical shadow. For prime leaves N(L_p) = p; for composed trees
N(T₁ ⊠ T₂) = N(T₁)N(T₂). So N(L₂ ⊠ L₃) = 6 and N(L₃ ⊠ L₂) = 6, but
L₂ ⊠ L₃ ≠ L₃ ⊠ L₂. Many different multiplicative structures can have the same ordinary
numerical shadow.

## 8. Meaning of composite numbers

Composite numbers are not primitive multiplicative objects. There is no fundamental
multiplicative object L₆. Instead, 6 appears as the shadow of trees such as L₂ ⊠ L₃ or
L₃ ⊠ L₂. Likewise 4 appears as N(L₂ ⊠ L₂), and 8 can appear as either
N((L₂ ⊠ L₂) ⊠ L₂) or N(L₂ ⊠ (L₂ ⊠ L₂)). The ordinary number 8 is the same shadow, but
the trees are distinct.

## 9. Distributivity is not fundamental

In ordinary arithmetic a(b+c) = ab+ac. But in this theory multiplication and addition
live in different structural worlds, so we do **not** impose
a ⊠ (b+c) = (a ⊠ b) + (a ⊠ c). The left side means "multiply by an already-added
object"; the right side means "multiply separately, then add the resulting structures".
Those are different constructions. At most, distributivity may appear after applying a
projection or forgetting map: N(a ⊠ (b+c)) = N((a ⊠ b) + (a ⊠ c)) belongs to the
collapsed ordinary theory, not the fundamental theory.

> Distributivity is emergent, not axiomatic.

## 10. Ordinary arithmetic as a quotient theory

Ordinary arithmetic is recovered by imposing equivalences on the richer tree theory.
First quotient by associativity: (T₁ ⊠ T₂) ⊠ T₃ ∼ T₁ ⊠ (T₂ ⊠ T₃). Then by
commutativity: T₁ ⊠ T₂ ∼ T₂ ⊠ T₁. After these identifications, a tree becomes just a
prime multiset L₂^α₂ L₃^α₃ L₅^α₅ ⋯, and applying the norm gives the ordinary integer
∏_p p^α_p. So the hierarchy is:

multiplication trees → prime multisets → ordinary integers.

Classical arithmetic is the coarse-grained shadow of the prime-leaf tree theory.

## 11. Degeneracy of an integer

Let n = ∏_p p^α_p and Ω(n) = Σ_p α_p (number of prime leaves with multiplicity). The
number of ordered, parenthesized multiplication trees with shadow n is

g_T(n) = C_{Ω(n)−1} · Ω(n)! / ∏_p α_p!

where C_k = (1/(k+1)) · binom(2k, k) is the k-th Catalan number.

Example: 6 = 2·3, Ω(6) = 2, g_T(6) = C₁ · 2!/1!1! = 2 (the trees L₂ ⊠ L₃ and L₃ ⊠ L₂).
For 12 = 2²·3, Ω(12) = 3, g_T(12) = C₂ · 3!/2!1! = 2·3 = 6: twelve has six distinct
multiplication-tree microstates.

## 12. Riemann zeta function in the theory

The ordinary ζ(s) = Σ 1/n^s is not the most fundamental zeta function here. The more
fundamental object is the **tree zeta function**

Z_T(s) = Σ_{T ∈ 𝒯} 1/N(T)^s = Σ_{n≥2} g_T(n)/n^s.

Z_T counts multiplicative structures; ordinary ζ counts projected integers. Ordinary
Riemann zeta appears after collapsing all trees with the same shadow.

## 13. Euler product

At the quotient level every multiplicative object is a prime multiset, each prime leaf
appears independently, and the geometric series gives

ζ(s) = ∏_p 1/(1 − p^{−s}).

The Euler product is the partition function of the quotient theory where only
prime-leaf multiplicities matter.

## 14. Prime zeta function connection

Define P(s) = Σ_p p^{−s}. A multiplication tree is either a prime leaf or a binary
composition of two multiplication trees. Therefore Z_T(s) = P(s) + Z_T(s)², so

Z_T(s) = (1 − √(1 − 4P(s)))/2.

The full tree zeta function is governed by the prime zeta function:
prime leaves → tree zeta → quotient → Riemann zeta.

## 15. Euler totient function

For a tree T we can define a pulled-back totient
φ_T(T) = N(T) ∏_{p ∈ supp(T)} (1 − 1/p).
For T = L₂ ⊠ L₃: φ_T(T) = 6 · (1/2) · (2/3) = 2, recovering the ordinary Euler totient
after projection. The tree theory could also define a richer, tree-sensitive totient
that distinguishes L₂ ⊠ L₃ from L₃ ⊠ L₂; classical φ cannot see that distinction.

## 16. Physical interpretation

The multiplication tree is like a microscopic state; the ordinary integer is an
observable or coarse-grained shadow. L₂ ⊠ L₃ and L₃ ⊠ L₂ are different microstates with
the same measured value 6. With E(T) = log N(T) we get 1/N(T)^s = e^{−sE(T)}, so
Z_T(s) = Σ_T e^{−sE(T)} looks like a partition function over multiplicative
microstates, with s an inverse-temperature-like parameter.

## 17. Core philosophical claim

- Prime leaves are fundamental multiplicative primitives.
- Composite numbers are not fundamental objects.
- Multiplication creates structured trees.
- Ordinary integers are shadows of these trees.
- Addition produces flat quantities.
- Distributivity, commutativity, and associativity are not fundamental; they appear
  only after quotienting or projecting the richer theory.

> Prime-leaf tree theory treats ordinary arithmetic as the coarse-grained shadow of a
> richer non-distributive, non-commutative, non-associative multiplicative structure
> built from prime leaves.

## 18. Open foundational questions

1. Are prime leaves labeled by known primes, or can abstract leaves derive their
   shadows later?
2. Should left-right order be fundamental? (So far, yes.)
3. Should parenthesization be fundamental? (So far, yes.)
4. How should addition interact with trees? A complete theory needs a precise bridge
   between 𝒜 and 𝒯.
5. What is the right generalized zeta function?
6. Can Riemann-type questions be lifted to Z_T?

## 19. Current working name

**Prime-Leaf Tree Theory** (alternatives: Multiplicative Tree Theory, Prime-Leaf
Arithmetic, Structured Multiplication Theory, Non-Distributive Prime-Leaf Theory, Tree
Arithmetic of Prime Leaves).
