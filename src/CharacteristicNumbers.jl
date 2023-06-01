module CharacteristicNumbers

export characteristic_number_generator
export det_polarized
export chromatic_polynomial, relative_chromatic_polynomial, all_invariants
export characteristic_number, characteristic_number_generator

include("tensors.jl")
export generic_rankr, tensor_n, cm_tensors

using HomotopyContinuation
using LinearAlgebra
using IterTools
using DataStructures
using Combinatorics
using Memoize

function chromatic_polynomial(T)
  a, n, _ = size(T)
  characteristic_number = characteristic_number_generator(T)
  return [characteristic_number(alpha=[a - 1 - k; zeros(Int64, n - 3); k]) for k in 0:a-1]
end

function relative_chromatic_polynomial(T)
  a, n, _ = size(T)
  characteristic_number = characteristic_number_generator(T)
  return [characteristic_number(beta=[a - 1 - k; zeros(Int64, n - 3); k]) for k in 0:a-1]
end

function all_invariants(T, ordinary=true, relative=true, startsols=1,
  xtol=1e-14, compile=true, show_progress=true, max_loops_no_progress=3)
  @assert ordinary || relative
  a, n, _ = size(T)
  characteristic_number = characteristic_number_generator(T)
  res = []
  for s in with_replacement_combinations([ordinary ? (1:n-1) : [1]; 
      relative ? (n+2:n+ n-1) : []], a-1)
    alpha = zeros(Int64,n-1)
    beta = zeros(Int64,n-1)
    for i in s
      if i == 1
        if ordinary
          alpha[1] += 1
        else
          beta[1] += 1
        end
      elseif i <= n-1
        alpha[i] += 1
      else
        beta[i-n] += 1
      end
    end
    if show_progress
      if ordinary && relative
        print(alpha," ",beta)
      elseif ordinary
        print(alpha)
      else
        print(beta)
      end
    end
    num = characteristic_number(alpha=alpha, beta=beta, startsols=startsols,
      xtol=xtol, compile=compile, show_progress=false, 
      max_loops_no_progress=max_loops_no_progress)
    if show_progress
      println(" ",num)
    end
    if ordinary && relative
      push!(res,((alpha,beta), num))
    elseif ordinary
      push!(res,(alpha, num))
    else
      push!(res,(beta, num))
    end
  end
  return res
end

# computes the (relative) characteristic number of T in CC^a\ot CC^n\ot \CC^n, 
# possibly in combination. 
# here alpha and beta are vectors of length n-1, and their total sum is a-1
# alpha encodes the multidegree wanted of ordinary characteristic numbers
# and beta encodes the multidegree wanted of relative characteristic numbers
# The value of alpha[1] and beta[1] are indestinguishable, and the result is only 
# dependant on the sum
# 
# The settings here are somewhat conservative. It may be possible to reduce 
# max_loops_no_progress or increase startsols. To be more conservative, increase
# max_loops_no_progress, perhaps to 5 or so
function characteristic_number(T; alpha=nothing, beta=nothing, startsols=1,
  xtol=1e-14, compile=true, show_progress=true, max_loops_no_progress=3)
  return characteristic_number_generator(T)(alpha=alpha, beta=beta, startsols=startsols,
    xtol=xtol, compile=compile, show_progress=show_progress, max_loops_no_progress=max_loops_no_progress)
end

# returns a function computing (relative) characteristic numbers of T such that redundant work is avoided when calling this function multiple times with different arguments
function characteristic_number_generator(T)
  a, n, m = size(T)
  @assert n == m

  @var x[1:a]
  @var y[1:n, 1:n]
  @var v[1:a]

  M(e) = dropdims(sum(e .* T, dims=1), dims=1)
  m = M(x)

  z = vcat(x, vec(y))
  getz(x0) = vcat(x0, vec(inv(M(x0))))

  function spanning_set_indices(fs)
    F = qr(((z0, exp) -> exp(z => z0)).([getz(randn(a)) for _ in 1:length(fs)], [
        exp for _ in 1:1, exp in fs]), ColumnNorm())
    r = findfirst(abs.(diag(F.R)) .< 1e-10)
    if r === nothing
      r = length(fs)
    else
      r = r - 1
    end
    return sort(F.p[1:r])
  end

  @memoize function relative_f_from_derivative(ix)
    nderivatives = length(ix)
    k = n - nderivatives
    if k < n - k
      # return det_polarized( [fill(m,k); [T[i,:,:] for i in ix] ])
      p = det(m)
      for i in ix
        p = differentiate(p, x[i])
      end
      return p
    else
      return det_polarized([y * T[i, :, :] for i in ix])
    end
  end

  @memoize function get_derivative_coordinates_relative(nderivatives)
    if nderivatives == 1
      spanning_set = [[i] for i in 1:a]
    else
      spanning_set = Set{Vector{Int64}}()
      for s in get_derivative_coordinates_relative(nderivatives - 1)
        for i in 1:a
          push!(spanning_set, sort([s; i]))
        end
      end
      spanning_set = sort(collect(spanning_set))
    end
    spanning_fs = [relative_f_from_derivative(s) for s in spanning_set]
    ixs = spanning_set_indices(spanning_fs)
    return [spanning_set[i] for i in ixs]
  end

  @memoize function get_fs(k; relative=false)
    if k == 1
      # this is the same for relative and nonrelative
      return x
    end
    if relative
      nderivatives = n - k
      return [relative_f_from_derivative(s)
              for s in get_derivative_coordinates_relative(nderivatives)]
    else
      if k < n - k
        spanning_fs = vec([det(m[ix, jx]) for (ix, jx) in product(subsets(1:n, k), subsets(1:n, k))])
      else
        # using det(m[ix,jx]) = det(m)*det(inv(m)[1:n - jx, 1:n - ix])
        # so the below are the same functions as above up to multiplication by det(m)
        spanning_fs = vec([det(y[ix, jx]) for (ix, jx) in product(subsets(1:n, n - k), subsets(1:n, n - k))])
      end
      ixs = spanning_set_indices(spanning_fs)
      spanning_fs = [spanning_fs[i] for i in ixs]
      return spanning_fs
    end
  end

  function characteristic_number(; alpha=nothing, beta=nothing, startsols=1,
    xtol=1e-14, compile=true, show_progress=true, max_loops_no_progress=3)
    if alpha === nothing
      alpha = fill(0, n - 1)
    end
    if beta === nothing
      beta = fill(0, n - 1)
    end
    @assert length(alpha) == n - 1
    @assert length(beta) == n - 1
    @assert sum(alpha) + sum(beta) == a - 1

    x0s = randn(a, startsols)
    v0 = vec(svd(x0s).U[:, 1])
    # v0 can be any vector, pick one take make the next line not blow up x0's too much
    x0s ./= sum(v0 .* x0s, dims=1)
    x0s = [x0s[:, j] for j in 1:size(x0s, 2)]
    eqs = [sum(v .* x) .- 1]

    append!(eqs, vec(y * m - I))
    z0s = [getz(x0) for x0 in x0s]

    pv = copy(v)
    p0 = v0

    function add_constraint(k, dim; relative=false)
      if dim == 0
        return
      end
      fs = get_fs(k, relative=relative)
      # println(k," ",relative)
      # display(fs)
      @var p[relative ? 2 : 1, k, 1:length(fs), 1:dim]
      append!(eqs, [sum(eq * v for (eq, v) in zip(fs, p[:, i])) for i in 1:dim])

      p0cur = transpose(linear_forms_vanishing_on_prefix(
        [eq(z => z0) for eq in fs, z0 in z0s], dim))
      append!(pv, vec(p))
      append!(p0, vec(p0cur))
    end

    for (k, dim) in enumerate(alpha)
      add_constraint(k, dim, relative=false)
    end
    for (k, dim) in enumerate(beta)
      add_constraint(k, dim, relative=true)
    end

    F = System(eqs, variables=z, parameters=pv)
    sols = count(z0 -> norm(Float64.(evaluate(F.expressions, F.variables => z0,
        F.parameters => p0))) < 1e-12, z0s)
    z0s = z0s[1:sols]
    # if show_progress
    #   println("starting with ", sols, " solutions")
    # end
    return nsolutions(monodromy_solve(F, z0s, p0, compile=compile, show_progress=show_progress,
      max_loops_no_progress=max_loops_no_progress, unique_points_atol=xtol))
    println("done")
  end
  return characteristic_number
end


# let E_k be the k! times the kth elementary symmetric function of the eigenvalues of 
# an n times n matrix. This function evaluates the polarized form of E_k, that is, E_k 
# considered as a symmetric multilinear map. The input is thus given as a vector of 
# n times n matrices of length k

# relations
# 1/k! det_polarized(X,..(k times)..,X,A,...,Z) = det X * det_polarized(X^-1 A, ... X^-1 Z)
# for instance both are equal to the coefficient of a..z in the expansion of
# det (X + aA + ... + zZ). To see it for the right side, set A = ... = Z
# and write det(X + aA) = det X * det(I + (a+..+z) X^-1 A) and the coefficinet of a..z is the elementary function of the eigenvalues of A. Then note that the expression with A..Z possibily different is multilinear and symmetric, so it must agree with the left hand side
function det_polarized(Xs, deduplicate_terms=true)
  n = length(Xs)

  function conjugate_by_transposition(p, i, j)
    q = [k == j ? i : (k == i ? j : k) for k in p]
    (q[i], q[j]) = (q[j], q[i])
    return q
  end

  equal = DisjointSets(Combinatorics.permutations(1:n))
  if deduplicate_terms
    for (i, X) in enumerate(Xs)
      j = findnext(m -> m == X, Xs, i + 1)
      if j !== nothing
        for p in equal
          union!(equal, p, conjugate_by_transposition(p, i, j))
        end
      end
    end
  end

  sizes = Dict()
  for p in equal
    r = find_root!(equal, p)
    sizes[r] = get(sizes, r, 0) + 1
  end
  tot = 0
  for (p, sz) in sizes
    term = levicivita(p) * sz
    seen = Set()
    for i in 1:n
      if i in seen
        continue
      end
      M = Xs[i]
      push!(seen, i)
      j = p[i]
      while j != i
        M = M * Xs[j]
        push!(seen, j)
        j = p[j]
      end
      term *= tr(M)
    end
    tot += term
  end
  return tot
end
  
function linear_forms_vanishing_on_prefix(vs, k)
  return transpose(qr(hcat(vs, randn(size(vs)[1], k))).Q[:, end-k+1:end])
end

end # module CharacteristicNumbers
