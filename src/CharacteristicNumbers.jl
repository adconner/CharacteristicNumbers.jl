module CharacteristicNumbers
export characteristic_number, chromatic_polynomial_coefficient, characteristic_numbers, chromatic_polynomial
export relative_characteristic_number, relative_chromatic_polynomial_coefficient, relative_characteristic_numbers, relative_chromatic_polynomial
export mat_tensors

using HomotopyContinuation
using LinearAlgebra
using IterTools
using DataStructures
using Combinatorics
using Memoize

using Debugger

export det_polarized
# when n is smaller than the size of the Xs, this is the polarization of the unique polynomial computing the nth elementary symmetric function of the eigen values of its input

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
  
function characteristic_number(T; alpha = nothing, beta = nothing, startsols=1, 
    xtol=1e-14, compile=true, show_progress=true, max_loops_no_progress=2, monodromy=true)
  return characteristic_number_generator(T)(alpha = alpha, beta = beta, startsols=startsols, 
    xtol=xtol, compile=compile, show_progress=show_progress, max_loops_no_progress=max_loops_no_progress, monodromy=monodromy)
end
  
export characteristic_number_generator
function characteristic_number_generator(T)
  a, n, m = size(T)
  @assert n == m
  
  @var x[1:a]
  @var y[1:n,1:n]
  @var v[1:a]
  
  M(e) = dropdims(sum(e .* T, dims=1), dims=1)
  m = M(x)

  z = vcat(x, vec(y))
  getz(x0) = vcat(x0, vec(inv(M(x0))))
  
  function spanning_set_indices(fs)
    F = qr( ((z0,exp) -> exp(z=>z0)).([getz(randn(a)) for _ in 1:length(fs)],[ 
        exp for _ in 1:1, exp in fs ]), ColumnNorm())
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
    k = n-nderivatives
    if k < n-k 
      # return det_polarized( [fill(m,k); [T[i,:,:] for i in ix] ])
      p = det(m)
      for i in ix
        p = differentiate(p,x[i])
      end
      return p
    else
      return det_polarized( [y*T[i,:,:] for i in ix] )
    end
  end

  @memoize function get_derivative_coordinates_relative(nderivatives)
    if nderivatives == 1
      spanning_set = [[i] for i in 1:a]
    else
      spanning_set = Set{Vector{Int64}}()
      for s in get_derivative_coordinates_relative(nderivatives-1)
        for i in 1:a
          push!(spanning_set,sort([s; i]))
        end
      end
      spanning_set = sort(collect(spanning_set))
    end
    spanning_fs = [ relative_f_from_derivative(s) for s in spanning_set ]
    ixs = spanning_set_indices(spanning_fs)
    return [spanning_set[i] for i in ixs]
  end
    
  @memoize function get_fs(k; relative=false)
    if k == 1 
      # this is the same for relative and nonrelative
      return x
    end
    if relative
      nderivatives = n-k
      return [relative_f_from_derivative(s)
        for s in get_derivative_coordinates_relative(nderivatives)]
    else
      if k < n-k
        spanning_fs = vec([ det(m[ix,jx]) for (ix,jx) in product(subsets(1:n,k),subsets(1:n,k)) ])
      else
        # using det(m[ix,jx]) = det(m)*det(inv(m)[1:n - jx, 1:n - ix])
        # so the below are the same functions as above up to multiplication by det(m)
        spanning_fs = vec([ det(y[ix,jx]) for (ix,jx) in product(subsets(1:n,n-k),subsets(1:n,n-k)) ])
      end
      ixs = spanning_set_indices(spanning_fs)
      spanning_fs = [spanning_fs[i] for i in ixs]
      return spanning_fs
    end
  end
  
  function characteristic_number(; alpha = nothing, beta = nothing, startsols=1, 
    xtol=1e-14, compile=true, show_progress=true, max_loops_no_progress=3)
    if alpha === nothing
      alpha = fill(0,n-1)
    end
    if beta === nothing
      beta = fill(0,n-1)
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
      fs = get_fs(k,relative=relative) 
      # println(k," ",relative)
      # display(fs)
      @var p[relative ? 2 : 1, k, 1:length(fs), 1:dim]
      append!(eqs,[ sum(eq*v for (eq,v) in zip(fs,p[:,i])) for i in 1:dim ])

      p0cur = transpose(linear_forms_vanishing_on_prefix(
        [eq(z=>z0) for eq in fs, z0 in z0s], dim))
      append!(pv,vec(p))
      append!(p0,vec(p0cur))
    end

    for (k,dim) in enumerate(alpha)
      add_constraint(k,dim,relative=false)
    end
    for (k,dim) in enumerate(beta)
      add_constraint(k,dim,relative=true)
    end
      
    F = System(eqs, variables=z, parameters=pv)
    sols = count(z0 -> norm(Float64.(evaluate(F.expressions, F.variables => z0,
        F.parameters => p0))) < 1e-12, z0s)
    z0s = z0s[1:sols]
    if show_progress
      println("starting with ", sols, " solutions")
    end
    return nsolutions(monodromy_solve(F, z0s, p0, compile=compile, show_progress=show_progress,
      max_loops_no_progress=max_loops_no_progress, unique_points_atol=xtol))
    println("done")
  end
  return characteristic_number
end
  
export characteristic_number_ss
function characteristic_number_ss(T; alpha = nothing, beta = nothing)
  a, n, m = size(T)
  @assert n == m
  if alpha === nothing
    alpha = fill(0,n-1)
  end
  if beta === nothing
    beta = fill(0,n-1)
  end
  @assert length(alpha) == n - 1
  @assert length(beta) == n - 1
  @assert sum(alpha) + sum(beta) == a - 1
  
  @var x[1:a]
  @var y[1:n,1:n]
  
  # v0 can be any vector, pick one take make the next line not blow up x0's too much
  eqs = [sum(randn(a) .* x) .- 1]
  
  M(e) = dropdims(sum(e .* T, dims=1), dims=1)
  m = M(x)

  append!(eqs, vec(y * m - I))
  
  z = vcat(x, vec(y))
  getz(x0) = vcat(x0, vec(inv(M(x0))))
  
  function add_constraint(k, dim; relative=false)
    if dim == 0
      return
    end
    k = n-k
    if relative
      @var t[1:a]
      mt = y*M(t)
    else
      @var t[1:n*n]
      mt = y*reshape(t,n,n)
    end
    D = det_polarized(fill(mt,k))
    
    eqdim_upperbound = min( binomial(n,k)^2, binomial(length(t)+k-1,k) )
    eqdim = rank( ((z0,t0) -> D(z=>z0, t=>t0)).([getz(randn(a)) for _ in 1:eqdim_upperbound],[ 
          randn(length(t)) for _ in 1:1, _ in 1:eqdim_upperbound ]) )
      
    Ds = [ D(t => randn(length(t))) for _ in 1:eqdim ]
    
    append!(eqs,[ sum(randn()*e for e in Ds) for _ in 1:dim])
  end

  for (k,dim) in enumerate(alpha)
    add_constraint(k,dim,relative=false)
  end
  for (k,dim) in enumerate(beta)
    add_constraint(k,dim,relative=true)
  end
    
  F = System(eqs, variables=z)
  return F
end

# We use the fact that Lambda^k V \cong Lambda^(n-k) V^* \ot \Lambda^n V as GL(V) modules, so linear combinations of k x k minors of an n x n matrix M is a linear combination of n-k x n-k minors of M^{-1} (up to a scale of det M). Specifically, a k x k minor is up to \pm det M the complementary n-k x n-k minor of (M^(-1))^t

export characteristic_number_old
# size(T) = (a,b,b), len(alpha) = b-1, sum(alpha) + sum(beta) = a-1
function characteristic_number_old(T, alpha; startsols=1, xtol=1e-14, compile=true,
  show_progress=true, max_loops_no_progress=2)
  a, b, c = size(T)
  @assert(length(alpha) == b - 1 && sum(alpha) == a - 1 && b == c)
  @var x[1:a]
  x0s = randn(a,startsols)
  v0 = vec(svd(x0s).U[:,1]) 
  # v0 can be any vector, pick one take make the next line not blow up x0's too much
  x0s ./= sum(v0 .* x0s,dims=1)
  x0s = [x0s[:,j] for j in 1:size(x0s,2)]

  eqs = [sum(v0 .* x) .- 1]

  M(e) = dropdims(sum(e.*T,dims=1),dims=1)
  m = M(x)

  @var y[1:b,1:b]
  append!(eqs,vec(y * m - I))

  z = vcat(x,vec(y))
  z0s = [vcat(x0,vec(inv(M(x0)))) for x0 in x0s]

  p = Vector{Variable}()
  p0 = Vector{Float64}()

  function minors_up_to(m,k)
    ms_dict = Dict{Tuple{Vector{Int64},Vector{Int64}},Expression}()
    function lsb(n)
      for i in 0:63
        if n & (1 << i) != 0
          return i
        end
      end
      return -1
    end
    for msize in 1:k
      lp2 = 1 << lsb(msize)
      for (mmi,(ix,jx)) in enumerate(product(subsets(1:b,msize),subsets(1:b,msize)))
        function minor(ix,jx,m,ms_dict)
          if length(ix) == 1
            return m[ix[1],jx[1]]
          end
          m1 = lp2 == msize ? div(lp2,2) : lp2
          return sum( [ (-1)^sum(indexin(e,ix)[1]-ei for (ei,e) in enumerate(ii)) *
                  ms_dict[(ii, jx[1:m1])] *
                  ms_dict[(setdiff(ix,ii), jx[m1+1:end])] for ii in subsets(ix,m1)] )
        end
        ms_dict[(ix,jx)] = minor(ix,jx,m,ms_dict)
      end
    end
    return ms_dict
  end

  sp = div(length(alpha),2)

  # lastnz = sp # change to last nonzero of alpha up to sp
  lastnz = maximum([0,[i for i in 1:sp if alpha[i] != 0]...])
  if lastnz > 0
    ms_dict = minors_up_to(m,lastnz)
    for (msize,k) in enumerate(alpha[1:lastnz])
      if k > 0
        ms = vec([ms_dict[(ix,jx)] for (ix,jx) in product(subsets(1:b,msize),subsets(1:b,msize))])
        @var mv[msize,1:length(ms),1:k]
        append!(eqs,vec(sum(mv.*ms,dims=1)))
        append!(p,vec(mv))
        mv0 = linear_forms_vanishing_on_prefix(hcat([Float64.(evaluate(ms,z=>z0)) for z0 in z0s]...),k)
        append!(p0,vec(transpose(mv0)))
      end
    end
  end

  firstnz = minimum([length(alpha)+1,[i for i in sp+1:length(alpha) if alpha[i] != 0]...])
  if firstnz < length(alpha)+1
    ms_dict = minors_up_to(y,b-firstnz)
    for (msize,k) in enumerate(alpha[end:-1:firstnz])
      if k > 0
        ms = vec([ms_dict[(ix,jx)] for (ix,jx) in product(subsets(1:b,msize),subsets(1:b,msize))])
        @var mv[b-msize,1:length(ms),1:k]
        append!(eqs,vec(sum(mv.*ms,dims=1)))
        append!(p,vec(mv))
        mv0 = linear_forms_vanishing_on_prefix(hcat([Float64.(evaluate(ms,z=>z0)) for z0 in z0s]...),k)
        append!(p0,vec(transpose(mv0)))
      end
    end
  end

  C = System(eqs,variables=z,parameters=p)
  sols = count(z0 -> norm(Float64.(evaluate(C.expressions,C.variables => z0, 
                                        C.parameters => p0))) < 1e-5,z0s)
  # println("starting with ", sols, " solutions")
  return nsolutions(monodromy_solve(C,z0s,p0,compile=compile,show_progress=show_progress, 
                                    max_loops_no_progress=max_loops_no_progress,unique_points_atol=xtol))
end

function linear_forms_vanishing_on_prefix(vs,k)
  return transpose(qr(hcat(vs,randn(size(vs)[1],k))).Q[:,end-k+1:end])
end

# size(T) = (a,b,b), len(alpha) = b-1, sum(alpha) = a-1
function characteristic_number_startsystem(T,alpha;compile=false,show_progress=false)
  a,b,_ = size(T)
  @var x[1:a]
  eqs = [sum(randn(length(x)) .* x) .- 1]
  M(e) = dropdims(sum(e.*T,dims=1),dims=1)
  m = M(x)
  for (msize,k) in enumerate(alpha)
    if k > 0
      ms = vec([det(m[ix,jx]) for (ix,jx) in product(subsets(1:b,msize),subsets(1:b,msize))])
      append!(eqs,vec(sum(randn(length(ms),k).*ms,dims=1)))
    end
  end
  C = System(eqs)
  sol = solve(C,compile=compile,show_progress=show_progress,tracker_options=(parameters=:conservative,))
  
  # newton_it(e) = e-jacobian(C,e)\evaluate(C,e)

  logsvals = [map(e->(f=svd(M(e)); log10(f.S[end-1] / f.S[1])) ,solutions(sol))..., -15]
  logsvals = filter(e->e >= -15,logsvals)
  sort!(logsvals,rev=true)
  cv = map(i-> (e=logsvals[i]; f=logsvals[i+1]; gap=e-f; 
                max(e-(-7),0) - max(f-(-7),0) + 2*gap),
                    1:length(logsvals)-1)
  # println(cv)
  res = if (length(cv) > 0) argmax(cv) else 0 end

  # println(res,"/",nsolutions(sol))
  # svals = sort(map(e->(f=svd(M(e)); f.S[end-1] / f.S[1]) ,solutions(sol)),rev=true)
  # println( svals[max(1,res-3):res] )
  # println( svals[res+1:min(end,res+4)] )
  return res
end

# size(T) = (a,b,b), len(alpha) = b-1, sum(alpha) = a-1
# these are mateusz's different notion of intersection_number, it doesnt compute
# the same thing as above!
function relative_characteristic_number(T,alpha;compile=false,show_progress=false)
  a,b,_ = size(T)
  @var x[1:a]
  eqs = [sum(randn(length(x)) .* x) .- 1]
  M(e) = dropdims(sum(e.*T,dims=1),dims=1)
  m = M(x)
  D_memo = Dict{Array{Int64,1},Expression}()
  D_memo[[]] = det(m)
  function D(beta)
    if !haskey(D_memo,beta)
      p = D(beta[1:end-1])
      D_memo[beta] = differentiate(p,x[beta[end]])
    end
    return D_memo[beta]
  end

  for (d,k) in enumerate(alpha)
    ds = [ D([e-i+1 for (i,e) in enumerate(ixs)]) for ixs in subsets(1: a + (b-d) - 1, b-d) ]
    # println(k," ",[([e-i+1 for (i,e) in enumerate(ixs)] ,e) for (ixs,e) in zip(subsets(1:a+(b-d)-1,b-d),ds) ] )
    # println(k," ",ds)
    append!(eqs,vec(sum(randn(length(ds),k).*ds,dims=1)))
  end

  display(eqs)
  C = System(eqs)
  sol = solve(C,compile=compile,show_progress=show_progress,tracker_options=(parameters=:conservative,))

  grad = [ D([i]) for i in 1:a ]
  logns = [[log10(norm(evaluate(grad,x => e))) for e in solutions(sol)]..., -15]
  logns = filter(e->e>=-15, logns)
  sort!(logns,rev=true)
  cv = map(i-> (e=logns[i]; f=logns[i+1]; gap=e-f; 
                max(e-(-7),0) - max(f-(-7),0) + 2*gap),
                    1:length(logns)-1)
  res = if (length(cv) > 0) argmax(cv) else 0 end

  # println(res,"/",nsolutions(sol))
  # ns = sort([norm(evaluate(grad,x => e)) for e in solutions(sol)],rev=true)
  # println( ns[max(1,res-3):res] )
  # println( ns[res+1:min(end,res+4)] )
  return res
end
  
function chromatic_polynomial(T)
  a,b,_ = size(T)
  characteristic_number = characteristic_number_generator(T)
  return [characteristic_number(alpha=[a-1-k , zeros(Int64,b-3)..., k]) for k in 0:a-1]
end
  
function relative_chromatic_polynomial(T)
  a,b,_ = size(T)
  characteristic_number = characteristic_number_generator(T)
  return [characteristic_number(beta=[a-1-k , zeros(Int64,b-3)..., k]) for k in 0:a-1]
end
  
export chromatic_polynomial_old
function chromatic_polynomial_old(T)
  a, b, _ = size(T)
  return [characteristic_number_old(T,[a-1-k,zeros(Int64,b-3)...,k]) for k in 0:a-1]
end
  
export relative_chromatic_polynomial_old
function relative_chromatic_polynomial_old(T)
  a,b,_ = size(T)
  return [relative_characteristic_number(T,[a-1-k , zeros(Int64,b-3)..., k]) for k in 0:a-1]
end

function characteristic_numbers(T)
  a,b,_ = size(T)
  # degree a-1 in b-1 variables
  out = []
  for ixs in subsets(1:b-1+a-1-1, a-1)
    alpha = zeros(Int64,b-1)
    for (j,i) in enumerate(ixs)
      alpha[i-j+1] += 1
    end
    push!(out,(alpha,characteristic_number(T,alpha)))
    println(join(map(string,out[end])," "))
  end
  return out
end

function relative_characteristic_numbers(T)
  a,b,_ = size(T)
  # degree a-1 in b-1 variables
  out = []
  for ixs in subsets(1:b-1+a-1-1, a-1)
    alpha = zeros(Int64,b-1)
    for (j,i) in enumerate(ixs)
      alpha[i-j+1] += 1
    end
    push!(out,(alpha,relative_characteristic_number(T,alpha)))
    # println(join(map(string,out[end])," "))
  end
  return out
end

include("misc.jl")
include("mindegree.jl")

end # module CharacteristicNumbers
