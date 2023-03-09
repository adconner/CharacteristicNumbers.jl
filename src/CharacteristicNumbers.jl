module CharacteristicNumbers
export characteristic_number,chromatic_polynomial_coefficient,characteristic_numbers,chromatic_polynomial
export relative_characteristic_number,relative_chromatic_polynomial_coefficient,relative_characteristic_numbers,relative_chromatic_polynomial
export mat_tensors

using HomotopyContinuation
using LinearAlgebra
using IterTools
using DataStructures

# size(T) = (a,b,b), len(alpha) = b-1, sum(alpha) = a-1
function characteristic_number(T,alpha;startsols=1,xtol=1e-14,compile=true,
    show_progress=true,max_loops_no_progress=5)
  a,b,c = size(T)
  @assert(length(alpha) == b-1 && sum(alpha) == a-1 && b == c)
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
  
function chromatic_polynomial_coefficient(T, k)
  a,b,_ = size(T)
  return characteristic_number(T,[a-1-k , zeros(Int64,b-3)..., k])
end
  
function relative_chromatic_polynomial_coefficient(T, k)
  a,b,_ = size(T)
  return relative_characteristic_number(T,[a-1-k , zeros(Int64,b-3)..., k])
end

function chromatic_polynomial(T)
  a,b,_ = size(T)
  return [chromatic_polynomial_coefficient(T, k) for k in 0:a-1]
end

function relative_chromatic_polynomial(T)
  a,b,_ = size(T)
  return [relative_chromatic_polynomial_coefficient(T, k) for k in 0:a-1]
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

end # module CharacteristicNumbers
