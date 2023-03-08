using HomotopyContinuation
using LinearAlgebra
using IterTools
using DataStructures

S23 = zeros(10,7,7)
S23[1,3,1] = -1
S23[1,1,3] = -1
S23[1,2,2] = 2
S23[2,1,5] = -1
S23[2,2,4] = 1
S23[2,4,2] = 1
S23[2,5,1] = -1
S23[3,1,6] = -1
S23[3,4,3] = 1
S23[3,3,4] = 1
S23[3,6,1] = -1
S23[4,1,6] = -1
S23[4,5,2] = 1
S23[4,2,5] = 1
S23[4,6,1] = -1
S23[5,1,7] = -1
S23[5,3,5] = 1
S23[5,5,3] = 1
S23[5,7,1] = -1
S23[6,5,5] = 2
S23[6,4,6] = -1
S23[6,6,4] = -1
S23[7,6,2] = 1
S23[7,7,1] = -1
S23[7,2,6] = 1
S23[7,1,7] = -1
S23[8,2,7] = -1
S23[8,3,6] = 1
S23[8,6,3] = 1
S23[8,7,2] = -1
S23[9,5,6] = 1
S23[9,7,4] = -1
S23[9,6,5] = 1
S23[9,4,7] = -1
S23[10,7,5] = -1
S23[10,5,7] = -1
S23[10,6,6] = 2

S14 = zeros(10,7,7)
S14[1,2,3] = 1
S14[1,4,1] = -1
S14[1,1,4] = -1
S14[1,3,2] = 1
S14[2,1,5] = -1
S14[2,2,4] = 1
S14[2,4,2] = 1
S14[2,5,1] = -1
S14[3,5,3] = -1
S14[3,3,5] = -1
S14[3,4,4] = 2
S14[4,1,6] = -1
S14[4,5,2] = 1
S14[4,2,5] = 1
S14[4,6,1] = -1
S14[5,4,5] = 1
S14[5,3,6] = -1
S14[5,6,3] = -1
S14[5,5,4] = 1
S14[6,3,7] = -1
S14[6,7,3] = -1
S14[6,5,5] = 2
S14[7,6,2] = 1
S14[7,7,1] = -1
S14[7,2,6] = 1
S14[7,1,7] = -1
S14[8,3,7] = -1
S14[8,7,3] = -1
S14[8,4,6] = 1
S14[8,6,4] = 1
S14[9,5,6] = 1
S14[9,7,4] = -1
S14[9,6,5] = 1
S14[9,4,7] = -1
S14[10,7,5] = -1
S14[10,5,7] = -1
S14[10,6,6] = 2

S22 = zeros(6,6,6)
S22[1,3,1] = -1
S22[1,1,3] = -1
S22[1,2,2] = 2
S22[2,1,5] = -1
S22[2,2,4] = 1
S22[2,4,2] = 1
S22[2,5,1] = -1
S22[3,1,6] = -1
S22[3,4,3] = 1
S22[3,3,4] = 1
S22[3,6,1] = -1
S22[4,1,6] = -1
S22[4,5,2] = 1
S22[4,2,5] = 1
S22[4,6,1] = -1
S22[5,6,2] = -1
S22[5,3,5] = 1
S22[5,5,3] = 1
S22[5,2,6] = -1
S22[6,5,5] = 2
S22[6,4,6] = -1
S22[6,6,4] = -1

S13 = zeros(6,6,6)
S13[1,2,3] = 1
S13[1,4,1] = -1
S13[1,1,4] = -1
S13[1,3,2] = 1
S13[2,1,5] = -1
S13[2,2,4] = 1
S13[2,4,2] = 1
S13[2,5,1] = -1
S13[3,5,3] = -1
S13[3,3,5] = -1
S13[3,4,4] = 2
S13[4,1,6] = -1
S13[4,5,2] = 1
S13[4,2,5] = 1
S13[4,6,1] = -1
S13[5,4,5] = 1
S13[5,3,6] = -1
S13[5,6,3] = -1
S13[5,5,4] = 1
S13[6,5,5] = 2
S13[6,4,6] = -1
S13[6,6,4] = -1

V2 = zeros(6,6,6)
V2[1,4,1] = -1
V2[1,1,4] = -1
V2[1,2,2] = 2
V2[2,2,3] = 1
V2[2,1,5] = -1
V2[2,3,2] = 1
V2[2,5,1] = -1
V2[3,1,6] = -1
V2[3,6,1] = -1
V2[3,3,3] = 2
V2[4,3,4] = 1
V2[4,4,3] = 1
V2[4,5,2] = -1
V2[4,2,5] = -1
V2[5,6,2] = -1
V2[5,3,5] = 1
V2[5,5,3] = 1
V2[5,2,6] = -1
V2[6,5,5] = 2
V2[6,4,6] = -1
V2[6,6,4] = -1


        

# L = zeros(15,10)
# L[1,1] = 1
# L[2,1] = 1
# L[3,1] = 1
# L[1,2] = -1
# L[4,2] = 1
# L[5,2] = 1
# L[4,3] = -1
# L[6,3] = 1
# L[7,3] = 1
# L[6,4] = -1
# L[8,4] = 1
# L[9,4] = 1
# L[2,5] = -1
# L[8,5] = -1
# L[10,5] = 1
# L[3,6] = -1
# L[11,6] = 1
# L[12,6] = 1
# L[5,7] = -1
# L[13,7] = 1
# L[14,7] = 1
# L[7,8] = -1
# L[11,8] = -1
# L[15,8] = 1
# L[9,9] = -1
# L[12,9] = -1
# L[13,9] = -1
# L[10,10] = -1
# L[14,10] = -1
# L[15,10] = -1
# L = transpose(L)

# function chromatic_number(L,k)
#   a,e = size(L)
#   @var x[1:a]
#   @var y[1:e]

#   eqs = reshape(reshape(x,1,a)*L, e) .* y - ones(e)

#   # startsols = 100
#   # x0s = randn(a,startsols)
#   # v0 = vec(svd(x0s).U[:,1]) 
#   # v0 can be any vector, pick one take make the next line not blow up x0's too much
#   # x0s ./= sum(v0 .* x0s,dims=1)
#   # x0s = [x0s[:,j] for j in 1:size(x0s,2)]

#   v0 = randn(a)
#   append!(eqs, [sum(v0 .* x) - 1])

#   append!(eqs, [sum(x.*randn(a)) for _ in 1:a-1-k ])
#   append!(eqs, [sum(y.*randn(e)) for _ in 1:k ])
#   C = System(eqs,variables=[x; y])
#   return C

# end

# # size(T) = (a,b,b), len(alpha) = b-1, sum(alpha) = a-1
# function top_intersection_number(T)
#   a,b,_ = size(T)
#   @var x[1:a]
#   c0 = randn(length(x))
#   eqs = [sum(c0 .* x) .- 1]
#   M(e) = dropdims(sum(e.*T,dims=1),dims=1)
#   ms(x) = vec([det(M(x)[ix,jx]) for (ix,jx) in product(subsets(1:b,b-1),subsets(1:b,b-1))])

#   @var c[1:length(ms(x))]
#   append!(eqs,[a-b for (a,b) in zip(c,ms(x))])
#   C = System(eqs;variables=[x...,c...])

#   sols = []
#   for i in 1:a-1
#     x0 = randn(length(x))
#     x0 ./= sum(c0 .* x0)
#     v = vec(ms(x0))
#     push!(sols,[x0...,v...])
#   end
#   L = collect(transpose(svd(hcat(sols...),full=true).U[:,length(sols)+1:end]))
#   L = LinearSubspace(L)

#   return C,sols,L
#   # s = monodromy_solve(C,sols,L)
#   return s
  
# end

# size(T) = (a,b,b), len(alpha) = b-1, sum(alpha) = a-1
function intersection_number_monodromy(T,alpha;startsols=1,xtol=1e-14,compile=true,show_progress=true,max_loops_no_progress=2)
  a,b,_ = size(T)
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
        mv0 = linear_forms_vanishing_on_prefix(hcat([Float64.(subs(ms,z=>z0)) for z0 in z0s]...),k)
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
        mv0 = linear_forms_vanishing_on_prefix(hcat([Float64.(subs(ms,z=>z0)) for z0 in z0s]...),k)
        append!(p0,vec(transpose(mv0)))
      end
    end
  end

  C = System(eqs,variables=z,parameters=p)
  sols = count(z0 -> norm(Float64.(subs(C.expressions,C.variables => z0, 
                                        C.parameters => p0))) < 1e-5,z0s)
  # println("starting with ", sols, " solutions")
  return nsolutions(monodromy_solve(C,z0s,p0,compile=compile,show_progress=show_progress, 
                                    max_loops_no_progress=max_loops_no_progress,unique_points_atol=xtol))
end

# size(T) = (a,b,b), len(alpha) = b-1, sum(alpha) = a-1
function intersection_number_monodromy_no_inv(T,alpha;startsols=1,xtol=1e-6,compile=false,show_progress=true,max_loops_no_progress=2)
  a,b,_ = size(T)
  @var x[1:a]
  x0s = randn(a,startsols)
  v0 = vec(svd(x0s).U[:,1]) 
  # v0 can be any vector, pick one take make the next line not blow up x0's too much
  x0s ./= sum(v0 .* x0s,dims=1)
  x0s = [x0s[:,j] for j in 1:size(x0s,2)]

  eqs = [sum(v0 .* x) .- 1]

  M(e) = dropdims(sum(e.*T,dims=1),dims=1)
  m = M(x)
  m0s = M.(x0s)
  params = Vector{Variable}()
  p0 = Vector{Float64}()
  vars = x
  vars0s = x0s
  ms_dict = Dict{Tuple{Vector{Int64},Vector{Int64}},Expression}()
  ms0_dicts = [Dict{Tuple{Vector{Int64},Vector{Int64}},Float64}() for i in 1:length(x0s)]
  function lsb(n)
    for i in 0:63
      if n & (1 << i) != 0
        return i
      end
    end
    return -1
  end
  for (msize,k) in enumerate(alpha)
    lp2 = 1 << lsb(msize)
    ms = Vector{Expression}()
    ms0s = [Vector{Float64}() for i in 1:length(x0s)]
    @var mm[msize,1:binomial(b,msize)^2] # may be unused
    for (mmi,(ix,jx)) in enumerate(product(subsets(1:b,msize),subsets(1:b,msize)))

      function minor(ix,jx,m,ms_dict)

        if length(ix) == 1
          return m[ix[1],jx[1]]
        end
        # return sum([ (-1)^(i+1)*m[e,jx[1]]*ms_dict[(vcat(ix[1:i-1],ix[i+1:end]), jx[2:end])] 
        #        for (i,e) in enumerate(ix)])
        m1 = lp2 == msize ? div(lp2,2) : lp2
        return sum( [ (-1)^sum(indexin(e,ix)[1]-ei for (ei,e) in enumerate(ii)) *
                ms_dict[(ii, jx[1:m1])] *
                ms_dict[(setdiff(ix,ii), jx[m1+1:end])] for ii in subsets(ix,m1)] )
      end

      mmcur = minor(ix,jx,m,ms_dict)
      mm0s = [minor(ix,jx,m0,ms0_dict) for (m0,ms0_dict) in zip(m0s,ms0_dicts)] 

      # this condition can be anything, decides when new variables for the
      # minors will be introduced
      if false
      # if msize == 2
      # if msize > 1 && msize == lp2 && rand() < 0.01
      # if msize > 1 && msize == lp2
        push!(eqs,mm[mmi] - mmcur)
        mmcur = mm[mmi]
        push!(vars,mm[mmi])
        for (vars0,mm0) in zip(vars0s,mm0s)
          push!(vars0,mm0)
        end
      end

      push!(ms,mmcur)
      ms_dict[(ix,jx)] = mmcur
      for (ms0,mm0) in zip(ms0s,mm0s)
        push!(ms0,mm0)
      end
      for (ms0_dict,mm0) in zip(ms0_dicts,mm0s)
        ms0_dict[(ix,jx)] = mm0
      end
    end
    if k > 0
      @var mv[msize,1:length(ms),1:k]
      append!(eqs,vec(sum(mv.*ms,dims=1)))
      append!(params,vec(mv))
      mv0 = linear_forms_vanishing_on_prefix(hcat(ms0s...),k)
      append!(p0,vec(transpose(mv0)))
    end
  end
  C = System(eqs,variables=vars,parameters=params)
  sols = count(x0 -> norm(Float64.(subs(C.expressions,C.variables => x0, 
                                        C.parameters => p0))) < 1e-5,vars0s)
  println("starting with ", sols, " solutions")
  return nsolutions(monodromy_solve(C,vars0s,p0,compile=compile,show_progress=show_progress, 
                                    max_loops_no_progress=max_loops_no_progress,unique_points_atol=xtol))
end

# size(T) = (a,b,b), len(alpha) = b-1, sum(alpha) = a-1
function intersection_number_monodromy_function_det(T,alpha;compile=true,show_progress=true,max_loops_no_progress=2,startsols=1)
  a,b,_ = size(T)
  @var x[1:a]
  x0s = randn(a,startsols)
  v0 = vec(svd(x0s).U[:,1]) 
  # v0 can be any vector, pick one take make the next line not blow up x0's too much
  x0s ./= sum(v0 .* x0s,dims=1)
  x0s = [x0s[:,j] for j in 1:size(x0s,2)]

  eqs = [sum(v0 .* x) .- 1]

  M(e) = dropdims(sum(e.*T,dims=1),dims=1)
  m = M(x)
  m0s = M.(x0s)
  params = Vector{Variable}()
  p0 = Vector{Float64}()
  for (msize,k) in enumerate(alpha)
    if k > 0
      Ms(m) = vec([poly_circuit_det(m[ix,jx]) for (ix,jx) in 
                   product(subsets(1:b,msize),subsets(1:b,msize))])
      ms = Ms(m)
      ms0s = Ms.(m0s)
      @var mv[msize,1:length(ms),1:k]
      append!(eqs,vec(sum(mv.*ms,dims=1)))
      append!(params,vec(mv))

      mv0 = linear_forms_vanishing_on_prefix(hcat(ms0s...),k)
      append!(p0,vec(transpose(mv0)))
    end
  end
  C = System(eqs,parameters=params)
  sols = count(x0 -> norm(Float64.(subs(C.expressions,C.variables => x0, 
                                        C.parameters => p0))) < 1e-5,x0s)
  println("starting with ", sols, " solutions")
  return nsolutions(monodromy_solve(C,x0s,p0,compile=compile,show_progress=show_progress, 
                                    max_loops_no_progress=max_loops_no_progress))
end

function linear_forms_vanishing_on_prefix(vs,k)
  return transpose(qr(hcat(vs,randn(size(vs)[1],k))).Q[:,end-k+1:end])
end

# size(T) = (a,b,b), len(alpha) = b-1, sum(alpha) = a-1
function intersection_number(T,alpha;compile=false,show_progress=false)
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
function intersection_number2(T,alpha;compile=false,show_progress=false)
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

function intersection_number_best_of(T,alpha,of=1)
  cnts = Dict{Int64,Int64}()
  while maximum([values(cnts)...,0]) < of
    cur = intersection_number_monodromy(T,alpha,startsols=1)
    cnts[cur] = get(cnts,cur,0)+1
  end
  if length(cnts) > of
    println("warning: needed ",length(cnts)," retries to get ",of," elements to agree: ",cnts)
  end
  return findfirst(p -> p == of, cnts)
end

function chromatic_numbers(T)
  a,b,_ = size(T)
  return [intersection_number_monodromy_no_inv(T, [a-1-k , zeros(Int64,b-3)..., k]) 
    for k in 0:a-1]
end

function chromatic_numbers2(T)
  a,b,_ = size(T)
  return [intersection_number2(T, [a-1-k , zeros(Int64,b-3)..., k]) for k in 0:a-1]
end

function characterestic_numbers(T)
  a,b,_ = size(T)
  # degree a-1 in b-1 variables
  out = []
  for ixs in subsets(1:b-1+a-1-1, a-1)
    alpha = zeros(Int64,b-1)
    for (j,i) in enumerate(ixs)
      alpha[i-j+1] += 1
    end
    push!(out,(alpha,intersection_number_monodromy_no_inv(T,alpha)))
    println(join(map(string,out[end])," "))
  end
  return out
end

function characteristic_numbers2(T)
  a,b,_ = size(T)
  # degree a-1 in b-1 variables
  out = []
  for ixs in subsets(1:b-1+a-1-1, a-1)
    alpha = zeros(Int64,b-1)
    for (j,i) in enumerate(ixs)
      alpha[i-j+1] += 1
    end
    push!(out,(alpha,intersection_number2(T,alpha)))
    # println(join(map(string,out[end])," "))
  end
  return out
end

function p_to_string(p)
  join([(c==1 ? "" : c==-1 ? "-" : string(c)*"*")*
        join(["t"*string(i)* (k==1 ? "" : "^"*string(k)) 
              for (i,k) in enumerate(alpha) if k > 0],"*") for (alpha,c) in p if c>0]," + ")
end

function generic_rankr(a,b,c,r)
  T = sum(reshape(randn(a),a,1,1).*reshape(randn(b),1,b,1).*reshape(randn(c),1,1,c) for i in 1:r)
end

function mateusz_tensor(n1,r1,n2,r2)
  T = generic_rankr(n2,n2,n2,r2)
  T[1:n1,1:n1,1:n1] += generic_rankr(n1,n1,n1,r1)
  return T
end

function go(b,a)
  p = nothing
  r = max(a,b)
  while true
    T = generic_rankr(a,b,b,r)
    # println(b," ",a," ",r)
    # q = characteristic_numbers(T)
    # q = characteristic_numbers2(T)
    q = p_to_string(intersection_polynomial(T))
    println(b," ",a," " ,r," ",q)
    if r > a && r > b && p == q
      break
    end
    p = q
    r += 1
  end
end

# like above, but prints results as they are obtained, and in an order more
# convenient for the computation (by r, then by alpha)
function go2(b,a)
  for ixs in subsets(1:a-1+b-1-1, b-1)
    alpha = zeros(Int64,b-1)
    for (j,i) in enumerate(ixs)
      alpha[i-j+1] += 1
    end

    r = max(a,b)
    last = nothing
    while true
      T = generic_rankr(a,b,b,r)
      num = intersection_number(T,alpha)
      println(b," ",a," ",r," ",alpha," ",num)
      if num == last
        break
      end
      last = num
      r += 1
    end
  end
end

# does like above, but for only top intersection number (which seems most
# discriminatory for r values)
function go3(b)
  r=Int64(ceil(3/2*b))
  last = nothing
  while true
    T = generic_rankr(b,b,b,r)
    alpha = [zeros(Int64,b-2)..., b-1]
    num = intersection_number(T,alpha)
    println(b," ",b," ",r," ",alpha," ",num)
    if num == last
      break
    end
    last = num
    r += 1
  end
end

function tensor_n(n)
  T = zeros(n,n,n)
  for i in 1:n
    T[i,i,i] = 1
  end
  return T
end

# makes sense when r >= a and r >= n-1
# (2,3,2) gives unstable numbers, mostly zero
# (k,n,n) gives (n-1 choose k-1)
function rank_number(a,n,r;startsols=32,xtol=1e-13)
  T = generic_rankr(a,n,n,r)
  while true
    sols = intersection_number_monodromy(T,[zeros(Int64,n-2)..., a-1],
                                        startsols=startsols,xtol=xtol)
    if sols > 0
      return sols
    end
    startsols = div(startsols,2)
    if startsols == 0
      break
    end
  end
  return sols
end

function rank_numbers(bnd)
  for n in 7:bnd
    # for a in 1:bnd
    for a in 1:2*n
    # for a in 1:n
    # for a in 1:n^2
      num = nothing
      r = max(a,n+1) # we know r==n
      # r = max(a,n-1)
      while true
        curnum = rank_number(a,n,r)
        println(a," ",n," ",r," ",curnum)
        # if r >= 2*n
        #   break
        # end
        if num == curnum
          break
        end
        num = curnum
        r += 1
      end
    end
  end
end

function rank_numbers_list(anrs)
  for (a,n,r) in anrs
    curnum = rank_number(a,n,r)
    println(a," ",n," ",r," ",curnum)
  end
end

# view these using permutedims(T,[2,3,1])
function mat_tensors()
  Ts = []
  T=zeros(4,3,3)
  T[1,1,1]=1
  T[2,1,2]=1
  T[2,2,1]=1
  T[3,2,2]=1
  T[4,3,3]=1
  push!(Ts,T)
  T=zeros(5,3,3)
  T[1,1,1]=1
  T[2,2,2]=1
  T[3,1,2]=1
  T[3,2,1]=1
  T[4,3,1]=1
  T[4,1,3]=1
  T[5,2,3]=1
  T[5,3,2]=1
  T[1:5,3,3].=1
  push!(Ts,T)
  T=zeros(8,4,4)
  T[1,1,1]=1
  T[2,2,1]=1
  T[2,1,2]=1
  T[3,2,2]=1
  T[4,3,2]=1
  T[4,2,3]=1
  T[5,3,3]=1
  T[6,4,4]=1
  T[7,3,4]=1
  T[7,4,3]=1
  T[8,1,4]=1
  T[8,4,1]=1
  push!(Ts,T)
  T=zeros(5,6,6)
  T[1,1,1]=1
  T[1,2,2]=-1
  T[2,2,2]=1
  T[2,3,3]=-1
  T[3,3,3]=1
  T[3,4,4]=-1
  T[4,4,4]=1
  T[4,5,5]=-1
  T[5,5,5]=1
  T[5,6,6]=-1
  push!(Ts,T)
  return Ts
end

function poly_circuit_det(m)
  n = size(m,1)
  mpows = [m]
  for i in 2:n
    push!(mpows,mpows[end]*m)
  end
  p = [tr(m) for m in mpows]
  principal_minors = [p[1]]
  # mm = [i < j-1 ? 0 : i == j-1 ? i : p[i-j+1] for (i,j) in product(1:n,1:n)]
  # result is det(mm), compute by laplace expansion
  for i in 2:n
    push!(principal_minors,
          sum((-1)^(i-j)*prod(j:i-1)*p[i-j+1]*(j > 1 ? principal_minors[j-1] : 1) for j in 1:i))
  end
  return principal_minors[end] # actually n!*det(m)
end

if !isinteractive()
  a,n,r = map(x->parse(Int64,x),ARGS)
  println(a," ",n," ",r," ",rank_number(a,n,r))
end

