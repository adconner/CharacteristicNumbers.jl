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
