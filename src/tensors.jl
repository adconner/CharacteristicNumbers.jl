function generic_rankr(a,b,c,r)
  T = sum(reshape(randn(a),a,1,1).*reshape(randn(b),1,b,1).*reshape(randn(c),1,1,c) for i in 1:r)
end

function tensor_n(n)
  T = zeros(n,n,n)
  for i in 1:n
    T[i,i,i] = 1
  end
  return T
end

# view these using permutedims(T,[2,3,1])
function cm_tensors()
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
