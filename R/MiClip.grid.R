MiClip.grid <-
function(a,b,c)
{
  b0=b[a==0]
  a1=a[a>0]
  b1=b[a>0]
  
  min_like=-1e300
  
  for (i in 1-c(1:20)*2*length(b1)/length(b))
  {
    if (i<0.5) {break}
    for (j in c(1:40)*c/200)
    {
      like=sum(log(i+(1-i)*(1-j)^b0))+sum(log((1-i)*j^a1*(1-j)^b1))
      if (like>min_like)
      {
        min_like=like
        pstr0=i
        prob=j
      }
    }
  } 
  
  fit=c(pstr0,prob)
  return(fit)
}
