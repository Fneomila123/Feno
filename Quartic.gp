allocatemem(50000000);
pbase = 5;
N = 7; 
acc = N+1; 
\\ Basic functions
padic(x) = x + O(pbase^(acc)) 
pu_pow(x,s) = exp(s * log(x)) 
\\ Cubic field parameters
f= 89;
if(gcd(pbase,f*eulerphi(f))!=1, print("SHIT!!!!"); break;,);
\\padic(x)=x + O(pbase^acc)
g = (teichmuller(padic(2)));
chi(n) = {if(n==1, return(padic(1)));
for(i=1,f-1, if((lift(znprimroot(f))^i)%f== n, return(g^i); break;));
}
\\chi(1) = 1+O(pbase^(acc));
sumtrunc = floor(N/eulerphi(2*f))+1;
chiarray = vector(2*f);
chivalues = vector(f);
for(i=1, f,chiarray[i]=chi(i); chiarray[i+f]=chiarray[i]; chivalues[i]=chiarray[i];);
\\Print out the values of chi
\\print(chivalues);
\\ The a_m(chi) coefficient vector
Aconst = -sum(j=1, f, j*chiarray[j] , 0)/f;
A(m)= {return(Aconst + sum(j=1, m-1, chiarray[j], 0) -
\end{verbatim}
\begin{verbatim}
 (2*sum(j=1, floor((m-1)/2), chiarray[j] , 0)) )}
a = vector(2*f);
a[1] = Aconst;
a[2] = Aconst+chiarray[1];
for(m=3, 2*f, a[m]=a[m-1] + chiarray[m-1] - 
\end{verbatim}
\begin{verbatim}
2*(floor((m-1)/2)-floor((m-2)/2))*chiarray[floor((m-1)/2)] );
q = floor(pbase^(sumtrunc*eulerphi(2*f)) / (2*f*pbase^N)); 
u = lift(Mod(1,2*f)/Mod(pbase,2*f)); \\ This is "varpi" in my paper
th(x, m) = lift( Mod(m,2*f*(pbase^N))+Mod((x-m)*((pbase*u)^N)-1 , 2*f*(pbase^N)) )+1;
mom(m) = sum(x=1, 2*f, if(th(x,m)<
(pbase^(sumtrunc*eulerphi(2*f))
-2*f*(pbase^N)*q), a[x] , padic(0) ) );
zetap(s, beta) = 
{return((sum(m=1, pbase^N, padic(if(gcd(pbase,m)==1,
mom(m)*(teichmuller(padic(m))^beta) * pu_pow(padic(m),-s), padic(0) )) , 0)));}
\\ Truncation of log(m)/log(1+p) modulo p^N
lambdaN(m) = lift(Mod((log(m + O(pbase^(N+1)))/log(pbase+1 + O(pbase^(N+1)))), pbase^N))
\\ Work out the j-th coefficient of the polynomial approximation modulo J_N
coeffappr(j, beta) = sum(m=1, pbase^N, 

if( gcd(pbase,m)==1 && j<=lambdaN(m), padic(binomial(lambdaN(m),j)) *

teichmuller(padic(m))^(beta) * mom(m), padic(0) ) );
\\ The interpolation polynomials F_N^(beta) -- might take ages to run!
polyzeta(beta) = sum(m=1, pbase^N, if(gcd(pbase,m)==1, mom(m)*
teichmuller(padic(m))^(beta) * (1+X)^(lambdaN(m)), 0));
laminv(beta) = {
for(i=0, 100, if(valuation(coeffappr(i, beta),pbase)==0, 
return(i-valuation(gcd(pbase,2^(beta+1)-1),pbase) ); break));
}
print("lambda(beta) = "laminv(1));
print("******************Voici les coef de F sont*************************")
pli= vector(11);
for(i=1,11, pli[i]=(coeffappr(i-1,1)));
print(pli);
laminv = laminv(1);
CapK = 10;
littleK = floor(CapK/laminv);
acc = CapK+1; \\ number of p-adic places to use
\\ Basic functions
padic(x) = x + O(pbase^(acc));
\\ Coefficients of f(T), to be0 inputted by hand, unfortunately
\\ Work out the coefficients b_n modulo p
bmodp = vector(CapK-laminv+1);
a=pli;
bmodp[1] = Mod(1,pbase)/a[1+laminv];
for(s=1, CapK-laminv, bmodp[s+1] = Mod(-1,pbase) * 
sum(i=1, s, a[1+laminv+i] * bmodp[1+s-i], 0) / a[1+laminv]);
\\ Work out the characteristic zero coefficients b_n
bvecp = vector(CapK-laminv+1);
for(i=0, CapK-laminv, bvecp[i+1] = padic(lift(bmodp[i+1])) );
\\ Iteratively compute better and better bvecp's
for(j=1, 30, bvecp[1] = (1/a[1+laminv]) * 
(1-sum(j=1, laminv, a[1+laminv-j]*bvecp[1+j],0)));
for(s=1, CapK-2*laminv, bvecp[1+s] = (-1/a[1+laminv]) * 
(sum(i=0, laminv+s, a[1+i]*bvecp[1+laminv+s-i], 0) - a[1+laminv]*bvecp[1+s])); 
\\ Compute c's from the a's and b's
cvecp = vector(1+laminv);
cvecp[1+laminv] = 1 +O(pbase^(littleK));
for(n=0, laminv-1, cvecp[n+1] = sum(i=0, n, a[i+1]*bvecp[n-i+1], 0)+O(pbase^(littleK)) );
 print(cvecp);
print(-cvecp[2]);
\end{verbatim}
