function zfsum = zfgradsum(g)
[r,c] = size(g);
pop1 = zeros(r,c);
pop2 = zeros(r,c);
nfm = r/2+1;
mfm = c/2+1;
for n = 1:r;
    for m = 1:c;
            pop1(n,m) = (g(n,m))*abs( sin(2*pi*((n-nfm)/r))  );
            pop2(n,m) = (g(n,m))*abs( sin(2*pi*((m-mfm)/c)) );
    end
end
xf_sum = sum(sum(pop1));
yf_sum = sum(sum(pop1));
zfsum = xf_sum+yf_sum;


end