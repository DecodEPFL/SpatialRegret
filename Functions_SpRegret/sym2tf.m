function Pz  = sym2tf(Pz_symbolic, Ts)
if nargin<2
    Ts = -1;
end

Pz_symbolic = simplify(Pz_symbolic);


n=size(Pz_symbolic,1);
m=size(Pz_symbolic,2);
Pz = tf(zeros(n,m));
for i=1:n
    for j=1:m
        [symNum,symDen] = numden(Pz_symbolic(i,j)); %Get num and den of Symbolic TF
        TFnum = sym2poly(symNum);    %Convert Symbolic num to polynomial
        TFden = sym2poly(symDen);    %Convert Symbolic den to polynomial
        Pz(i,j) = tf(TFnum/TFden(1),TFden/TFden(1), Ts);
    end
end
end