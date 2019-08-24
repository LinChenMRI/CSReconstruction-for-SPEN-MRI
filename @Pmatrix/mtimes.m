function res = mtimes(a,b)
% res = mtimes(Pmatrix, x)

if a.adjoint
    res=a.P'*b;
%     res=a.P\b;
else
    res=a.P*b;
end
    
