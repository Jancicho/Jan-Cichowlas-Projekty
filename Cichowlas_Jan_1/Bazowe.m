function [x, d] = Bazowe(A, b)
% To samo, co CB_mod, tyle że funkcjami bazowymi

x = A\b;
d = det(A);
end