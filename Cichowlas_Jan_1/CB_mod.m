function [x, d] = CB_mod(A, b)
% Funkcja rozwiązuje układ równań Ax = b zmodyfikowaną metodą 
% Cholesky'ego- Banachiewicza oraz podaje wyznacznik maceirzy A.
% Autor - Jan Cichowlas
% Funkcja przyjmuje jako argumenty macierz A wymiaru nxn 
% oraz wektor b długości n.
% Funkcja zwraca rozwiązanie x równania Ax = b.

A_diag0 = diag(A, 0); % v przechowuje główną przekątną A
A_diag1 = diag(A, 1); % w przechowuje drugą przekątną A
% Uwaga: v i w zawierają w sobie wszystkie informacje na temat macierzy
% trójdiagonalnej symetrycznej A

n = length(A_diag0); % wymiar macierzy
U_diag0(n) = sqrt(A_diag0(1)); % przechowuje główną przekątną macierzy U
U_diag1(n-1) = A_diag1(1) / U_diag0(n); % przechowuje drugą przekątną macierzy U

% Algorytm jest analogiczny zwykłemu, temu stosowanemu w rozkładzie LL*.
% Musimy tylko "odbierać" elementy w odwrotnej kolejności. Dokładny opis
% modyfikacji standardowej metody w raporcie.

for k = 2:n
    U_diag0(n-k+1) = sqrt(A_diag0(k) - U_diag1(n-k+1)^2);
    if k == n
        break
    end
    U_diag1(n-k) = A_diag1(k) / U_diag0(n-k+1);
end

% Mając U, rozwiążmy równanie Uy = b:
y(n) = b(n) / U_diag0(n);
for k = (n-1):-1:1
    y(k) = (b(k) - U_diag1(k) * y(k+1)) / U_diag0(k);
end

% Mając U* oraz y, rozwiążmy równanie U*x = y:
x = zeros(n, 1); % aby x nie zmieniał rozmiaru w każdej iteracji
x(1) = y(1) / U_diag0(1);
for k = 2:n
    x(k) = (y(k) - U_diag1(k-1) * x(k-1)) / U_diag0(k);
end

% Pozostało jedynie policzyć żądany wyznacznik
d = prod(U_diag0)^2; % wyznacznik A to kwadrat iloczynu głównej przekątnej U

end
