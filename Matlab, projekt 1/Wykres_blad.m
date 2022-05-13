size = 3:200;
blad_rozw = zeros(1, 198);
blad_wyznacznik = zeros(1, 198);

for k = 3:200
    v = 10*ones(1, k);
    w = -ones(1, k-1);
    A = diag(v) + diag(w, -1) + diag(w, 1);
    b = transpose(1:k);
    [x_1, d_1] = CB_mod(A, b);
    [x_2, d_2] = Bazowe(A, b);
    r_x = x_1 - x_2;
    r_d = d_1 - d_2;
    blad_rozw(k-2) = norm(r_x);
    blad_wyznacznik(k-2) = norm(r_d);
end

figure, 
subplot(1, 2, 1)
plot(size, blad_rozw, "g-")
title("Błędy rozwiązań dla macierzy różnych rozmiarów")
xlabel("rozmiar macierzy")
ylabel("bład")
subplot(1, 2, 2)
plot(size, blad_wyznacznik, "r*")
set(gca, 'YScale', 'log')
title("Błędy wyznaczników dla macierzy różnych rozmiarów")
xlabel("rozmiar macierzy")
ylabel("bład (skala logarytmiczna)")
