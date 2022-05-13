size = 3:200;
t_mod = zeros(1, 198);
t_stand = zeros(1, 198);

for k = 3:200
    v = 10*ones(1, k);
    w = -ones(1, k-1);
    A = diag(v) + diag(w, -1) + diag(w, 1);
    b = transpose(1:k);
    tic
    CB_mod(A, b);
    t_mod(k-2) = toc
    tic
    Bazowe(A, b);
    t_stand(k-2) = toc
end

figure, 
plot(size, t_stand, "g-", size, t_mod, "r*")
legend("Wbudowane funkcje Matlaba", "Zmodyfikowana metoda C-B")
xlabel("rozmiar macierzy")
ylabel("czas")



