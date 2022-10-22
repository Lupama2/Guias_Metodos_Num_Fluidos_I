function gamma = gamma_matriz(inciso, M)
    %Calculo la matriz gamma necesaria para resolver la ecuación de onda

    L = 4; %Longitud del intervalo espacial
    Delta_x = L/(M+1);

    %Defino coeficientes para la resolución numérica
    a = c(0,inciso)/(2*Delta_x);
    b = c(4,inciso)/(2*Delta_x);
    d = 1/Delta_x^2;

    %Defino la matriz sparse gamma de tamaño (2M + 2) x (2M + 2)
    size_ = 2*M + 2;

    #1er cuadrante
    gamma_1_1 = sparse([1],[1],-3*a, size_,size_,0);
    gamma_1_2 = sparse([1],[2],4*a, size_,size_,0);
    gamma_1_3 = sparse([1],[3],-a, size_,size_,0);
    gamma_1_4 = sparse([M+2],[M],-b, size_,size_,0);
    gamma_1_5 = sparse([M+2],[M+1],4*b, size_,size_,0);
    gamma_1_6 = sparse([M+2],[M+2],-3*b, size_,size_,0);
    %Sumo
    gamma_1 = gamma_1_1 + gamma_1_2 + gamma_1_3 + gamma_1_4 + gamma_1_5 + gamma_1_6;

    #2do cuadrante
    gamma_2 = sparse([2:M+1],[1+ (M+2):M + (M+2)], ones(1,M), size_,size_,0);

    #3er cuadrante
    gamma_3 = sparse([1],[1],0,size_,size_,0); %Creo un sparse de ceros
    for ii = 1:M
        cii = c(ii*Delta_x, inciso);
        gamma_3_1 = sparse([ii + (M+2)],[ii], d*cii^2, size_,size_, 0);
        gamma_3_2 = sparse([ii + (M+2)],[ii+1],-2*cii^2*d, size_, size_, 0);
        gamma_3_3 = sparse([ii + (M+2)],[ii+2],d*cii^2, size_, size_, 0);
        gamma_3 = gamma_3 + gamma_3_1 + gamma_3_2 + gamma_3_3;
    endfor

    gamma_3;

    #4to cuadrante: son nulos

    gamma = gamma_1 + gamma_2 + gamma_3;

endfunction