function diverge_ = busco_divergencia(inciso, t_ini, t_max, M, N, L, umbral)
    %Resuelve la ecuación de onda para el inciso inciso y devuelve true si la solución diverge o false si no
    %t_ini: tiempo inicial
    %t_max: tiempo máximo
    %M: cantidad de puntos en el intervalo espacial
    %N: cantidad de puntos en el intervalo temporal
    %L: longitud del intervalo espacial

    %Inicializo
    diverge_ = false;

    %Discretización
    Delta_x = L/(M+1);
    Delta_t = (t_max - t_ini)/N;

    %Creo el vector solución de tamaño 2N + 2
    z = zeros(2*M + 2,1);

    #Calculo la matriz gamma
    gamma_ = gamma_matriz(inciso, M);

    #Grafico la matriz
    % spy(gamma_matriz)
    % pause(2)

    %A partir de la matriz gamma defino la función f(z) tal que dz/dt = f(z)
    %Esta está definida como f(z) = gamma * z
    function dzdt = f_vec(z,t, gamma_)
        dzdt = gamma_*z;
    endfunction

    %Aplico el método RK4 y ploteo la solución u en cada paso de tiempo
    %Inicializo la condición inicial
    z = z0_ini(M);
    z_matriz(1,:) = z.'; 
    t = 0;

    for jj = 1:N
        % porcentaje = jj/N*100 %imprimo el porcentaje del cálculo
        t = t + Delta_t; 
        z = RK4(z, t, Delta_t, @f_vec, gamma_);

        %Me fijo si diverge
        if diverge(z, umbral) == true
            diverge_ = true;
            return
        endif
    endfor
endfunction