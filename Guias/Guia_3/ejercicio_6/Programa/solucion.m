
function [t_array, z_matriz] = solucion(inciso, t_ini, t_max, M, N, L, plotear, intervalo_plot)
    %Resuelve la ecuación de onda para el inciso inciso y devuelve una matriz en la que cada fila corresponde a un instante de tiempo (tamaño N + 1) y cada columna a una posición espacial
    %t_ini: tiempo inicial
    %t_max: tiempo máximo
    %M: cantidad de puntos en el intervalo espacial
    %N: cantidad de puntos en el intervalo temporal
    %L: longitud del intervalo espacial
    %plotear: si es true, plotea la solución a todo tiempo
    %intervalo_plot: intervalo de tiempo entre cada plot

    %Discretización
    Delta_x = L/(M+1)
    Delta_t = (t_max - t_ini)/N

    %Creo el vector solución de tamaño 2N + 2
    z = zeros(2*M + 2,1);
    z_matriz = zeros(N+1, 2*M + 2);

    %Creo el vector temporal
    t_array = zeros(N+1);
    for ii = 1:(N+1)
        t_array(ii) = (ii-1)*Delta_t;
    endfor


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

    #Defino arrays horizontales para el gráfico
    xu_array = zeros(length(z(1:M+2)),1);
    xv_array = zeros(length(z(M+3:end)),1);
    for ii = 1:length(xu_array)
        xu_array(ii) = (ii-1)*Delta_x;
    endfor
    for ii = 1:length(xv_array)
        xv_array(ii) = (ii)*Delta_x;
    endfor

    %Defino arrays para graficar c(x)
    x_array = linspace(0,L,M+2);
    c_array = zeros(1,M+2);
    for ii = 1:(M+2)
        c_array(ii) = c(x_array(ii), inciso);
    endfor

    for jj = 1:N
        %Grafico
        
        if plotear == true
            plot(xu_array, z(1:M+2)); hold on;
            % plot(xu_array, z(1:M+2), "b", xv_array, z(M+3:end), "r"); hold on;

            % Grafico c(x)
            plot(x_array, c_array); hold off;

            % Mismo plot pero con legends. El ploteo es muchísimo más lento
            % plot(xu_array, z(1:M+2),";u(x,t);"); hold on;
            % plot(x_array, c_array, ";c(x);"); hold off;

            axis([0 4 -0.5 1.6])
            pause(intervalo_plot)




        endif

        t = t + Delta_t; 
        jj; #Imprimo el valor de jj
        z = RK4(z, t, Delta_t, @f_vec, gamma_);
        z_matriz(jj,:) = z.';
    endfor
endfunction