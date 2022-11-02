% Ejercicio 6 de la guía 3 de Métodos Numéricos en Fluidos I
% Chehade Pablo

% El objetivo es resolver una ecuación en derivadas parciales numéricamente a través de RK4 y algún método espacial, analizando posteriormente las características de los mismos

a = 1; %Este statement es necesario para que octave no piense que este archivo es un archivo de funciones

function u0_ = u0_CI(x)
    %Condición inicial: velocidad en x a t = 0 
    u0_ = exp(-200*(x-0.25)^2);
endfunction

function c_ = c(x, inciso)
    %Velocidad del sonido en función de la posición para el inciso inciso
    switch (inciso)
        case "a"
            %Arenisca porosa
            c_ = 1;
        case "b"
            %Transición a arenisca impermeable
            c_ = 1.25 - 0.25*tanh(40*(0.75-x));
        case "c"
            %Arenisca impermeable
            c_ = 1.5;
        case "d"
            %Nave alienígena sepultada
            c_ = 1.5 - exp(-300*(x-1.75)^2);
        case "e"
            %Verificación de la CB para c'(0) != 0
            c_ = x^2 + x + 1;
        otherwise
            error("La velocidad c no fue asignada")
    endswitch
endfunction


% El método RK4 se encuentra definido en RK4.m


function z0 = z0_ini(M)
    %Inicializo el vector z0
    L = 4; %Longitud del intervalo espacial
    Delta_x = L/(M+1);

    z0 = zeros(2*M + 2,1);
    for ii = 1:(M + 2)
        z0(ii) = u0_CI((ii-1)*Delta_x);
    endfor


endfunction



%Calculo u(x,t) para inciso
calculo = false;

if calculo == true


    inciso = "d";
    t_ini = 0;
    t_max = 8;

    %Discretización
    M = 390;
    N = 1600;

    % M = 4/Delta_x - 1
    % N = 8/Delta_t

    L = 4; %Longitud del intervalo espacial
    plotear = true;
    intervalo_plot = 0.00001;

    z_matriz = solucion(inciso, t_ini, t_max, M, N, L, plotear, intervalo_plot);

endif

#Calculo u(x,t) para cada inciso todos con distinta discretización y exporto los datos

calculo = false;

if calculo == true

    incisos_array = ["a", "a", "b", "c", "d"];

    #Tiempos inicial y final
    t_ini = 0;
    t_max = 8;

    %Discretización
    M_array = [3900, 390, 390,390, 390];
    N_array = [5300, 1600, 1600,1600, 1600];

    L = 4; %Longitud del intervalo espacial
    plotear = false;
    intervalo_plot = 0.0001;

    for kk = 1:length(incisos_array)
        inciso = incisos_array(kk);
        M = M_array(kk);
        N = N_array(kk);
        #Calculo la solución para el inciso inciso
        [t_array, z_matriz] = solucion(inciso, t_ini, t_max, M, N, L, plotear, intervalo_plot);
        #Guardo los datos
        datos = z_matriz;
        archivo = strcat("graficos/datos/solucion_",incisos_array(kk),"_M", num2str(M),".csv");
        csvwrite(archivo, datos);
    endfor

    
endif



















%Resulta que la condición en c'(0) no afecta que la onda pueda salir del dominio sin reflejarse! Raro

%Cuán bien está definido el límite dado por el nro de onda modificado?
%Para ello, propongo un conjunto de valores de N y busco M tal que la solución diverja mediante el método de bisección?
%El límite dependerá de c(x)?

% Dado Delta_x = L/(M+1), busco que Delta_t = (t_max-t_min)/(N) varié en (0.1, 2)*Delta_x
% Entonces, para dado M = L/Delta_x - 1, busco N = (t_max-t_min)/Delta_t tal que Delta_t varíe entre (0.1, 2)*Delta_x. Para cada Delta_t me dijo si la solución diverge

umbral = 2;
function diverge = diverge(z, umbral)
    %Me fijo si para dado z la función diverge
    %Defino como "divergencia" aquella situación en la que u(x,t) supera cierto valor umbral
    for ii=1:(length(z)/2 + 1) %Así lo que está dentro me da M+2
        if abs(z(ii)) > umbral
            diverge = true;
            return
        endif
    endfor
    %Si no diverge, asigno false
    diverge = false;

endfunction


calcular = false;

if calcular == true
    incisos_array = ["a", "b", "c", "d"];

    %Defino el intervalo de variación de Delta_x y Delta_t
    size_ = 200;
    Delta_x_lim = 0.25; %Límite máximo necesario para que la condición inicial esté bien representada
    Delta_x_array = linspace(1/size_,Delta_x_lim,size_-1);
    Delta_t_array = 2*Delta_x_array;

    for kk=1:4
        inciso = incisos_array(kk)

        t_ini = 0;
        t_max = 8;
        L = 4; %Longitud del intervalo espacial


        %Creo el vector de booleanos que dirá para dado N y M si diverge o no
        m_test = length(Delta_x_array);
        n_test = m_test; %Nro de test para dado M. Quiero que quede una matriz cuadrada

        divergencia = zeros(m_test,n_test);
        M_array = round(L./Delta_x_array - 1);
        N_array = round((t_max-t_ini)./Delta_t_array);

        % Delta_x_array = L./(M_array+1);
        % Delta_t_array = 2*Delta_x_array;

        contador = 0;
        for jj = 1:m_test
            for ii = 1:n_test
                contador = contador + 1;
                porcentaje = (contador)/(m_test*n_test)*100
                
                M = M_array(jj);
                N = N_array(ii); % round((t_max-t_ini)/Delta_t);
                %Calculo Delta_x
                Delta_x = Delta_x_array(jj);
                %Defino el intervalo en el que variará Delta_t
                Delta_t = Delta_t_array(ii);
                %Para cada Delta_t calculo N
                divergencia(jj,ii) = busco_divergencia(inciso, t_ini, t_max, M, N, L, umbral);

            endfor
        endfor

        %Guardo los datos
        datos = [Delta_x_array.', Delta_t_array.'];
        archivo = strcat("graficos/datos/divergencia_MyN_",incisos_array(kk),".csv");
        csvwrite(archivo, datos);

        datos = [divergencia];
        archivo = strcat("graficos/datos/divergencia_matrizbool_",incisos_array(kk),".csv");
        csvwrite(archivo, datos);
    endfor


    #Guardo los valores máximos de c(x) para todos los incisos
    L = 4;
    x_array = linspace(0,L,10000);
    cmax_array = zeros(length(incisos_array),1);
    for ii = 1:length(incisos_array)
        inciso = incisos_array(ii);
        c_array = zeros(length(x_array),1);
        for jj  = 1: length(x_array)
            c_array(jj) = c(x_array(jj), inciso);
        endfor
        cmax_array(ii) = max(c_array);
    endfor

    %Guardo datos
    datos = cmax_array
    csvwrite("graficos/datos/divergencia_cmax.csv", datos);
endif