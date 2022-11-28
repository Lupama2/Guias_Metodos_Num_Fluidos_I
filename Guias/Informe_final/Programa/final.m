% Los parámetros a modificar en el tp final son

% * nombre de los archivos salientes. VERIFICAR si se reportan los valores de u y v para todo tiempo o solo al final, los necesito para todo tiempo
% 


%Testeo la función
% Re = 5000.0;
% n1 = 80;
% dt = 0.5;
% nsimpler = 1;
% termino_advectivo = "D";
% metodo_temporal = "EI";
% tol_estacionario = 1e-6;
% Ndeltat=100; #numero de pasos de tiempo
% archivo_velocidades_centrales = "graficos/datos/velocentral.csv";
% archivo_evolucion_variables = "graficos/datos/evolucion.csv";
% archivo_parametros = "graficos/datos/parametros.csv";

% function U0_top = U0_top(x)
%     U0_top = 1.0;
% endfunction

% Chehade_NSdc2(Re, @U0_top, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)




#Resultado típico:
% fid = 3
% ans = Tiempo de Calculo en segundos para
% Ndeltat = 80
% ans = iteraciones
% ans = 42.577
% ans = 0
% fid = 3
% ans = 0

#Luego de varias veces:
% fid = 3
% ans = Tiempo de Calculo en segundos para
% Ndeltat = 80
% ans = iteraciones
% ans = 42.019
% ans = 0
% fid = 3
% ans = 0

#Al inicio:
% fid = 3
% ans = Tiempo de Calculo en segundos para
% Ndeltat = 80
% ans = iteraciones
% ans = 40.959
% ans = 0
% fid = 3
% ans = 0

#Al inicio nuevamente:
% fid = 3
% ans = Tiempo de Calculo en segundos para
% Ndeltat = 80
% ans = iteraciones
% ans = 39.309
% ans = 0
% fid = 3
% ans = 0




1; #Esto está puesto para que octave no interprete que este es un archivo de funciones







function U0_top = U0_top_cte(t)
    U0_top = 1.0;
endfunction


#Inciso a:

% Verifico si el resultado depende de deltat para el caso particular Re=1000, n1 = 21. Para distintos valores de deltat, calculo u_central y v_central y los guardo en un archivo. Consideraciones
% * Uso Ndeltat lo suficientemente grande como para que converja con cualquier dt

calcular = false;
tol_estacionario = 1e-5;


if calcular == true
    "Inciso a"
    Re = 1000
    n1 = 20
    % dt_array = [1, 2]
    dt_array = [0.005, 0.01,0.02,0.05,0.1,0.2,0.5,1,2, 5, 10, 20];

    #Abro el archivo en el que voy a guardar los datos
    archivo_velcentral = "graficos/datos/a_velcentral.csv";
    fid = fopen(archivo_velcentral, "w");


    for ii=1:length(dt_array)
        dt = dt_array(ii)

        Ndeltat = 1000000;

        
        nsimpler = 1;
        termino_advectivo = "D"
        metodo_temporal = "E";

        archivo_velocidades_centrales = "graficos/datos/velocentral_basura_a.csv";
        archivo_evolucion_variables = "graficos/datos/evolucion_basura_a.csv";
        archivo_parametros = "graficos/datos/parametros_basura_a.csv";

        Chehade_NSdc2(Re, @U0_top_cte, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)

        #Abro el archivo de velocidades y extraigo los valores en el centro

        data = importdata(archivo_velocidades_centrales);
        posicion_centro_izq = (length(data(:,1)))/2
        posicion_centro_der = (length(data(:,1)))/2 + 1
        % data(posicion_centro,1)
        a_ucentral(ii) = (data(posicion_centro_izq, 2) + data(posicion_centro_der, 2))/2;
        a_vcentral(ii) = (data(posicion_centro_izq, 3) + data(posicion_centro_der, 3))/2;

        #Guardo el dato
        fprintf(fid, "%e %e %e\n", dt, a_ucentral(ii),a_vcentral(ii));

    endfor

    #Cierro el archivo
    fclose(fid);

    #Exporto los datos
    
    % datos = [dt_array.', a_ucentral.', a_vcentral.'];
    % csvwrite(archivo_velcentral, datos);
endif





#Inciso a2:

% Buscamos no perder tiempo y agarrar de prepo el mejor dt posible para cada Re y cada n1. Para ello, partimos de un dt y corremos con Nmax = 10. Si vemos que va a converger, dividimos dt por 2. Si no vemos que vaya a converger, lo multiplicamos por 2 y nos quedamos con ese. ¿Cómo definimos si "va a converger"? Si el cociente entre las derivadas de las velocidades y la tolerancia del estacionario disminuye en el tiempo. Para esto, calculamos tal cociente a cada tiempo y calculamos la pendiente del ajuste lineal. Si la pendiente es negativa, "va a converger"


function pendiente = pendiente_ajuste_lineal(t_array, y_array)
    %Calcula la pendiente del ajuste lineal de los datos
    %t_array: array de tiempos
    %y_array: array de valores de la variable que se quiere ajustar
    %pendiente: pendiente del ajuste lineal

    %Calculo la pendiente del ajuste lineal
    p = polyfit(t_array, y_array, 1);
    pendiente = p(1);

endfunction


function m = va_a_converger(Re, U0_top, n1, dt, tol_estacionario, Ndeltat_max, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)
    #Criterio de convergencia en base a las primeras Ndeltat_max evaluaciones

    #Hace la cuenta
    Chehade_NSdc2(Re, @U0_top_cte, n1, dt, tol_estacionario, Ndeltat_max, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)

    #Abro el archivo de evolución de variables y extraigo los valores de las derivadas de u y v a distintos tiempos. La pendiente se calcula como la suma cuadrática en función del nro de pasos.
    data = importdata(archivo_evolucion_variables);
    % t_array = data(:,1);
    dudt = data(:,3);
    dvdt = data(:,4);
    x_array = linspace (0, length(dudt), length(dudt));
    #Calculo para cada t la suma cuadrática
    for ii = 1:length(dudt)
        dveldt(ii) = (dudt(ii)^2 + dvdt(ii)^2)^0.5;
    endfor

    #Calculo la pendiente
    % x_array = x_array(2:length(x_array));
    % dveldt = dveldt(2:length(dveldt));
    m = pendiente_ajuste_lineal(x_array, log(abs(dveldt)))

    if m<0
        convergencia = true;
    else
        convergencia = false;
    endif

endfunction

% #Para testear lo aplicamos al caso anterior
% Re = 5000
% n1 = 80

% dt_guess = 4;
% Ndeltat = 20;
% tol_estacionario = 1e-6;

% #Calculo:
% nsimpler = 1;
% termino_advectivo = "D"
% metodo_temporal = "E";

% archivo_velocidades_centrales = "graficos/datos/velocentral_basura.csv";
% archivo_evolucion_variables = "graficos/datos/evolucion_basura.csv";
% archivo_parametros = "graficos/datos/parametros_basura.csv";

% va_a_converger(Re, @U0_top_cte, n1, dt_guess, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)

% #dt = 0.5: 21 s para llegar a 12.494
% #dt = 1: 21 s para llegar a 11.177
% #dt = 2: 21 s para llegar a 12.494
% #dt = 4: x s para llegar a 13.074





function best_dt = busco_best_dt(m_old, Re, U0_top, n1, dt_guess, tol_estacionario, Ndeltat_max, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)
    #Para dadas características de la resolución busco el mejor paso de tiempo (aquel que converja más rápido).
    #Sí o sí me tengo que asegurar de entrada que dt_guess sea lo suficientemente pequeño para que haya convergencia y que sea lenta

    #Me dijo si para dt_guess hay convergencia
    dt = dt_guess;
    m_new = va_a_converger(Re, U0_top, n1, dt, tol_estacionario, Ndeltat_max, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);

    if m_new < m_old
        #Si la pendiente nueva es más grande en módulo que la anterior, entonces duplicá dt y fijate devuelta
        dt_guess = dt*2
        best_dt = busco_best_dt(m_new, Re, U0_top, n1, dt_guess, tol_estacionario, Ndeltat_max, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)
    else
        best_dt = dt/2
    endif

endfunction

#Para testear lo aplicamos al caso anterior
% Re = 1000
% n1 = 21


% dt_guess = 0.005;
% Ndeltat = 10;
% tol_estacionario = 1e-6;

% #Calculo:
% nsimpler = 1;
% termino_advectivo = "D";
% metodo_temporal = "E";

% archivo_velocidades_centrales = "graficos/datos/velocentral_basura.csv";
% archivo_evolucion_variables = "graficos/datos/evolucion_basura.csv";
% archivo_parametros = "graficos/datos/parametros_basura.csv";

% busco_best_dt(0,Re, @U0_top_cte, n1, dt_guess, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)















function comparacion_Ghuia(Re_array, U0_top, n1_array, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_ucentral, archivo_vcentral, dt_guess, Ndeltat_max, buscar, inciso)
    #Esta función sirve para comparar los resultados de la función Chehade_NSdc2 con los resultados de la función Chehade_NSdc2_Ghuia para los parámetros de entrada de Chehade_NSdc2
    #Re_array: array con los valores de Reynolds a probar
    #U0_top: función que devuelve la velocidad en la parte superior del dominio
    #n1_array: array con los valores de n1 a probar
    #...
    #archivo_ucentral: archivo donde se guardan los valores de u en el centro del dominio para cada valor de Re y n1
    #archivo_vcentral: archivo donde se guardan los valores de v en el centro del dominio para cada valor de Re y n1

    #Para cada caso busca automáticamente el mejor dt si buscar == true. Caso contrario, usa directamente dt_guess.

    #Creo las matrices que guardarán los valores de U y V en el centro para cada Re y n1
    b_ucentral = zeros(length(n1_array), length(Re_array));
    b_vcentral = zeros(length(n1_array), length(Re_array));

    for ii = 1:length(n1_array)
        for jj=1:length(Re_array)
            n1 = n1_array(ii)
            Re = Re_array(jj)

            archivo_velocidades_centrales = strcat("graficos/datos/", inciso, "_",termino_advectivo, "_",num2str(n1), "_",num2str(Re), "_velocentral.csv");
            archivo_evolucion_variables = strcat("graficos/datos/", inciso,  "_",termino_advectivo, "_",num2str(n1), "_",num2str(Re),"_evolucion.csv");
            archivo_parametros = strcat("graficos/datos/", inciso, "_",termino_advectivo, "_",num2str(n1), "_",num2str(Re), "_parametros.csv");

            if buscar == true
            #Busco el mejor dt:
                dt = busco_best_dt(0,Re, @U0_top, n1, dt_guess, tol_estacionario, Ndeltat_max, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)
            else
                dt = dt_guess;
            endif

            Chehade_NSdc2(Re, @U0_top_cte, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);

            #Abro el archivo de velocidad central y extraigo los valores de U en y = 0.5 y de V en 
            data = importdata(archivo_velocidades_centrales);
            posicion_centro_izq = (length(data(:,1)))/2
            posicion_centro_der = (length(data(:,1)))/2 + 1
            data(posicion_centro_izq, 1)
            data(posicion_centro_der, 1)
            % data(posicion_centro,1)
            b_ucentral(ii,jj) = (data(posicion_centro_izq, 2) + data(posicion_centro_der, 2))/2;
            b_vcentral(ii,jj) = (data(posicion_centro_izq, 3) + data(posicion_centro_der, 3))/2;
        
        % fid_u = fopen(archivo_ucentral, "w");
        % fid_v = fopen(archivo_vcentral, "w");
        
        % fprintf(fid, "%e %e %e\n", b_ucentral(ii,:));
        % fprintf(fid, "%e %e %e\n", b_ucentral(ii,:));

        #Guardo ambas matrices en un .csv. Lo pongo adentro del for para que valla guardando después de cada Re
        csvwrite(archivo_ucentral, b_ucentral);
        csvwrite(archivo_vcentral, b_vcentral);
        
        endfor
    endfor


endfunction


#Inciso b

calcular = false;
tol_estacionario = 1e-5;

if calcular == true
    inciso = "b"
    #Parámetros que puedo llegar a modificar
    Ndeltat= 80000; #80000; #numero de pasos de tiempo

    % n1_array = [80];
    % Re_array = [5000];
    n1_array = [80];
    Re_array = [5000];#,1000,5000];
    archivo_ucentral = "graficos/datos/b_ucentral.csv";
    archivo_vcentral = "graficos/datos/b_vcentral.csv";

    #Parámetros que seguro no modifique
    nsimpler = 1;
    termino_advectivo = "D";
    metodo_temporal = "E";
    
    dt_guess = 2;#0.005;
    Ndeltat_max = 10;
    buscar = false;

    comparacion_Ghuia(Re_array, @U0_top_cte, n1_array, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_ucentral, archivo_vcentral, dt_guess, Ndeltat_max, buscar, inciso)
endif

% Test en inciso b para determinar dt: el objetivo es calcular dt tal que el programa converja para Re alto y n1 alto. Si esto ocurre, creo que va a converger para cualquier otro Re y n1.
% * Para dt = 0.5 llegó a 1.7185e5 en k = 71
% * Para dt = 1 llegó a 2.3871e4 en k = 71
% * Para dt = 2 llegó a 8169 en k = 71
% * Para dt = 4 llegó a 2.1021e+24 en k = 71 (DIVERGE)
% Resulta que el caso anterior diverge para n1 = 81 y Re = 100.
% A partir de ahora trabajo este último caso en particular
% * Para dt = 1 diverge
% * Para dt = 0.2 llega a 1.7650e+05 en k = 35. Al parecer converge
% Testeo que todos los casos con n1 alto converjan:
% * n1 = 81, Re = 1000 llega a 6.1018e+05 en k = 72
% * n1 = 81, Re = 5000 llega a 2.4054e+06 en k = 14
% Al parecer converge en todos los casos para esas condiciones.




#Inciso c
#Hago lo mismo que el b pero con distinto esquema para el término advectivo.
#AÚN NO PROGRAMÉ EL QUICK así que ejecuto dos veces con el UP1
calcular = false;
tol_estacionario = 1e-5;

if calcular == true
    inciso = "c"
    termino_advectivo_array = ["Q"];# ["U","Q"]

    #Parámetros que puedo llegar a modificar
    % dt = 2; #Dentro de comparacion_Guia se busca el mejor dt posible
    Ndeltat= 80000; #80000; #numero de pasos de tiempo

    n1_array = [80];#,80];
    Re_array = [100]#,1000,5000];

    #Parámetros que seguro no modifique
    nsimpler = 1;
    metodo_temporal = "E";

    dt_guess = 0.35;
    Ndeltat_max = 10;
    buscar = false;

    for kk = 1:length(termino_advectivo_array)

        termino_advectivo = termino_advectivo_array(kk);
        archivo_ucentral = strcat("graficos/datos/c_ucentral_",termino_advectivo,".csv");
        archivo_vcentral = strcat("graficos/datos/c_vcentral_",termino_advectivo,".csv");

        comparacion_Ghuia(Re_array, @U0_top_cte, n1_array, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_ucentral, archivo_vcentral, dt_guess, Ndeltat_max, buscar, inciso)

    endfor
endif

% Hago lo mismo para UP1 y QUICK.

% Test en inciso c UP1 para determinar dt: el objetivo es calcular dt tal que el programa converja para Re alto y n1 alto. Si esto ocurre, creo que va a converger para cualquier otro Re y n1.
% * Para dt = 0.1, n1 = 81, Re = 100 converge pero Re = 1000 no
% * Para dt = 0.05 ocurre lo mismo
% * Para dt = 0.01 ocurre lo mismo (inicialmente parece converger pero luego diverge)
% * Para dt = 0.01 ocurre lo mismo (inicialmente parece converger pero luego diverge)
% * Para dt = 0.001 converge al parecer. Aunque para correr esto debería dejar la compu toda la noche creo.





#Inciso d
calcular = false;
tol_estacionario = 1e-5;

if calcular == true

    #Defino el mejor esquema advectivo
    termino_advectivo_mejor_esquema = "D";

    #Defino parámetros
    Re_array = [1,1000];
    n1_array = [10,20,40,60,80];

    Ndeltat=80000; #numero máximo de pasos de tiempo

    #Parámetros que seguro no modifique
    nsimpler = 1;
    metodo_temporal = "E";
    
    buscar = true;
    dt_guess = 0.005;
    Ndeltat_max = 10;

    #Recorro ambos Re y para cada uno de ellos voy calculando los errores
    for ii = 1:length(Re_array)
        Re = Re_array(ii)


        #Defino los vectores en los que guardaré los errores
        error_uparticular = zeros(1,length(n1_array));
        error_vparticular = zeros(1,length(n1_array));

        #Defino los valores solución
        n1_sol = 80;
        % Los valores particulares u(0.5, 0.2) y v(0.5,0.2) los voy a guardar en las variables
        % c_usol_particular
        % c_vsol_particular

        #Defino los nombres de los archivos
        archivo_velocidades_centrales = "graficos/datos/d_velocentral_basura.csv";
        archivo_evolucion_variables = "graficos/datos/d_evolucion_basura.csv";
        archivo_parametros = "graficos/datos/d_parametros_basura.csv";

        #Hago la cuenta

        #Busco el mejor dt:
        if buscar == true
            dt = busco_best_dt(0,Re, @U0_top_cte, n1_sol, dt_guess, tol_estacionario, Ndeltat_max, nsimpler, termino_advectivo_mejor_esquema, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);
        else
            dt = dt_guess;
        endif

        Chehade_NSdc2(Re, @U0_top_cte, n1_sol, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo_mejor_esquema, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);


        #Calculo los valores u(0.5, 0.2) y v(0.5,0.2) interpolando los valores próximos

        #Abro el archivo de velocidad central
        data = importdata(archivo_velocidades_centrales);
      
        c_usol_particular = (data((n1_sol/5),2) + data((n1_sol/5) + 1,2) )/2;
        c_vsol_particular = (data((n1_sol/5),3) + data((n1_sol/5) + 1,3) )/2;


        #Ya tengo los valores de la solución, calculemos los errores para distintos n1
        
        #Recorro los distintos n1
        for jj=1:length(n1_array)
            n1 = n1_array(jj)
            #Hago la cuenta
            if buscar == true
                dt = busco_best_dt(0,Re, @U0_top_cte, n1, dt_guess, tol_estacionario, Ndeltat_max, nsimpler, 'U', metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);
            else
                dt = dt_guess;
            endif

            Chehade_NSdc2(Re, @U0_top_cte, n1, dt, tol_estacionario, Ndeltat, nsimpler, "U", metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);

            #Extraigo el valor de u(0.5, 0.2) y v(0.5,0.2) de la aproximación
            data = importdata(archivo_velocidades_centrales);
            c_uaprox_particular = (data((n1/5),2) + data((n1/5) + 1,2) )/2;
            c_vaprox_particular = (data((n1/5),3) + data((n1/5) + 1,3) )/2;

            #Calculo el error
            error_uparticular(jj) = abs(c_usol_particular - c_uaprox_particular);
            error_vparticular(jj) = abs(c_vsol_particular - c_vaprox_particular);
        endfor

        #Guardo los errores
        archivo_error_uparticular = strcat("graficos/datos/d_error_uparticular_Re",num2str(Re),".csv");
        archivo_error_vparticular = strcat("graficos/datos/d_error_vparticular_Re",num2str(Re),".csv");
        csvwrite(archivo_error_uparticular, error_uparticular);
        csvwrite(archivo_error_vparticular, error_vparticular);
    
    endfor
        
endif



#Inciso e
#ESTÁ CALCULANDO MAL LA VELOCIDAD CENTRAL


function ucentral = u_central_estacionario(Re, U0_top, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)
    #Esta función evoluciona el sistema hasta el estado estacionario y devuelve el valor de u en el centro en el estado estacionario

    Chehade_NSdc2(Re, U0_top, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)


    #Abro el archivo de velocidad central y extraigo los valores de U en y = 0.5 y de V en 
    data = importdata(archivo_velocidades_centrales);



    % data(posicion_centro,1)

    posicion_centro_izq = (length(data(:,1)))/2
    posicion_centro_der = (length(data(:,1)))/2 + 1

    ucentral = (data(posicion_centro_izq, 2) + data(posicion_centro_der, 2))/2;

endfunction 





#Testeo:

% Ndeltat=8000;

% Re = 1000;
% n1 = 21
% dt = 0.05; #0.001 Con este valor converge para DC2 y Re = 1000
% termino_advectivo = "D";
% tol_estacionario = 1e-3;

% nsimpler = 1;
% metodo_temporal = "EI";

% archivo_velocidades_centrales = "graficos/datos/velocentral_basura_e.csv";
% archivo_evolucion_variables = "graficos/datos/evolucion_basura_e.csv";
% archivo_parametros = "graficos/datos/parametros_basura_e.csv";

% u_central_estacionario(Re, @U0_top_cte, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)



function n1 = n1_minimo(Re, U0_top, n1_guess, dt_guess, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros, e_tol, u_central_guia)

    % Planteo un n1_guess. Sí o sí tengo que empezar de un valor para el que la solución ya converja.

    # Busco el mejor dt
    Ndeltat_max = 10;
    best_dt = busco_best_dt(0, Re, U0_top, n1_guess, dt_guess, tol_estacionario, Ndeltat_max, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)
    % best_dt = dt_guess #Si no quiero buscar el mejor dt

    # Evoluciono hasta el estado estacionario y me fijo si bajo esa condición se logra un error menor a la tolerancia
    ucentral = u_central_estacionario(Re, U0_top, n1_guess, best_dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)

    error_e = abs((ucentral - u_central_guia)/u_central_guia)

    if (error_e  < e_tol)
        "Disminuyo n1_guess"
        n1_guess = n1_guess - 2
        n1 = n1_minimo(Re, U0_top, n1_guess, dt_guess, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros, e_tol, u_central_guia)
    else
        % Caso negativo: aumento en 2 n1_guess y ya tengo el n1
        n1 = n1_guess + 2
        return
    endif

endfunction


#TEST:

% Re = 1000;
% termino_advectivo = "D";
% tol_estacionario = 1e-3;
% u_central_guia = -0.06080;

% e_tol = 0.05;
% n1_guess = 41
% dt_guess = 0.05; #0.001 Con este valor converge para DC2 y Re = 1000
% Ndeltat=8000;

% nsimpler = 1;
% metodo_temporal = "EI";
% archivo_velocidades_centrales = "graficos/datos/velocentral_basura_e.csv";
% archivo_evolucion_variables = "graficos/datos/evolucion_basura_e.csv";
% archivo_parametros = "graficos/datos/parametros_basura_e.csv";

% n1 = n1_minimo(Re, @U0_top_cte, n1_guess, dt_guess, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros, e_tol, u_central_guia)


calcular = false;
if calcular == true
    "Inciso e"
    Re_array = [1000];
    termino_advectivo_array = ['U'] #["D","U"];
    tol_estacionario = 1e-5;
    u_central_guia_array = [-0.20581,-0.06080];

    e_tol = 0.05;
    n1_guess = 80; #Para 100 U es menor a 60. Con 1000 tiene un error de 0.7 para 60
    dt_guess = 0.005#0.005; #0.001 Con este valor converge para DC2 y Re = 1000
    Ndeltat= 80000;

    nsimpler = 1;
    metodo_temporal = "E";
    archivo_velocidades_centrales = "graficos/datos/e_velocentral_basura_e.csv";
    archivo_evolucion_variables = "graficos/datos/e_evolucion_basura_e.csv";
    archivo_parametros = "graficos/datos/e_parametros_basura_e.csv";

    for i=1:length(Re_array)
        Re = Re_array(i)
        u_central_guia = u_central_guia_array(i)
        for j=1:length(termino_advectivo_array)

            #Inicio cronómetro
            t_crono_ini = time();

            #Calculo n1 mínimo
            termino_advectivo = termino_advectivo_array(j)
            n1 = n1_minimo(Re, @U0_top_cte, n1_guess, dt_guess, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros, e_tol, u_central_guia)

            #Detengo cronómetro
            costo_rta = time() - t_crono_ini;


            #Me interesa guardar
            #*termino advectivo
            #*Re
            #*n1
            #*deltat
            #*Ndeltat
            #*costo_rta

            data = importdata(archivo_parametros);
            Ndeltat_rta = data(3);
            dt_rta = data(4);

            #Guardo en un archivo

            archivo_costo_comp = strcat("graficos/datos/e_costo_computacional_Re", num2str(Re),"_termino_adv", termino_advectivo, ".csv");
            archivo_costo_comp = fopen(archivo_costo_comp, "w");
            fprintf(archivo_costo_comp,"%e %e %e %e %e %e\n", termino_advectivo, Re, n1, dt_rta, Ndeltat_rta, costo_rta);

            
        endfor
    endfor

endif







#Inciso f: análogo al anterior pero variando dt y lsimpler
#YA ESTÁ IMPLEMENTADO. HAY QUE CORRERLO SOLO
calcular = false;

if calcular == true
    inciso = "f"
    
    tol_estacionario = 1e-5;
    nsimpler_array = [1,2,3];


    Re = [1000];
    n1 = [30];
    termino_advectivo = "D";
    metodo_temporal = "E";

    buscar = true;
    dt_guess = 0.005;
    Ndeltat_max = 10;

    Ndeltat = 80000;    

    for i=1:length(nsimpler_array)

        #Inicio cronómetro
        t_crono_ini = time();


        #Busco el mejor dt y evoluciono para asegurar convergencia
        nsimpler = nsimpler_array(i)

        archivo_ucentral = strcat("graficos/datos/f_ucentral_nsimpler", num2str(nsimpler), ".csv");
        archivo_vcentral = strcat("graficos/datos/f_vcentral_nsimpler", num2str(nsimpler), ".csv");


        comparacion_Ghuia(Re, @U0_top_cte, n1, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_ucentral, archivo_vcentral, dt_guess, Ndeltat_max, buscar, inciso)



        #Detengo cronómetro
        costo_rta = time() - t_crono_ini;

        #Guardo en un archivo las velocidades a tiempo final
        #Abro el archivo de velocidad
        archivo_velocidades_centrales = strcat("graficos/datos/", inciso, "_",termino_advectivo, "_",num2str(n1), "_",num2str(Re),"_velocentral.csv");
        data = importdata(archivo_velocidades_centrales);
        #Guardo los datos agregando la terminación con el nsimpler
        
        archivo_velocidades_centrales = strcat("graficos/datos/", inciso, "_",termino_advectivo, "_",num2str(n1), "_",num2str(Re), "_lsimpler",num2str(nsimpler),"_velocentral.csv");
        csvwrite(archivo_velocidades_centrales, data);


        #Me interesa guardar el tiempo de ejecución de cada caso y el paso de tiempo
        #Abro el archivo de parámetros y extraigo el costo y el paso de tiempo
        archivo_parametros = strcat("graficos/datos/", inciso, "_",termino_advectivo, "_",num2str(n1(1)), "_",num2str(Re(1)), "_parametros.csv")
        data = importdata(archivo_parametros);
        Ndeltat_rta = data(3);
        dt_rta = data(4);

        #Guardo la info en un archivo distinto para cada lsimpler
        archivo_costo_comp = strcat("graficos/datos/f_costo_computacional_nsimpler", num2str(nsimpler), ".csv");
        archivo_costo_comp = fopen(archivo_costo_comp, "w");
        fprintf(archivo_costo_comp,"%e %e %e\n", Ndeltat_rta, dt_rta, costo_rta);

    endfor

endif


#Inciso g: cambiar condición de contorno
calcular = true;

if calcular == true
    inciso = "g"
    tol_estacionario = 1e-5; #Esta condición no se debería alcanzar. Se debería detener el programa por haber llegado a Nmax.

    function U0_top = U0_top_cos(t)
        U0_top = cos(t);
    endfunction

    #Testeo
    termino_advectivo = "D";
    n1 = 30;
    Re = 1000;
    dt = 0.4;
    tmax = 30;#60;
    Ndeltat = round(tmax/dt)
    nsimpler_array = [3];#[1,3];
    metodo_temporal_array = ['C']#['E'; 'C'];


    for i=1:length(nsimpler_array)
        for j=1:length(metodo_temporal_array)

            nsimpler = nsimpler_array(i)
            metodo_temporal = metodo_temporal_array(j)


            archivo_velocidades_centrales = strcat("graficos/datos/g_velocentral_nsimpler", num2str(nsimpler), "_metodotemporal", metodo_temporal, ".csv");
            archivo_evolucion_variables = strcat("graficos/datos/g_evolucion_nsimpler", num2str(nsimpler), "_metodotemporal", metodo_temporal, ".csv");
            archivo_parametros = strcat("graficos/datos/g_parametros_nsimpler", num2str(nsimpler), "_metodotemporal", metodo_temporal, ".csv");


            #Evoluciono
            Chehade_NSdc2(Re, @U0_top_cos, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)


        endfor
    endfor

endif