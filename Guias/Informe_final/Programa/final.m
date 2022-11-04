% Los parámetros a modificar en el tp final son

% * nombre de los archivos salientes. VERIFICAR si se reportan los valores de u y v para todo tiempo o solo al final, los necesito para todo tiempo
% 


% %Testeo la función
% Re = 100.0;
% n1 = 20;
% dt = 0.5;
% nsimpler = 1;
% termino_advectivo = "DC2";
% metodo_temporal = "EI";
% tol_estacionario = 1e-6;
% Ndeltat=80; #numero de pasos de tiempo
% archivo_velocidades_centrales = "graficos/datos/velocentral.txt";
% archivo_evolucion_variables = "graficos/datos/evolucion.txt";
% archivo_parametros = "graficos/datos/parametros.txt";

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



#Tuve que eliminar el clear del inicio del código. No sé en qué afectará
#Terminar de programar QUICK

% Dudas:
% 1. Identifiqué bien cómo hacer solo un UP1, un DC2 y un QUICK?
% 2. En QUICK solo tengo que modificar el loop que dice "#Interior desde i=2 hasta n1-1 y desde j=2 hasta n1:", no?
% 3. En ese loop viejo se multiplica por dx*uwest(i,j) para directamente calcular un término en el cálculo de la integral 2 (ver mis notas)? O sea que para implementar QUICK eso no lo tengo que cambiar, no?
% 4. Está bien dónde cambié la velocidad Utop1? Está cerca de la línea 620
% 5. Cuando ejecuto el código aparece el resultado
% % ans = 39.309
% Pero no siempre aparece el mismo valor. El código no es determinista?
% 6. Euler implícito solo se usa en la integral 1, no?
% 7. No entiendo cómo aplicar CN a partir de EI. Por qué tengo que usar la mitad del paso de tiempo?
% * La cuenta para CN la puse cuando terminó la iteración simples, está bien?
% *La matriz b de CN solo afecta a las matrices bu y bv en el PASO 1 SIMPLER, no?
% *El valor de ans (14.796 o 42.115) da distinto si conecto o no la notebook a lacorriente. Es el tiempo que le toma al programa hacer la cuenta
% *Están bien los resultados que encontré?
% *Puedo usar números impares para n1 de modo de no tener que interpolar el valor en 0.5?
% *En el inciso c me diverge usando UP1, puede ser que lo haya implementado mal? O será el valor de Deltat?


#En algún momento tengo que sumar su efecto en la ecuación diferencial.

1; #Esto está puesto para que octave no interprete que este es un archivo de funciones

function U0_top = U0_top_cte(x)
    U0_top = 1.0;
endfunction

function comparacion_Ghuia(Re_array, U0_top, n1_array, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_ucentral, archivo_vcentral)
    #Esta función sirve para comparar los resultados de la función Chehade_NSdc2 con los resultados de la función Chehade_NSdc2_Ghuia para los parámetros de entrada de Chehade_NSdc2
    #Re_array: array con los valores de Reynolds a probar
    #U0_top: función que devuelve la velocidad en la parte superior del dominio
    #n1_array: array con los valores de n1 a probar
    #...
    #archivo_ucentral: archivo donde se guardan los valores de u en el centro del dominio para cada valor de Re y n1
    #archivo_vcentral: archivo donde se guardan los valores de v en el centro del dominio para cada valor de Re y n1

    #Creo las matrices que guardarán los valores de U y V en el centro para cada Re y n1
    b_ucentral = zeros(length(n1_array), length(Re_array));
    b_vcentral = zeros(length(n1_array), length(Re_array));

    for ii = 1:length(n1_array)
        for jj=1:length(Re_array)
            n1 = n1_array(ii)
            Re = Re_array(jj)

            archivo_velocidades_centrales = "graficos/datos/velocentral_basura.txt";
            archivo_evolucion_variables = "graficos/datos/evolucion_basura.txt";
            archivo_parametros = "graficos/datos/parametros_basura.txt";

            Chehade_NSdc2(Re, @U0_top_cte, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);

            #Abro el archivo de velocidad central y extraigo los valores de U en y = 0.5 y de V en 
            data = importdata(archivo_velocidades_centrales);
            posicion_centro = (length(data(:,1))+1)/2;
            % data(posicion_centro,1)
            b_ucentral(ii,jj) = data(posicion_centro, 2);
            b_vcentral(ii,jj) = data(posicion_centro, 3);
        endfor
    endfor

    #Guardo ambas matrices en un .csv
    csvwrite(archivo_ucentral, b_ucentral);
    csvwrite(archivo_vcentral, b_vcentral);

endfunction




#Inciso b

#MANTENER EN FALSE, YA HICE UNA CORRIDA LARGA QUE TARDÓ HORAS
calcular = false;

if calcular == true
    #Parámetros que puedo llegar a modificar
    dt = 0.05;
    Ndeltat=80000; #numero de pasos de tiempo

    n1_array = [21,41,81];
    Re_array = [100,1000,5000];
    archivo_ucentral = "graficos/datos/b_ucentral.csv";
    archivo_vcentral = "graficos/datos/b_vcentral.csv";

    #Parámetros que seguro no modifique
    nsimpler = 1;
    termino_advectivo = "DC2";
    metodo_temporal = "EI";
    tol_estacionario = 1e-6;

    comparacion_Ghuia(Re_array, @U0_top_cte, n1_array, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_ucentral, archivo_vcentral)
endif

#Inciso c
#Hago lo mismo que el b pero con distinto esquema para el término advectivo.
#AÚN NO PROGRAMÉ EL QUICK así que ejecuto dos veces con el UP1
calcular = false;
if calcular == true
    termino_advectivo_array = ["UP1, UP1"];

    #Parámetros que puedo llegar a modificar
    dt = 0.001;
    Ndeltat=2000; #numero de pasos de tiempo

    n1_array = [21];#[21,41,81];
    Re_array = [1000];#[100,1000,5000];

    #Parámetros que seguro no modifique
    nsimpler = 1;
    metodo_temporal = "EI";
    tol_estacionario = 1e-6;

    for kk = length(termino_advectivo_array)

        termino_advectivo = termino_advectivo_array(kk);
        archivo_ucentral = strcat("graficos/datos/c_ucentral_",termino_advectivo_array(kk),".csv");
        archivo_vcentral = strcat("graficos/datos/c_vcentral_",termino_advectivo_array(kk),".csv");

        comparacion_Ghuia(Re_array, @U0_top_cte, n1_array, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_ucentral, archivo_vcentral)

    endfor
endif

#Inciso d
#FALTA DEFINIR CUÁL ES EL MEJOR ESQUEMA y VER CÓMO CLAVAR UN VALOR EN T = 0.2 PARA NO TENER QUE HACER INTERPOLACIÓN
calcular = true;

if calcular == true

    #Defino el mejor esquema advectivo
    termino_advectivo_mejor_esquema = "DC2";

    #Defino parámetros
    Re_array = [1,1000];
    n1_array = [10,20,40,60,80];
    dt = 0.5;
    Ndeltat=3; #numero máximo de pasos de tiempo

    #Parámetros que seguro no modifique
    nsimpler = 1;
    metodo_temporal = "EI";
    tol_estacionario = 1e-6;

    #Recorro ambos Re y para cada uno de ellos voy calculando los errores
    for ii = 1:length(Re_array)
        Re = Re_array(ii)

        #Defino los vectores en los que guardaré los errores
        error_uparticular = zeros(1,length(n1_array));
        error_vparticular = zeros(1,length(n1_array));

        #Defino los valores solución
        n1_sol = 80-2;
        % Los valores particulares u(0.5, 0.2) y v(0.5,0.2) los voy a guardar en las variables
        % c_usol_particular
        % c_vsol_particular

        #Defino los nombres de los archivos
        archivo_velocidades_centrales = "graficos/datos/velocentral_basura.txt";
        archivo_evolucion_variables = "graficos/datos/evolucion_basura.txt";
        archivo_parametros = "graficos/datos/parametros_basura.txt";

        #Hago la cuenta
        Chehade_NSdc2(Re, @U0_top_cte, n1_sol, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo_mejor_esquema, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);

        #Extraigo el valor de u(0.5, 0.2) y v(0.5,0.2) de la solución

        #Abro el archivo de velocidad central
        data = importdata(archivo_velocidades_centrales)
        posicion_0_2 = (length(data(:,1))+2)/5
        data(posicion_0_2,1)
        c_usol_particular = data(posicion_0_2, 2);
        c_vsol_particular = data(posicion_0_2, 3);
        pause(100)
        #Recorro los distintos n1
        for jj=1:length(n1_array)
            n1 = n1_array(jj)
            #Hago la cuenta
            Chehade_NSdc2(Re, @U0_top_cte, n1, dt, tol_estacionario, Ndeltat, nsimpler, "DC2", metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);

            #Extraigo el valor de u(0.5, 0.2) y v(0.5,0.2) de la aproximación
            data = importdata(archivo_velocidades_centrales)
            c_uaprox_particular = data(posicion_0_2, 2);
            c_uaprox_particular = data(posicion_0_2, 3);

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