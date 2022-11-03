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

    for ii = 1:3
        for jj=1:3
            n1 = n1_array(ii)
            Re = Re_array(jj)

            archivo_velocidades_centrales = "graficos/datos/velocentral_basura.txt";
            archivo_evolucion_variables = "graficos/datos/evolucion_basura.txt";
            archivo_parametros = "graficos/datos/parametros_basura.txt";

            Chehade_NSdc2(Re, @U0_top_cte, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros);

            #Abro el archivo de velocidad central y extraigo los valores de U en y = 0.5 y de V en x = 0.5
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

calcular = true;

if calcular == true
    #Parámetros que puedo llegar a modificar
    dt = 0.5;
    Ndeltat=800; #numero de pasos de tiempo

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
calcular = true;
if calcular == true
    termino_advectivo_array = ["UP1, UP1"];

    #Parámetros que puedo llegar a modificar
    dt = 0.5;
    Ndeltat=80; #numero de pasos de tiempo

    n1_array = [21,41,81];
    Re_array = [100,1000,5000];

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
#FALTA DEFINIR CUÁL ES EL MEJOR ESQUEMA
Re_array = [1,1000]
n1_sol = 80

