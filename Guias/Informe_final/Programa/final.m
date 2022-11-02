% Los parámetros a modificar en el tp final son
% * tol_estacionario. Condición de estado estacionario (agregar una función que elija Ndeltat de modo de correr el código hasta que se llegue a una determinada condición)
% * t_final. Condición adicional por las dudas. Estaría bueno imprimir en la terminal qué condición se cumplió para terminar la simulación
% * nombre de los archivos salientes. VERIFICAR si se reportan los valores de u y v para todo tiempo o solo al final, los necesito para todo tiempo
% * metodo_temporal: "EI" o "CN". 


%Testeo la función
Re = 100.0;
n1 = 20;
dt = 0.5;
nsimpler = 1;
termino_advectivo = "DC2";

function U0_top = U0_top(x)
    U0_top = 1.0;
endfunction

Chehade_NSdc2(Re,n1, dt,nsimpler, termino_advectivo, @U0_top)
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


% Duda:
% 1. Identifiqué bien cómo hacer solo un UP1, un DC2 y un QUICK?
% 2. En QUICK solo tengo que modificar el loop que dice "#Interior desde i=2 hasta n1-1 y desde j=2 hasta n1:", no?
% 3. En ese loop viejo se multiplica por dx*uwest(i,j) para directamente calcular un término en el cálculo de la integral 2 (ver mis notas)? O sea que para implementar QUICK eso no lo tengo que cambiar, no?
% 4. Está bien dónde cambié la velocidad Utop1? Está cerca de la línea 620
% 5. Cuando ejecuto el código aparece el resultado
% % ans = 39.309
% Pero no siempre aparece el mismo valor. El código no es determinista?
% 6. Euler implícito solo se usa en la integral 1, no?
% 7. No entiendo cómo aplicar CN a partir de EI. Por qué tengo que usar la mitad del paso de tiempo?