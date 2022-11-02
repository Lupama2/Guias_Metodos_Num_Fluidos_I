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



#Tuve que eliminar el clear del inicio del código. No sé en qué afectará
#Terminar de programar QUICK