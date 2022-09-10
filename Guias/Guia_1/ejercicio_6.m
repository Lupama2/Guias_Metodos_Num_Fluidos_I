#
# Resuelvo la ecuación 
#
#



#---------------------------------------------
#Datos del problema
#---------------------------------------------
XL = 0.0; #X inicial
XR = 1.0; #X final
SIZE_DOMAIN = (XR-XL); #Tamanho del dominio
UL = 0.; #Condicion Dirichlet izquierda
UR = 1.; #Condicion Dirichlet derecha
#Nn numero de nodos
#K nro de términos considerados en la suma

#---------------------------------------------
#Datos para la discretizacion
#---------------------------------------------


#término extra en la EDO: y'' - y = g(x)
function gx = g(x,K)
    suma = 0;
    for jj=1:K
        suma = suma + (1 + (jj*pi)*(jj*pi))*sin(jj*pi*x);
    endfor
    gx = -suma;
endfunction

#---------------------------------------------
#SOLUCIÓN EXACTA
#---------------------------------------------

#Defino la función que es solución exacta
function y = y(x,K)
    C1 = 1/(e - power(e,-1));
    C2 = -C1;
    y_h = C1*exp(x) + C2*exp(-x);
    y_p = 0;
    for jj=1:K
        y_p = y_p + sin(jj*pi*x);
    endfor
    y = y_h + y_p;
endfunction

#Defino la función error
function err = error(y_num, x, K, norma)
    #Para dado array calcula la diferencia entre la solución exacta y la aproximación
    #y_num: solución numérica
    #x: valores en los que se calcula la solución numérica
    #norma: norma de la diferencia
    err = norm(y_num - y(x,K), norma);
endfunction


#---------------------------------------------
#MÉTODO DE DIFERENCIAS FINITAS CENTRADAS DE 2DO ORDEN
#---------------------------------------------
function [DC_x, DC_ysol] = DC_num_sol(Nn,K,UL,UR,SIZE_DOMAIN)
    # Función que resuelve el problema con diferencias centradas y devuelve los puntos en los que
    # se evaluó la solución y los valores de la solución aproximada en dichos puntos
    # Nn: nro de nodos
    # K: nro de términos considerados en la suma de la EDO
    # UL: Condicion Dirichlet izquierda
    # UR: Condicion Dirichlet derecha
    # x: Array de 1/h hasta Nn/h
    # h: Distancia entre nodos
    # h2: h*h

    #Defino el vector x donde aproximaré la solución
    h = SIZE_DOMAIN / (Nn+1); #Distancia entre nodos
    h2 = h*h;
    x=h*(1:Nn); #Array de 1/h hasta Nn/h

    #Creo la matriz A

    DC_S1 = sparse([1:Nn], [1:Nn], (-2/h2 -1)*ones(1,Nn), Nn, Nn, 0);
    DC_S2 = sparse([1:Nn-1], [2:Nn], (1/h2)*ones(1,Nn-1), Nn, Nn, 0);
    DC_S3 = sparse([2:Nn], [1:Nn-1],(1/h2)*ones(1,Nn-1), Nn, Nn, 0);
    A = DC_S1 + DC_S2 + DC_S3;

    #Vector solución
    DC_ysol = zeros(Nn,1);

    #Termino independiente  A * ysol = b

    #b = zeros(Nn,1);
    b = g(x,K); #x es un vector. Se puede meter vectores como argumentos de funciones!
    b(1) = b(1) - UL/h2;
    b(Nn) = b(Nn) - UR/h2;

    #Resolucion del problema

    DC_ysol=A\b.';  # SIEMPRE USAR EL COMANDO "\" PARA RESOLVER SISTEMAS LINEALES - MAS ECONOMICO
    #El operador ".'" calcula la transpuesta. Esto es necesario porque el vector b es un vector fila.
    DC_x = x;

endfunction

%Para hallar la solución numérica:
%[DC_x, DC_ysol] = DC_num_sol(Nn,K,UL,UR,x,h,h2);
%Ploteo:
%plot(DC_x,DC_ysol,".","markersize", 10);
%pause(10);


#Testeo el funcionamiento de la función anterior. El siguiente código muestra que la solución converge
#para N grande a la solución exacta

% [DC_x, DC_ysol] = DC_num_sol(23,6,UL,UR,SIZE_DOMAIN);
% plot(DC_x,DC_ysol,".","markersize", 10);
% hold on;
% [DC_x, DC_ysol] = DC_num_sol(15,6,UL,UR,SIZE_DOMAIN);
% plot(DC_x,DC_ysol,".","markersize", 10);
% hold on;
% [DC_x, DC_ysol] = DC_num_sol(55,6,UL,UR,SIZE_DOMAIN);
% plot(DC_x,DC_ysol,".","markersize", 10);
% hold on;
% [DC_x, DC_ysol] = DC_num_sol(105,6,UL,UR,SIZE_DOMAIN);
% plot(DC_x,DC_ysol,".","markersize", 10);
% hold on;
% [DC_x, DC_ysol] = DC_num_sol(205,6,UL,UR,SIZE_DOMAIN);
% plot(DC_x,DC_ysol,".","markersize", 10);
% hold on;
% ysol = y(DC_x, 6);
% plot(DC_x,ysol,".","markersize", 10);
% pause(10)


#---------------------------------------------
#MÉTODO DE PADÉ
#---------------------------------------------
function [P_x, P_ysol] = P_num_sol(Nn,K,UL,UR,SIZE_DOMAIN)
    # Función que resuelve el problema con el método de Padé
    # Nn: nro de nodos
    # K: nro de términos considerados en la suma de la EDO
    # UL: Condicion Dirichlet izquierda
    # UR: Condicion Dirichlet derecha
    # x: Array de 1/h hasta Nn/h
    # h: Distancia entre nodos
    # h2: h*h


    #Defino el vector x donde aproximaré la solución
    h = SIZE_DOMAIN / (Nn+1); #Distancia entre nodos
    h2 = h*h;
    x=h*(1:Nn); #Array de 1/h hasta Nn/h

    #Creo la matriz A2
    #Nodos internos:
    P_S1 = sparse([2:Nn-1], [2:Nn-1], (10/12)*ones(1,Nn-2), Nn, Nn, 0);
    P_S2 = sparse([2:Nn-1], [1:Nn-2], (1/12)*ones(1,Nn-2), Nn, Nn, 0);
    P_S3 = sparse([2:Nn-1], [3:Nn],(1/12)*ones(1,Nn-2), Nn, Nn, 0);
    #Nodos de frontera:
    P_F1 = sparse([1], [1], (1), Nn, Nn, 0);
    P_F2 = sparse([Nn], [Nn], (1), Nn, Nn, 0);

    A2 = P_S1 + P_S2 + P_S3 + P_F1 + P_F2;

    #Creo la matriz B2
    #Nodos internos:
    P_S4 = sparse([2:Nn-1], [2:Nn-1], (-2/h2)*ones(1,Nn-2), Nn, Nn, 0);
    P_S5 = sparse([2:Nn-1], [1:Nn-2], (1/h2)*ones(1,Nn-2), Nn, Nn, 0);
    P_S6 = sparse([2:Nn-1], [3:Nn],(1/h2)*ones(1,Nn-2), Nn, Nn, 0);

    #Nodos de frontera:
    P_F3 = sparse([1,1], [1,2], [-2/h2,1/h2], Nn, Nn, 0);
    P_F4 = sparse([Nn,Nn], [Nn-1,Nn], [1/h2,-2/h2], Nn, Nn, 0);

    B2 = P_S4 + P_S5 + P_S6 + P_F3 + P_F4;

    # Vectores de solution
    P_ysol = zeros(Nn,1);

    # Termino independiente  A2*ysol'' = B2*ysol + c
    c = zeros(Nn,1);
    c(1) = UL/h2;
    c(Nn) = UR/h2;

    #Resuelvo el problema
    P_ysol = (B2-A2)\(A2*g(x,K).' - c);
    P_x = x;

endfunction

#Testeo el funcionamiento de la función anterior. El siguiente código muestra que la solución converge
#para N grande a la solución exacta.

% [P_x, P_ysol] = P_num_sol(23,6,UL,UR,SIZE_DOMAIN);
% plot(P_x,P_ysol,".","markersize", 10);
% hold on;
% [P_x, P_ysol] = P_num_sol(15,6,UL,UR,SIZE_DOMAIN);
% plot(P_x,P_ysol,".","markersize", 10);
% hold on;
% [P_x, P_ysol] = P_num_sol(55,6,UL,UR,SIZE_DOMAIN);
% plot(P_x,P_ysol,".","markersize", 10);
% hold on;
% [P_x, P_ysol] = P_num_sol(105,6,UL,UR,SIZE_DOMAIN);
% plot(P_x,P_ysol,".","markersize", 10);
% hold on;
% ysol = y(P_x, 6);
% plot(P_x,ysol,".","markersize", 10);
% pause(10)


#Comparo la solución con ambos métodos para N grande:
% [P_x, P_ysol] = P_num_sol(105,6,UL,UR,SIZE_DOMAIN);
% plot(P_x,P_ysol,".","markersize", 10);
% hold on;
% [DC_x, DC_ysol] = DC_num_sol(205,6,UL,UR,SIZE_DOMAIN);
% plot(DC_x,DC_ysol,".","markersize", 10);
% hold on;
% pause(10)


#---------------------------------------------
#INCISO e
#---------------------------------------------
K_e = 6;
Nn_e = 23;

#Gráfico de solución exacta y ambas aproximaciones:
Nn_sol = 200;
h_sol = SIZE_DOMAIN / (Nn_sol+1); #Distancia entre nodos
h2_sol = h_sol*h_sol;
xsol=h_sol*(1:Nn_sol); #Array de 1/h hasta Nn/h

ysol = y(xsol, K_e);

[DC_x, DC_ysol]  = DC_num_sol(Nn_e,K_e,UL,UR,SIZE_DOMAIN);
[P_x, P_ysol]  = P_num_sol(Nn_e,K_e,UL,UR,SIZE_DOMAIN);

plotear = false;

if(plotear == true)
    plot(xsol,ysol,";Exacta;","linewidth", 2);
    hold on
    plot(DC_x,DC_ysol,".","markersize", 10);
    hold on
    plot(P_x,P_ysol,".","markersize", 10);
    pause(10)
endif

#Grafico el error en el punto central x = 0.5 en función de h
#El punto central es x((Nn+1)/2)

function [DF_errorx0, P_errorx0] = error_en_x(x0, Nn, K, UL, UR, SIZE_DOMAIN)
    #Calcula el error en x = x0 para el método de diferencias finitas y el método de Padé
    #Nn: nro de nodos
    #K: K de la EDO
    #POR LO PRONTO SÓLO FUNCIONA PARA X0 = 0.5


    #Calculo la aproximación de diferencias finitas en x0
    [DC_x, DC_ysol]  = DC_num_sol(Nn, K,UL,UR,SIZE_DOMAIN);
    norma = 2;
    DF_errorx0 = error(DC_ysol((Nn+1)/2), x0, K, norma);

    #Calculo la aproximación de padé en X0
    [P_x, P_ysol]  = P_num_sol(Nn ,K ,UL,UR,SIZE_DOMAIN);  
    norma = 2;
    P_errorx0 = error(P_ysol((Nn+1)/2), x0, K, norma);


endfunction


Nn_array = [5,7,9,11,15,21,41,81,161,321, 641, 1281, 2561, 5121, 10241, 20481,  40961]; %tienen que ser valores impares
%Le tomamucho tiempo cargar:
%Nn_array = [5,7,9,11,15,21,41,81,161,321, 641, 1281, 2561, 5121, 10241, 20481,  40961, 81921, 163841, 327681, 655361, 1310721]; %tienen que ser valores impares
h_array = zeros(length(Nn_array),1).';
DC_error = zeros(length(Nn_array),1).';
P_error = zeros(length(Nn_array),1).';
for ii=1:length(Nn_array)
    [DC_error(ii), P_error(ii)] = error_en_x(0.5, Nn_array(ii), K_e, UL, UR, SIZE_DOMAIN);
    h_array(ii) = SIZE_DOMAIN / (Nn_array(ii)+1); 
endfor

plotear = false;
if (plotear == true)
    plot(h_array,DC_error,".","markersize", 10)
    hold on
    plot(h_array,P_error,".","markersize", 10)
    pause(10)
endif

#Determino el orden de truncamiento
plotear = false;
if (plotear == true)
    #Recta con dependencia de orden 2
    y2 = zeros(length(Nn_array),1).';
    for ii=1:length(Nn_array)
        y2(ii) = h_array(ii)*h_array(ii);
    endfor
    loglog(h_array, y2,";Orden2;","linewidth", 2)
    hold on

    #Recta con dependencia de orden 4
    y4 = zeros(length(Nn_array),1).';
    for ii=1:length(Nn_array)
        y4(ii) = power(h_array(ii),4);
    endfor
    loglog(h_array, y4,";Orden4;","linewidth", 2)
    hold on


    loglog(h_array,DC_error,".","markersize", 10)
    hold on
    loglog(h_array,P_error,".","markersize", 10)
    pause(10)
endif


#---------------------------------------------
#INCISO f
#---------------------------------------------

function [DC_Nnmin, P_Nnmin] = Nm_minimo(tol)
    #Calcula el nro de puntos mínimos Nn necesarios para obtener un error menor a tol.
    #A priori es muy ineficiente: va calculando con Nn sucesivos hasta llegar a un error menor a tol
    

endfunction



