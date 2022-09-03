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
K = 6; #Nro de términos considerados en la suma

#---------------------------------------------
#Datos para la discretizacion
#---------------------------------------------

Nn = 23; #Numero de nodos



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
function y = y(x)
    y = 0*x; #No sé aún cómo hacer esto
endfunction

#Defino la función error
function err = error(y_num, x, norma)
    #Para dado array calcula la diferencia entre la solución exacta y la aproximación
    #y_num: solución numérica
    #x: valores en los que se calcula la solución numérica
    #norma: norma de la diferencia
    err = norm(y_num - y(x), norma);
endfunction


#---------------------------------------------
#MÉTODO DE DIFERENCIAS FINITAS CENTRADAS DE 2DO ORDEN
#---------------------------------------------
function DC_ysol = DC_num_sol(Nn,K,UL,UR,SIZE_DOMAIN)
    # Función que resuelve el problema con diferencias centradas
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

    # Vectores de solution
    DC_ysol = zeros(Nn,1);

    # Termino independiente  A * ysol = b

    b = zeros(Nn,1);
    b = g(x,K); #x es un vector. Se puede meter vectores como argumentos de funciones!
    b(1) = b(1) - UL/h2;
    b(1) = b(Nn) - UR/h2;

    #Resolucion del problema

    DC_ysol=A\b.';  # SIEMPRE USAR EL COMANDO "\" PARA RESOLVER SISTEMAS LINEALES - MAS ECONOMICO
    #El operador ".'" calcula la transpuesta. Esto es necesario porque el vector b es un vector fila.




endfunction

%Para hallar la solución numérica:
%DC_ysol = DC_num_sol(Nn,K,UL,UR,x,h,h2);
%Ploteo:
%plot(x,DC_ysol,".","markersize", 10);
%pause(10);


#---------------------------------------------
#MÉTODO DE PADÉ
#---------------------------------------------
function P_ysol = P_num_sol(Nn,K,UL,UR,SIZE_DOMAIN)
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
    P_S4 = sparse([2:Nn-1], [2:Nn-1], (-2)*ones(1,Nn-2), Nn, Nn, 0);
    P_S5 = sparse([2:Nn-1], [1:Nn-2], (1)*ones(1,Nn-2), Nn, Nn, 0);
    P_S6 = sparse([2:Nn-1], [3:Nn],(1)*ones(1,Nn-2), Nn, Nn, 0);

    #Nodos de frontera:
    P_F3 = sparse([1,1], [1,2], [-2,1], Nn, Nn, 0);
    P_F4 = sparse([Nn,Nn], [Nn-1,Nn], [1,-2], Nn, Nn, 0);

    B2 = P_S4 + P_S5 + P_S6 + P_F3 + P_F4;

    # Vectores de solution
    P_ysol = zeros(Nn,1);

    # Termino independiente  A2*ysol'' = B2*ysol + c
    c = zeros(Nn,1);
    c(1) = UL/h2;
    c(Nn) = UR/h2;

    #Resuelvo el problema
    P_ysol = (B2-A2)\(A2*g(x,K).' - c);

endfunction

%Para hallar la solución numérica:
%P_ysol = P_num_sol(Nn,K,UL,UR,x,h,h2);
#Ploteo:
%plot(x,P_ysol,".","markersize", 10);
%pause(10)


#---------------------------------------------
#INCISO e
#---------------------------------------------
K = 6;
Nn = 23;

#Gráfico de solución exacta y ambas aproximaciones:

Nnsol = 100; #Numero de nodos
hsol = SIZE_DOMAIN / (Nnsol+1); #Distancia entre nodos
xsol =hsol*(1:Nnsol); #Array de 1/h hasta Nn/h
ysol = y(xsol);

Nn = 81;

DC_ysol = DC_num_sol(Nn,K,UL,UR,SIZE_DOMAIN);
P_ysol = P_num_sol(Nn,K,UL,UR,SIZE_DOMAIN);

plotear = true;

if(plotear == true)

    h = SIZE_DOMAIN / (Nn+1);
    x=h*(1:Nn);

    plot(xsol,ysol,";Exacta;","linewidth", 2);
    hold on
    plot(x,DC_ysol,".","markersize", 10,x,DC_ysol,";DC;","linewidth", 1);
    hold on
    plot(x,P_ysol,".","markersize", 10,x,P_ysol,";P;","linewidth", 1);
    pause(10)
endif

#Grafico el error en el punto central x = 0.5 en función de h
#El punto central es x((Nn+1)/2)


Nn_array = [5,7,9,11,15,21,41,81,161,321]; %tienen que ser valores impares
h_array = zeros(length(Nn_array),1).';
DC_error = zeros(length(Nn_array),1).';
P_error = zeros(length(Nn_array),1).';
for ii=1:length(Nn_array)
    DC_ysol = DC_num_sol(Nn_array(ii),K,UL,UR,SIZE_DOMAIN);
    P_ysol = P_num_sol(Nn_array(ii),K,UL,UR,SIZE_DOMAIN);
    DC_error(ii) = y(0.5) - DC_ysol((Nn_array(ii)+1)/2);
    P_error(ii) = y(0.5) - P_ysol((Nn_array(ii)+1)/2);

    h_array(ii) = SIZE_DOMAIN / (Nn_array(ii)+1); 
endfor


%plot(h_array,DC_error,".");
%hold on
plot(h_array,P_error,".");
pause(10)


