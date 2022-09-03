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
h = SIZE_DOMAIN / (Nn+1); #Distancia entre nodos
h2 = h*h;
x=h*(1:Nn); #Array de 1/h hasta Nn/h

#---------------------------------------------
#SOLUCIÓN EXACTA
#---------------------------------------------






#---------------------------------------------
#MÉTODO DE DIFERENCIAS FINITAS CENTRADAS DE 2DO ORDEN
#---------------------------------------------

#Creo la matriz A

DC_S1 = sparse([1:Nn], [1:Nn], (-2/h2 -1)*ones(1,Nn), Nn, Nn, 0);
DC_S2 = sparse([1:Nn-1], [2:Nn], (1/h2)*ones(1,Nn-1), Nn, Nn, 0);
DC_S3 = sparse([2:Nn], [1:Nn-1],(1/h2)*ones(1,Nn-1), Nn, Nn, 0);
A = DC_S1 + DC_S2 + DC_S3;

# Vectores de solution
DC_ysol = zeros(Nn,1);

# Termino independiente  A * ysol = b

function gx = g(x,K)
    suma = 0;
    for jj=1:K
        suma = suma + (1 + (jj*pi)*(jj*pi))*sin(jj*pi*x);
    endfor
    gx = -suma;
endfunction


b = zeros(Nn,1);
b = g(x,K); #x es un vector. Se puede meter vectores como argumentos de funciones!
b(1) = b(1) - UL/h2;
b(1) = b(Nn) - UR/h2;


#Resolucion del problema

DC_ysol=A\b.';  # SIEMPRE USAR EL COMANDO "\" PARA RESOLVER SISTEMAS LINEALES - MAS ECONOMICO
#El operador ".'" calcula la transpuesta. Esto es necesario porque el vector b es un vector fila.

#Ploteo

#plot(xsol,sol)
%plot(x,DC_ysol,".","markersize", 10);






#---------------------------------------------
#MÉTODO DE PADÉ
#---------------------------------------------

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
P_ysol = (B2-A2)\(A2*b.' - c);

#Ploteo
plot(x,DC_ysol,".","markersize", 10);
pause(10)



