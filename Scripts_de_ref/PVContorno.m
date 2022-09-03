#clear
#
# Resolvemos y''(x)+4*y(x)+3*y(x)=0 en [0,1] Utilizando Diferencias Centradas 2 - Condiciones Borde Dirichlet en ambos extremos
#
#---------------------------------------------
#Datos del problema
#---------------------------------------------

XL = 0.0; #X inicial
XR = 1.0; #X final
SIZE_DOMAIN = (XR-XL); #Tamanho del dominio
UL = 0.; #Condicion Dirichlet izquierda
UR = 1.; #Condicion Dirichlet derecha

#---------------------------------------------
#Datos para la discretizacion
#---------------------------------------------

nn = 26; #Numero de nodos
dx = SIZE_DOMAIN / (nn+1); #Distancia entre nodos
dx2 = dx*dx;
x=dx*(1:nn);

#---------------------------------------------
#Reservar Vectores y Matrices
#---------------------------------------------

# Vamos a trabajar con los nodos internos y pasar al lado derecho las condiciones
# de contorno => habra (nn) incognitas.
# Como es un problema lineal las matrices no dependen de la solucion
# Todo multiplicado por dx2

nunk = nn;

# SIEMPRE USAR MATRICES SPARSE O RALAS PARA NO OPERAR CON CEROS

S1 = sparse([1:nunk], [1:nunk], (-2+3*dx2)*ones(1,nunk), nunk, nunk, 0);
S2 = sparse([1:nunk-1], [2:nunk], (2*dx+1)*ones(1,nunk-1), nunk, nunk, 0);
S3 = sparse([2:nunk], [1:nunk-1], (-2*dx+1)*ones(1,nunk-1), nunk, nunk, 0);
# Sparce Matrix
A = S1+S2+S3;

# Vectores de solution
usol = zeros(nn,1);

# Termino independiente  A * usol = C
for ni=1:nn
  C(ni,1) = 0;
endfor
C(1,1) = -(-2*dx+1)*UL;
C(nn,1) = -(2*dx+1)*UR;

#---------------------------------------------
#Solucion (nodal) exacta
#---------------------------------------------

for ni=1:nn
  solex(ni,1) =(UL*e^2-UR*e^3)/(e^2-1)*exp(-3*ni*dx)+(-UL+UR*e^3)/(e^2-1)*exp(-ni*dx);
endfor

#---------------------------------------------
#Resolucion del problema
#---------------------------------------------

usol=A\C;  # SIEMPRE USAR EL COMANDO "\" PARA RESOLVER SISTEMAS LINEALES - MAS ECONOMICO
#calculo del error

error2=norm(usol-solex);
errorinf=norm(usol-solex,inf);

#---------------------------------------------
#Solucion  exacta con muchos puntos
#---------------------------------------------
xsol=1/1001*(1:1000);
dxsol=1/1001;
for ni=1:1000
  sol(ni,1) = (UL*e^2-UR*e^3)/(e^2-1)*exp(-3*ni*dxsol)+(-UL+UR*e^3)/(e^2-1)*exp(-ni*dxsol);;
endfor

#---------------------------------------------
#Ploteo
#---------------------------------------------

#plot(xsol,sol)
plot(x,usol,".","markersize", 10,xsol,sol,";Exacta;","linewidth", 2)

pause(10);