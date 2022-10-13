#clear all

#  Ecuacion de tercer orden con solucion real
#  y'''-alfa*y''+w^2*y'-alfa*w^2*y=0
#
# Exponentes de las 3 ecuaciones modelos (alfa, iw conjugados)
alfa=-0.28;
w=1.0; 

# Condición inicial a partir de A, B y C, A*exp(alfa*t)+B*exp(wi*t)+C*exp(-wi*t)
A=.1;
B=.1;
C=B;  ##### no modificar
u0=[A+2*B A*alfa A*alfa^2-2*B*w^2];


# Linea temporal
Tini=0;    # Tiempo inicial
Tfin=400;   # Tiempo final
n=401;     # Instantes de tiempo donde se evaluará la solución
t=linspace(Tini, Tfin, n)';

# Solucion Exacta 
nex=1500;
tex=linspace(Tini, Tfin, nex);
for i=1:nex
 uex(i)=A*exp(alfa*tex(i))+2*B*cos(w*tex(i));
endfor


# Función para calcular las derivadas de las funciones incógnita,
# a partir de un estado de las variables y del tiempo
function xdot = myf (x, t, alfa, w)
  xdot(1) = x(2);
  xdot(2) = x(3);
  xdot(3) = alfa*x(3)-w^2*x(2)+alfa*w^2*x(1);
endfunction

# Función que calcula la solución completa,
# a partir de la función derivada y las condiciones iniciales.
# El vector t indica los puntos en los que se calculará la solución
# func: string con el nombre de la función
# u0: vector (fila) con las condiciones iniciales
# t: vector (columna) con los tiempos

#########################################################
#                  Euler Explicito
##########################################################
function u = EE (func, u0, t, alfa, w)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial

# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
    delta1 = feval(func,u(i-1,:)', t(i-1), alfa, w);
    u(i,:) = u(i-1,:) + h * delta1;

  endfor
endfunction

#########################################################
#                  Euler Implicito
##########################################################
function u = EI (func, u0, t, alfa, w)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial



# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
# Definimos la matriz implicita
  A=eye(m);
  B=[1 -h 0; 0 1 -h; -alfa*w^2*h w^2*h 1-h*alfa];
#   Metodo iterativo sistemas no lineales
#    j=1;
#    uaprox=u(i-1,:);
#    errores=1;
#    while (j<50 & errores > 10^-6)  
#      delta1 = feval(func,uaprox', t(i), alfa, w);
#      u(i,:) = u(i-1,:) + h * delta1;
#      j=j+1;
#      errores=norm(u(i,:)-uaprox);
#      uaprox=u(i,:);
#    endwhile
    u(i,:)=(B\A*u(i-1,:)')';

  endfor
endfunction

#########################################################
#                  C-N
##########################################################
function u = CN (func, u0, t, alfa, w)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial

# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
#    Caso no lineal
#    j=1;
#    uaprox=u(i-1,:);
#    errores=1;
#    while (j<50 & errores > 10^-6)  
#      delta1 = feval(func,uaprox', t(i), alfa, w);
#      u(i,:) = u(i-1,:) + h/2 *(delta1+feval(func,u(i-1,:)', t(i-1), alfa, w));
#      j=j+1;
#      errores=norm(u(i,:)-uaprox);
#      uaprox=u(i,:);
#    endwhile

# Definimos la matriz implicita
    A=[1 h/2 0; 0 1 h/2; alfa*w^2*h/2 -w^2*h/2 1+h/2*alfa];
    B=[1 -h/2 0; 0 1 -h/2; -alfa*w^2*h/2 w^2*h/2 1-h/2*alfa];
    u(i,:)=(B\A*u(i-1,:)')';

  endfor
endfunction

#########################################################
#                  RK2
##########################################################
function u = RK2 (func, u0, t, alfa, w)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial

# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
    k1 = h*feval(func,u(i-1,:)', t(i-1), alfa, w);
    k2 = h*feval(func,u(i-1,:)'+k1'/2, t(i-1)+h/2, alfa, w);
    u(i,:) = u(i-1,:) + k2;
  
  endfor
endfunction

#########################################################
#                  RK4
##########################################################
function u = RK4 (func, u0, t, alfa, w)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial

# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
    k1 = h*feval(func,u(i-1,:)', t(i-1), alfa, w);
    k2 = h*feval(func,u(i-1,:)'+k1'/2, t(i-1)+h/2, alfa, w);
    k3 = h*feval(func,u(i-1,:)'+k2'/2, t(i-1)+h/2, alfa, w);
    k4 = h*feval(func,u(i-1,:)'+k3', t(i-1)+h, alfa, w);
    u(i,:) = u(i-1,:) + (k1+k4)/6 + (k2+k3)/3;
  
  endfor
endfunction

#########################################################
#                  AB-2
##########################################################
function u = AB2 (func, u0, t, alfa, w)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial
# Usamos RK4 para calcular el punto faltante
  i=2;
  h = t(i)-t(i-1);
  k1 = h*feval(func,u(i-1,:)', t(i-1), alfa, w);
  k2 = h*feval(func,u(i-1,:)'+k1'/2, t(i-1)+h/2, alfa, w);
  k3 = h*feval(func,u(i-1,:)'+k2'/2, t(i-1)+h/2, alfa, w);
  k4 = h*feval(func,u(i-1,:)'+k3', t(i-1)+h, alfa, w);
  u(i,:) = u(i-1,:) + (k1+k4)/6 + (k2+k3)/3;
#
# Recorremos los siguientes puntos
  for i = 3:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
    u(i,:) = u(i-1,:) + 3*h/2*feval(func,u(i-1,:)', t(i-1), alfa, w)-h/2*feval(func,u(i-2,:)', t(i-2), alfa, w);
  
  endfor
endfunction

# ####################  MAIN ##################################################################

#
# Llamamos al integrador 
  u=RK2("myf",u0, t, alfa, w);

plot(t,u(:,1),'linestyle','-','or','color','r','linewidth',1,tex,uex,'linewidth',2);
#hold on
legend("RK2","Exacta")
title("h=0.125")
axis([250 400  -0.2 0.4])
