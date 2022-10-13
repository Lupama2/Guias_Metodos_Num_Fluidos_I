clear all

#  Ecuacion de Bernoulli no lineal ODE
#  y'''-alfa*y''+w^2*y'-alfa*w^2*y=0
#

# Condición inicial 
Cini=-.9;
u0=[(Cini+1)^(2/3)];

# Linea temporal
Tini=0;    # Tiempo inicial
Tfin=30;   # Tiempo final
n=66;     # Instantes de tiempo donde se evaluará la solución
t=linspace(Tini, Tfin, n)';

# Solucion Exacta 
tex=linspace(Tini, Tfin, 30*n);
for i=1:30*n
 uex(i)=(Cini*exp(-3/2*tex(i))+1)^(2/3);
endfor


# Función para calcular las derivadas de las funciones incógnita,
# a partir de un estado de las variables y del tiempo
function xdot = myf (x, t)
#  xdot(1) = x(1)+x(1)^2*exp(t);   #           2/3*t*x(1)*(x(1)^3-1)/(1+t^2);
  xdot(1) = -x(1)+x(1)^(-1/2);
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
function u = EE (func, u0, t)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial

# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
    delta1 = feval(func,u(i-1,:)', t(i-1));
    u(i,:) = u(i-1,:) + h * delta1;

  endfor
endfunction

#########################################################
#                  Euler Implicito
##########################################################
function u = EI (func, u0, t)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial



# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
# Definimos la matriz implicita

#   Metodo iterativo sistemas no lineales
    j=1;
    uaprox=u(i-1,:)+h*feval(func,u(i-1,:)', t(i));
    errores=1;
    while (j<1050 & errores > 10^-12)  
      delta1 = feval(func,uaprox', t(i));
      u(i,:) = u(i-1,:) + h * delta1;
      j=j+1;
      errores=norm(u(i,:)-uaprox);
      uaprox=u(i,:);
    endwhile

  endfor
endfunction

#########################################################
#                  C-N
##########################################################
function u = CN (func, u0, t)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial

# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
#    Caso no lineal
    j=1;
    uaprox=u(i-1,:);
    errores=1;
    while (j<1000 & errores > 10^-12)  
      delta1 = feval(func,uaprox', t(i));
      u(i,:) = u(i-1,:) + h/2 *(delta1+feval(func,u(i-1,:)', t(i-1)));
      j=j+1;
      errores=norm(u(i,:)-uaprox);
      uaprox=u(i,:);
    endwhile

# Definimos la matriz implicita
#    A=[1 h/2 0; 0 1 h/2; alfa*w^2*h/2 -w^2*h/2 1+h/2*alfa];
#    B=[1 -h/2 0; 0 1 -h/2; -alfa*w^2*h/2 w^2*h/2 1-h/2*alfa];
#    u(i,:)=(B\A*u(i-1,:)')';

  endfor
endfunction

#########################################################
#                  RK2
##########################################################
function u = RK2 (func, u0, t)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial

# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
    k1 = h*feval(func,u(i-1,:)', t(i-1));
    k2 = h*feval(func,u(i-1,:)'+k1'/2, t(i-1)+h/2);
    u(i,:) = u(i-1,:) + k2;
  
  endfor
endfunction

#########################################################
#                  RK4
##########################################################
function u = RK4 (func, u0, t)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial

# Recorremos los siguientes puntos
  for i = 2:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
    k1 = h*feval(func,u(i-1,:)', t(i-1));
    k2 = h*feval(func,u(i-1,:)'+k1'/2, t(i-1)+h/2);
    k3 = h*feval(func,u(i-1,:)'+k2'/2, t(i-1)+h/2);
    k4 = h*feval(func,u(i-1,:)'+k3', t(i-1)+h);
    u(i,:) = u(i-1,:) + (k1+k4)/6 + (k2+k3)/3;
  
  endfor
endfunction

#########################################################
#                  AB-2
##########################################################
function u = AB2 (func, u0, t)
  n = size(t)(1); # Cantidad de puntos en los que se requiere la solución
  m = size(u0)(2);# Cantidad de funciones incógnita
  u=zeros(n,m);   # matriz solución
  u(1,:) = u0;   # Primer punto: condición inicial
# Usamos RK4 para calcular el punto faltante
  i=2;
  h = t(i)-t(i-1);
  k1 = h*feval(func,u(i-1,:)', t(i-1));
  k2 = h*feval(func,u(i-1,:)'+k1'/2, t(i-1)+h/2);
  k3 = h*feval(func,u(i-1,:)'+k2'/2, t(i-1)+h/2);
  k4 = h*feval(func,u(i-1,:)'+k3', t(i-1)+h);
  u(i,:) = u(i-1,:) + (k1+k4)/6 + (k2+k3)/3;
#
# Recorremos los siguientes puntos
  for i = 3:n
# Evaluamos la función func, pasándole como argumentos el punto anterior y el tiempo anterior
    h = t(i)-t(i-1);
    u(i,:) = u(i-1,:) + 3*h/2*feval(func,u(i-1,:)', t(i-1))-h/2*feval(func,u(i-2,:)', t(i-2));
  
  endfor
endfunction



#
# Llamamos al integrador 
u=EI("myf",u0, t);
plot(t,u(:,1),"linestyle",'-','or','color','r','linewidth',3,tex,uex,'linewidth',1);
#loglog(inter,erroree,inter,errorei,inter,errorcn,inter,y1,'linewidth',2,inter,y2,'linewidth', 2)
legend("EI", "Exacto")
#axis([0 5 0 10])
title("h=0.46154")
