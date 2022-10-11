# Ejercicio 4 - Guía 2 - Métodos numéricos para fluidos I
# Pablo Chehade
# Última edición: 2022/10/11
# Se busca resolver el problema del péndulo simple y péndulo doble

# Inciso a: péndulo simple
# EDO:
# tita''(t) = -g/l * sin(tita(t))
# o, escrito como sistema de ecuaciones
# u' = v
# v' = -g/l * sin(u)
# o, vectorialmente
# dy_vec/dt = f_vec(y_vec,t), donde y_vec = (u,v)



#---------------------------------------------
#Datos del problema
#---------------------------------------------
global g = 10;
global l = 1;

#---------------------------------------------
#Datos de la discretización
#---------------------------------------------
global h = 0.1;
global Nmax_bis = 10000;
global tol_bis = 0.0001;
global alea = 0; #Este es un nro aleatorio necesario para ejecutar el método de bisección. No debería tener impacto en el resultado final
global a = -10;
global b = 10;

function f_vec =  f_vec(y_vec, t)
    g = 10;
    l = 1;
    f_vec = [y_vec(2), -g/l * sin(y_vec(1))];
endfunction

function y_new2 = Euler_implicito(y_new, y_old, t_old, h)
    # Dado y_n calcula y_{n+1} mediante el método de Euler implícito
    # y_{n+1} = y_n + h * f_vec(y_{n+1},t_{n+1})
    #y_new es el vector [y1_{n+1}, y2_{n+1}]
    #y_old es el vector [y1_n, y2_n]
    #t_old es el tiempo t_n
    #h es el paso de tiempo
    t_new = t_old + h;
    y_new2 = y_old + h*f_vec(y_new,t_new);
endfunction

function y_new2 = CrackNicholson(y_new, y_old, t_old, h)

endfunction

function y_new2 = RK4(y_new, y_old, t_old, h)

endfunction

function y_new2 = CN_LeapFrog(y_new, y_old, t_old, h)

endfunction





#Test del método de bisección
#Busco resolver la ec sqrt(x) = 2, sol x = 4
#La planteo como sqrt(x) - 2 = 0 con límites a = 0, b = 8

% function f_ = f(x)
%     f_ = sqrt(x) - 2;
% endfunction
% a = 2;
% b = 8;
% #Parámetros de corte:
% Nmax = 10;
% tol = 0.01;
% biseccion(@f,a,b, Nmax, tol)



#Aplico el método de EI
y_ini = [pi/2,0]; #condiciones iniciales
t_old = 0;

n = 100;
tita_array = zeros(n,1);
tita_punto_array = zeros(n,1);
#Aplico el método n veces
for ii = 1:n
    [tita_array(ii), tita_punto_array(ii)] = biseccion_gral(@Euler_implicito,y_ini,t_old,h, a, b, Nmax_bis, tol_bis, alea);
    y_ini = [tita_array(ii), tita_punto_array(ii)];
endfor

#Ploteo tita_array vs tiempo:
t_array = linspace(0, n*h, n)
tita_array
plot(t_array,tita_array,";Aprox;","linewidth", 2);
pause(1000)























#Inciso b: péndulo doble
#