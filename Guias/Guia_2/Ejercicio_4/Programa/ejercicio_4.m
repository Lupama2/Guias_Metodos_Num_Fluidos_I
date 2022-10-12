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

global Nmax_bis = 10000;
global tol_bis = 0.00001;
global alea = 0; #Este es un nro aleatorio necesario para ejecutar el método de bisección. No debería tener impacto en el resultado final
global a = -10;
global b = 10;

function f_vec =  f_vec(y_vec, t)
    g = 10;
    l = 1;
    f_vec = [y_vec(2), -g/l * sin(y_vec(1))];
endfunction

function y_new1 = Euler_implicito(y_new, y_old, t_old, h)
    # Dado y_n calcula y_{n+1} mediante el método de Euler implícito
    # y_{n+1} = y_n + h * f_vec(y_{n+1},t_{n+1})
    #y_new es el vector [y1_{n+1}, y2_{n+1}]
    #y_old es el vector [y1_n, y2_n]
    #t_old es el tiempo t_n
    #h es el paso de tiempo
    t_new = t_old + h;
    y_new1 = y_old + h*f_vec(y_new,t_new);
endfunction

function y_new1 = CrackNicholson(y_new, y_old, t_old, h)
    # Dado y_n calcula y_{n+1} mediante el método de Crack-Nicholson
    # y_{n+1} = y_n + h/2 * (f_vec(y_n,t_n) + f_vec(y_{n+1},t_{n+1}))
    #y_new es el vector [y1_{n+1}, y2_{n+1}]
    #y_old es el vector [y1_n, y2_n]
    #t_old es el tiempo t_n
    #h es el paso de tiempo
    t_new = t_old + h;
    y_new1 = y_old + h/2*(f_vec(y_old,t_old) + f_vec(y_new,t_new));
endfunction

function y_new1 = RK4(y_old, t_old, h)
    # Dado y_n calcula y_{n+1} mediante el método de Runge-Kutta de orden 4
    #y_new1 es el vector [y1_{n+1}, y2_{n+1}]. No se coloca como input porque se trata de un método explícito
    #y_old es el vector [y1_n, y2_n]
    #t_old es el tiempo t_n
    #h es el paso de tiempo
    t_new = t_old + h;
    k1 = f_vec(y_old,t_old);
    k2 = f_vec(y_old + h/2*k1,t_old + h/2);
    k3 = f_vec(y_old + h/2*k2,t_old + h/2);
    k4 = f_vec(y_old + h*k3,t_old + h);
    y_new1 = y_old + h/6*(k1 + 2*k2 + 2*k3 + k4);
endfunction

function y_new1 = LeapFrog(y_old,y_old_old, t_old,h)
    # Dado y_n calcula y_{n+1} mediante el método de Leap-Frog (método explícito)
    #y_new es el vector [y1_{n+1}, y2_{n+1}]
    #y_old es el vector [y1_n, y2_n]
    #y_old_old es el vector [y1_n-1, y2_n-1]
    #t_old es el tiempo t_n
    #h es el paso de tiempo
    t_new = t_old + h;
    y_new1 = y_old_old + 2*h*f_vec(y_old,t_old);
endfunction



function [tita_array, tita_punto_array] = Metodo_implicito(metodo, y_ini, t_ini, h, n, a, b, Nmax_bis, tol_bis, alea)
    #Aplico el método implícito al problema y retorno tita y tita_punto para todo t discreto entre 0 y h*n
    #metodo es un puntero al método implícito
    #y_ini es el vector [y1_0, y2_0]
    #t_ini es el tiempo t_0
    #h es el paso de tiempo
    #n es el número de pasos de tiempo
    #a y b son los límites de la bisección
    #Nmax_bis es el número máximo de iteraciones de la bisección
    #tol_bis es la tolerancia de la bisección
    #alea es un nro aleatorio necesario para aplicar el método de la bisección. No debería tener impacto alguno

    tita_array = zeros(n,1);
    tita_punto_array = zeros(n,1);
    #Aplico el método n veces
    t_old = t_ini;
    y_ini_copia = y_ini;
    for ii = 1:n
        [tita_array(ii), tita_punto_array(ii)] = biseccion_gral(metodo,y_ini_copia,t_old,h, a, b, Nmax_bis, tol_bis, alea);
        y_ini_copia = [tita_array(ii), tita_punto_array(ii)];
    endfor
endfunction

function [tita_array, tita_punto_array] = Metodo_explicito(metodo, y_ini, t_ini, h, n)
    #Aplico un método explícito al problema y retorno tita y tita_punto para todo t discreto entre 0 y h*n
    #metodo es un puntero al metodo explícito
    #y_ini es el vector [y1_0, y2_0]
    #t_ini es el tiempo t_0
    #h es el paso de tiempo
    #n es el número de pasos de tiempo

    tita_array = zeros(n,1);
    tita_punto_array = zeros(n,1);
    #Aplico el método n veces
    t_old = t_ini;
    y_ini_copia = y_ini;
    for ii = 1:n
        y_new = metodo(y_ini_copia,t_ini,h);
        tita_array(ii) = y_new(1);
        tita_punto_array(ii) = y_new(2);
        #No se puede hacer directamente [tita_array(ii), tita_punto_array(ii)] = metodo(y_ini_copia,t_ini,h); porque la función metodo registra un único output, no varios.
        y_ini_copia = [tita_array(ii), tita_punto_array(ii)];
    endfor

endfunction

function [tita_array, tita_punto_array] = CN_LeapFrog(y_ini, t_ini, h, n, a, b, Nmax_bis, tol_bis, alea)
    #2 pasos con C-N y uno con Leap-Frog, es decir, se avanza de a 3h y por lo tanto n debe ser múltiplo de 3

    #Controlo que n sea múltiplo de 3
    if(mod(n,3)!=0)
        error("n debe ser multiplo de 3")
    endif

    tita_array = zeros(n,1);
    tita_punto_array = zeros(n,1);

    y_ini_copia = y_ini;
    t_old = t_ini;

    for jj = 1:(n/3)
        ii = 3*(jj-1) + 1;
    #C-N:
        [tita_array(ii), tita_punto_array(ii)] = Metodo_implicito(@CrackNicholson, y_ini_copia, t_old, h, 1, a, b, Nmax_bis, tol_bis, alea);
        y_ini_copia = [tita_array(ii), tita_punto_array(ii)];
        t_old = t_old + h;
        [tita_array(ii+1), tita_punto_array(ii+1)] = Metodo_implicito(@CrackNicholson, y_ini_copia, t_old, h, 1, a, b, Nmax_bis, tol_bis, alea);
        y_ini_copia = [tita_array(ii+1), tita_punto_array(ii+1)];
        t_old = t_old + h;

        #Leap-Frog:
        y_new= LeapFrog([tita_array(ii+1), tita_punto_array(ii+1)],[tita_array(ii), tita_punto_array(ii)], t_old,h);
        tita_array(ii+2) = y_new(1);
        tita_punto_array(ii+2) = y_new(2);
        y_ini_copia = [tita_array(ii+2), tita_punto_array(ii+2)];
        t_old = t_old + h;
    endfor

endfunction

#Aplico los métodos
y_ini = [pi/2,0]; #condiciones iniciales
t_ini = 0;
n = 99;
h = 0.1;
plotear = false;

if plotear == true
    [tita_array_EI, tita_punto_array_EI] = Metodo_implicito(@Euler_implicito, y_ini, t_ini, h, n, a, b, Nmax_bis, tol_bis, alea);
    [tita_array_CN, tita_punto_array_CN] = Metodo_implicito(@CrackNicholson, y_ini, t_ini, h, n, a, b, Nmax_bis, tol_bis, alea);
    [tita_array_RK4, tita_punto_array_RK4] = Metodo_explicito(@RK4, y_ini, t_ini, h, n);
    [tita_array_CNLF, tita_punto_array_CNLF] = CN_LeapFrog(y_ini, t_ini, h, n, a, b, Nmax_bis, tol_bis, alea);


    #Ploteo tita_array vs tiempo:
    t_array = linspace(0, n*h, n);
    plot(t_array,tita_array_EI,";Euler Implícito;","linewidth", 2);
    hold on;
    plot(t_array,tita_array_CN,";Crack - Nicholson;","linewidth", 2);
    hold on;
    plot(t_array,tita_array_RK4,";RK4;","linewidth", 2);
    hold on;
    plot(t_array,tita_array_CNLF,";CN + Leap Frog;","linewidth", 2);

    pause(1000)
endif

function error_fase
    #Dado t_array, tita_array y tita_punto_array, calculo el error de fase
endfunction


function erroramplitud
    #Dado t_array, tita_array y tita_punto_array, calculo el error de amplitud
endfunction





















#Inciso b: péndulo doble
#