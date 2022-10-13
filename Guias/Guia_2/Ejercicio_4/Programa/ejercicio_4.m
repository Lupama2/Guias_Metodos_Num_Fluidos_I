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
        t_old = t_old + h;
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
        y_new = metodo(y_ini_copia,t_old,h);
        tita_array(ii) = y_new(1);
        tita_punto_array(ii) = y_new(2);
        #No se puede hacer directamente [tita_array(ii), tita_punto_array(ii)] = metodo(y_ini_copia,t_ini,h); porque la función metodo registra un único output, no varios.
        y_ini_copia = [tita_array(ii), tita_punto_array(ii)];
        t_old = t_old + h;
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

function Tau = Periodo_exacto(l,g,tita,tol_periodo)
    #Calcula el período del péndulo simple con la fórmula exacta para dado l y g. Para esto se emplea una suma infinita truncada hasta cierto término. El truncamiento está dado por la tolerancia aceptada tol_periodo, donde el error se calcula como la diferencia entre considerar o no el siguiente término en la suma.
    #l es la longitud del péndulo
    #g es la aceleración de la gravedad
    #tol_periodo es la tolerancia de truncamiento, da la cantidad de términos considerados en la suma infinita
    #tita es la amplitud angular

    T0 = 2*pi*sqrt(l/g);
    Tau = 0;
    error_trunc = tol_periodo + 1;
    nn = 0;
    while(error_trunc > tol_periodo)
        #Calculo el siguiente término de la suma (adimensional)
        termino_suma = ( factorial(2*nn) / (2^(2*nn)*(factorial(nn))^2) )^2 * sin(tita/2)^(2*nn);
        #Actualizo el valor de nn
        nn = nn + 1;
        #Calculo el error de truncamiento
        error_trunc = T0*termino_suma;
        #Actualizo el valor de Tau
        Tau = Tau + error_trunc;
    endwhile

endfunction

% function e_fase = error_fase(t_array, tita_array,tita_punto_array, Tau)
%     #Dado t_array, tita_array y tita_punto_array, calculo el error de fase (tan−1(tita′/tita))
%     #a t = Tau período exacto, fase_exacta = 0, entonces el error estará dado por la fase

%     #tol_fase es el error en la determinación del tiempo tal que t = Tau. Si es muy pequeño podría no encontrarse el tiempo que cumpla la condición.

%     #Busco el tiempo t en t_array tal que t = Tau (no considero múltiplos)
%     % e_fase = 10;
%     for ii=1:length(t_array)
%         #Calculo la fase a ese tiempo
%         if(abs(t_array(ii) - Tau) < 0.001) #considero que el error en la determinación de este nro está dado solamente por un error de punto flotante
%             e_fase = atan(tita_punto_array(ii)/tita_array(ii));
%         endif
%     endfor

% endfunction


function e_fase = error_fase(t_array, tita_array,tita_punto_array, Tau)
    e_fase = atan(tita_punto_array(length(tita_punto_array))/tita_array(length(tita_punto_array)));
endfunction

function ampli = amplitud(tita_punto, tita, l, g)
    ampli = 1/2 * l^2 * tita_punto ^2 - g*l*cos(tita);
endfunction

% function e_ampli = erroramplitud(t_array, tita_array,tita_punto_array, Tau, l, g)
%     #Dado t_array, tita_array y tita_punto_array, calculo el error de amplitud (1/2l^2tita′^2 − gl cos(tita)) en t = Tau período exacto

%     amplitud_exacta = amplitud(tita_punto_array(1), tita_array(1), l, g);
%     #Busco el tiempo t en t_array tal que t = Tau (no considero múltiplos)
%     e_ampli = 10;
%     for ii=1:length(t_array)
%         #Calculo el error de amplitud
%         if(abs(t_array(ii) - Tau)< 0.001)
%             e_ampli = amplitud_exacta - amplitud(tita_punto_array(ii), tita_array(ii), l, g);
%         endif
%     endfor
% endfunction

function e_ampli = erroramplitud(t_array, tita_array,tita_punto_array, Tau, l, g)
    amplitud_exacta = amplitud(tita_punto_array(1), tita_array(1), l, g);
    e_ampli = amplitud_exacta - amplitud(tita_punto_array(length(tita_punto_array)), tita_array(length(tita_punto_array)), l, g);
endfunction

#Aplico los métodos y calculo para cada uno de ellos el error de fase y el error de amplitud
plotear = true;

if plotear == true

    y_ini = [pi/2,0]; #condiciones iniciales
    t_ini = 0;

    #Calculo el período
    tol_periodo = 0.00000001;
    Tau = Periodo_exacto(l,g,y_ini(1),tol_periodo)

    #Elijo n múltiplo de 3
    n_array = 6*[1,2,5,10,20,50,100,200,500,1000,2000]; #debe ser impar para que Tau esté contenido en t_array
    #n_array = 33*[1,3,9,27,81,243, 729]; #debe ser impar para que Tau esté contenido en t_array
    #Elijo h tal que el período esté contenido en t_array una cantidad k de veces
    k = 1;
    h_array = zeros(length(n_array),1);
    for ii = 1:length(n_array)
        h_array(ii) = k*Tau/n_array(ii);
    endfor

    e_fase = zeros(length(n_array),4); #Las filas serán los distintos métodos en el orden EI, CN, RK4, CNFL
    e_amplitud = zeros(length(n_array),4);

    for ii=1:length(n_array)
        t_array = linspace(0, n_array(ii)*h_array(ii), n_array(ii));

        [tita_array_EI, tita_punto_array_EI] = Metodo_implicito(@Euler_implicito, y_ini, t_ini, h_array(ii), n_array(ii), a, b, Nmax_bis, tol_bis, alea);
        e_fase(ii,1) = error_fase(t_array, tita_array_EI,tita_punto_array_EI, Tau);
        e_amplitud(ii,1) = erroramplitud(t_array, tita_array_EI,tita_punto_array_EI, Tau, l, g);

        [tita_array_CN, tita_punto_array_CN] = Metodo_implicito(@CrackNicholson, y_ini, t_ini, h_array(ii), n_array(ii), a, b, Nmax_bis, tol_bis, alea);
        e_fase(ii,2) = error_fase(t_array, tita_array_CN,tita_punto_array_CN, Tau);
        e_amplitud(ii,2) = erroramplitud(t_array, tita_array_CN,tita_punto_array_CN, Tau, l, g);

        [tita_array_RK4, tita_punto_array_RK4] = Metodo_explicito(@RK4, y_ini, t_ini, h_array(ii), n_array(ii));
        e_fase(ii,3) = error_fase(t_array, tita_array_RK4,tita_punto_array_RK4, Tau);
        e_amplitud(ii,3) = erroramplitud(t_array, tita_array_RK4,tita_punto_array_RK4, Tau, l, g);

        [tita_array_CNLF, tita_punto_array_CNLF] = CN_LeapFrog(y_ini, t_ini, h_array(ii), n_array(ii), a, b, Nmax_bis, tol_bis, alea);
        e_fase(ii,4) = error_fase(t_array, tita_array_CNLF,tita_punto_array_CNLF, Tau);
        e_amplitud(ii,4) = erroramplitud(t_array, tita_array_CNLF,tita_punto_array_CNLF, Tau, l, g);

        #Ploteo tita_array vs tiempo:
        % e_fase(ii,:)

        % t_array = linspace(0, n_array(ii)*h_array(ii), n_array(ii));
        % plot(t_array,tita_array_EI,";Euler Implícito;","linewidth", 2); hold on;
        % plot(t_array,tita_array_CN,";Crack - Nicholson;","linewidth", 2); hold on;
        % plot(t_array,tita_array_RK4,";RK4;","linewidth", 2); hold on;
        % plot(t_array,tita_array_CNLF,";CN + Leap Frog;","linewidth", 2); hold off;
        % pause(4)
        

    endfor

    % #Ploteo los errores vs h_array:
    % plot(h_array,e_fase(:,1),";Euler Implícito;","linewidth", 2);
    % hold on;
    % plot(h_array,e_fase(:,2),";Crack - Nicholson;","linewidth", 2);
    % hold on;
    % plot(h_array,e_fase(:,3),";RK4;","linewidth", 2);
    % hold on;
    % plot(h_array,e_fase(:,4),";CN + Leap Frog;","linewidth", 2);
    % pause(2)
    % hold off;

    % plot(h_array,e_amplitud(:,1),";Euler Implícito;","linewidth", 2);
    % hold on;
    % plot(h_array,e_amplitud(:,2),";Crack - Nicholson;","linewidth", 2);
    % hold on;
    % plot(h_array,e_amplitud(:,3),";RK4;","linewidth", 2);
    % hold on;
    % plot(h_array,e_amplitud(:,4),";CN + Leap Frog;","linewidth", 2);
    % pause(2)
    % hold off;

    #Ploteo los errores vs h_array en escala log-log
    loglog(h_array,abs(e_fase(:,1)),";Euler Implícito;","linewidth", 2);
    hold on;
    loglog(h_array,abs(e_fase(:,2)),";Crack - Nicholson;","linewidth", 2);
    hold on;
    loglog(h_array,abs(e_fase(:,3)),";RK4;","linewidth", 2);
    hold on;
    loglog(h_array,abs(e_fase(:,4)),";CN + Leap Frog;","linewidth", 2);
    pause(5)
    hold off;

    loglog(h_array,abs(e_amplitud(:,1)),";Euler Implícito;","linewidth", 2);
    hold on;
    loglog(h_array,abs(e_amplitud(:,2)),";Crack - Nicholson;","linewidth", 2);
    hold on;
    loglog(h_array,abs(e_amplitud(:,3)),";RK4;","linewidth", 2);
    hold on;
    loglog(h_array,abs(e_amplitud(:,4)),";CN + Leap Frog;","linewidth", 2);
    pause(5)
    hold off;

endif

#Agrego rectas con distintas dependencias
#Estaría bueno hacer una corrida grande y exportar los datos anteriores para luego graficarlos




#Inciso b: péndulo doble

function fb_vec =  fb_vec(y_vec, t)
    f1 = y_vec(3);
    f2 = y_vec(4);
    f3 = (y_vec(4)^2*sin(y_vec(2)-y_vec(1)) - 20*sin(y_vec(1)) + y_vec(3)^2*sin(y_vec(2)-y_vec(1))*cos(y_vec(2)-y_vec(1)) + 10*sin(y_vec(2))*cos(y_vec(2)-y_vec(1))  )  /  (  2 - cos(y_vec(2)-y_vec(1))^2  );
    f4 = -y_vec(3)^2*sin(y_vec(2)-y_vec(1)) -10*sin(y_vec(2)) - f3*cos(y_vec(2)-y_vec(1));

    % (v2^2*sin(u2-u1) - 20*sin(u1) + v1^2*sin(u2-u1)*cos(u2-u1) + 10*sin(u2)*cos(u2-u1)  )  /  (  2 - cos(u2-u1)^2  )
    % -v1^2*sin(u2-u1) -10*sin(u2) - (v2^2*sin(u2-u1) - 20*sin(u1) + v1^2*sin(u2-u1)*cos(u2-u1) + 10*sin(u2)*cos(u2-u1)  )  /  (  2 - cos(u2-u1)^2  )*cos(u2-u1)

    fb_vec = [f1,f2,f3,f4];
endfunction

function y_new1 = RK4b(y_old, t_old, h)
    # Dado y_n calcula y_{n+1} mediante el método de Runge-Kutta de orden 4
    #y_new1 es el vector [y1_{n+1}, y2_{n+1}]. No se coloca como input porque se trata de un método explícito
    #y_old es el vector [y1_n, y2_n]
    #t_old es el tiempo t_n
    #h es el paso de tiempo
    t_new = t_old + h;
    k1 = fb_vec(y_old,t_old);
    k2 = fb_vec(y_old + h/2*k1,t_old + h/2);
    k3 = fb_vec(y_old + h/2*k2,t_old + h/2);
    k4 = fb_vec(y_old + h*k3,t_old + h);
    y_new1 = y_old + h/6*(k1 + 2*k2 + 2*k3 + k4);
endfunction

function [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble(y_ini, t_ini, h, n)
    #Aplico un método explícito al pendulo doble y retorno tita y tita_punto para todo t discreto entre 0 y h*n
    #y_ini es el vector [y1_0, y2_0, y3_0, y4_0] = [tita1_0, tita2_0, tita1_punto_0, tita2_punto_0]
    #t_ini es el tiempo t_0
    #h es el paso de tiempo
    #n es el número de pasos de tiempo

    tita1_array = zeros(n,1);
    tita2_array = zeros(n,1);
    tita1_punto_array = zeros(n,1);
    tita2_punto_array = zeros(n,1);

    #Aplico el método n veces
    t_old = t_ini;
    y_ini_copia = y_ini;
    for ii = 1:n
        y_new = RK4b(y_ini_copia,t_old,h);
        tita1_array(ii) = y_new(1);
        tita2_array(ii) = y_new(2);
        tita1_punto_array(ii) = y_new(3);
        tita2_punto_array(ii) = y_new(4);
        y_ini_copia = [tita1_array(ii), tita2_array(ii),tita1_punto_array(ii), tita2_punto_array(ii)];
        t_old = t_old + h;
    endfor

endfunction

#Test de la solución para las condiciones de la cátedra
y_ini = [pi/2, pi/2, 0, 0];
t_ini = 0;
t_max = pi;
n = 1000; #discretización
h = (t_max - t_ini)/n;

[tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble(y_ini, t_ini, h, n);

#Resulta ser el mismo valor que la cátedra
[tita1_array(length(tita1_array)), tita2_array(length(tita1_array)), tita1_punto_array(length(tita1_array)), tita2_punto_array(length(tita1_array))]



#Estudio la solución para distintas condiciones iniciales hasta t = 10*pi
#Grafico las diferencias entre tita1, tita2, tita1_punto y tita2_punto para distintas condiciones iniciales respecto a una de ellas

plotear = false;
if plotear == true
    t_ini = 0;
    t_max = 10*pi;
    n = 1000; #discretización
    h = (t_max - t_ini)/n;
    t_array = linspace(0, n*h, n);

    a = pi/2;
    [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble([pi/2, a, 0, 0], t_ini, h, n);
    tita_array_A = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];

    a = 1.00001*pi/2;
    [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble([pi/2, a, 0, 0], t_ini, h, n);
    tita_array_B = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];

    a = 0.99999*pi/2;
    [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble([pi/2, a, 0, 0], t_ini, h, n);
    tita_array_C = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];

    #Grafico tita1 y tit2 de todos
    plot(t_array,tita_array_A(:,1),";a = pi/2;","linewidth", 1);
    hold on;
    plot(t_array,tita_array_B(:,1),";a = 1.00001*pi/2;","linewidth", 1);
    hold on;
    plot(t_array,tita_array_C(:,1),";a = 0.99999*pi/2;","linewidth", 1);
    pause(10);
    hold off;

    #Grafico las diferencias respecto a A
    % plot(t_array,tita_array_A(:,1) - tita_array_B(:,1),";a = 1.00001*pi/2;","linewidth", 1);
    % hold on;
    % plot(t_array,tita_array_A(:,1) - tita_array_C(:,1),";a = 0.99999*pi/2;","linewidth", 1);
    % pause(10);


endif



#Grafico las diferencias entre tita1, tita2, tita1_punto, tita2_punto para t = 5*pi para distintos valores de h
plotear = false;

if plotear == true
    t_ini = 0;
    t_max = 5*pi;
    n_array = [50,100,200,500,1000, 2000, 5000, 10000]; #discretización
    #n_array = [50,100,200,500,1000, 2000, 5000, 10000, 20000, 50000]; #discretización
    h_array = zeros(length(n_array),1);
    diferencias_B = zeros(length(n_array),4);
    diferencias_C = zeros(length(n_array),4);

    for ii=1:length(n_array)
        n = n_array(ii);
        h_array(ii) = (t_max - t_ini)/n;
        h = h_array(ii);

        a = pi/2;
        [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble([pi/2, a, 0, 0], t_ini, h, n);
        tita_array_A = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];

        a = 1.00001*pi/2;
        [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble([pi/2, a, 0, 0], t_ini, h, n);
        tita_array_B = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];

        a = 0.99999*pi/2;
        [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble([pi/2, a, 0, 0], t_ini, h, n);
        tita_array_C = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];

        diferencias_B(ii,:) = tita_array_A(length(tita_array_A(:,1)),:) - tita_array_B(length(tita_array_B(:,1)),:);
        diferencias_C(ii,:) = tita_array_A(length(tita_array_A(:,1)),:) - tita_array_C(length(tita_array_C(:,1)),:);
    endfor

    #Grafico las diferencias
    plot(h_array, diferencias_B(:,1),";B - tita1;","linewidth", 1); hold on;
    plot(h_array, diferencias_B(:,2),";B - tita2;","linewidth", 1); hold on;
    plot(h_array, diferencias_B(:,3),";B - tita1_punto;","linewidth", 1); hold on;
    plot(h_array, diferencias_B(:,4),";B - tita2_punto;","linewidth", 1); hold on;
    plot(h_array, diferencias_C(:,1),";C - tita1;","linewidth", 1); hold on;
    plot(h_array, diferencias_C(:,2),";C - tita2;","linewidth", 1); hold on;
    plot(h_array, diferencias_C(:,3),";C - tita1_punto;","linewidth", 1); hold on;
    plot(h_array, diferencias_C(:,4),";C - tita2_punto;","linewidth", 1); hold off;
    pause(100);
endif


#Grafico con la función comet (VER EJEMPLO DE CÁTEDRA)
t_ini = 0;
t_max = 10*pi;
n = 1000; #discretización
h = (t_max - t_ini)/n;
t_array = linspace(0, n*h, n);

a = pi/2;
[tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble([pi/2, a, 0, 0], t_ini, h, n);
tita_array_A = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];

a = 1.00001*pi/2;
[tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble([pi/2, a, 0, 0], t_ini, h, n);
tita_array_B = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];

a = 0.99999*pi/2;
[tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble([pi/2, a, 0, 0], t_ini, h, n);
tita_array_C = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];


title ("Trayectorias");
for l=1:n
    pendulo1_A = [sin(tita_array_A(l,1)), -cos(tita_array_A(l,1))];
    pendulo1_B = pendulo1_A + [sin(tita_array_A(l,2)), -cos(tita_array_A(l,2))];
    plot(pendulo1_A(1), pendulo1_A(2), 'r'); hold on;
    plot(pendulo1_B(1), pendulo1_B(2), 'r'); hold off;

    % plot(sin(uexacta(1,l)),-cos(uexacta(1,l)),sin(uee(1,l)),-cos(uee(1,l)),sin(uei(1,l)),-cos(uei(1,l)),sin(ucn(1,l)),-cos(ucn(1,l)),'linewidth',1.5);
    % legend("Exacto","EE","EI","CN");
    axis([-3 3 -3 3]); 
    pause(0.1);
endfor



























#Calculo el error de amplitud
function ampli_b = amplitud_b(tita_array, g)
    #Calcula la energía del sistema
    #tita_array = [tita1, tita2, tita1_punto, tita2_punto]
    m1 = 1;
    m2 = 1;
    l1 = 1;
    l2 = 1;

    tita1 = tita_array(1);
    tita2 = tita_array(2);
    tita1_punto = tita_array(3);
    tita2_punto = tita_array(4);

    T = 1/2*(m1+m2)*l1^2*tita1_punto^2 + 1/2*m2*l2^2*tita2_punto^2 + m2*l1*l2*tita1_punto*tita2_punto*cos(tita2-tita1);

    V = -(m1+m2)*g*l1*cos(tita1) - m2*g*l2*cos(tita2);

    ampli_b = T + V;
endfunction

function e_ampli_b_array = erroramplitud_b(t_array, tita_array,g)
    #Dado t_array, tita_array y tita_punto_array, calculo el error de amplitud (1/2l^2tita′^2 − gl cos(tita)) para todo tiempo

    tita1_array = tita_array(:,1);
    tita2_array = tita_array(:,2);
    tita1_punto_array = tita_array(:,3);
    tita2_punto_array = tita_array(:,4);

    amplitud_exacta = amplitud_b(tita_array(1,:), g);
    #Busco el tiempo t en t_array tal que t = Tau (no considero múltiplos)
    e_ampli_b_array = zeros(length(t_array),1);
    for ii=1:length(t_array)
        #Calculo el error de amplitud
        e_ampli_b_array(ii) = amplitud_exacta - amplitud_b(tita_array(ii,:), g);
    endfor
endfunction


#Calculo el error de amplitud para todo tiempo

plotear == false;
if plotear == true
    y_ini = [pi/2, pi/2, 0, 0];
    t_ini = 0;
    t_max = pi;
    n = 1000; #discretización
    h = (t_max - t_ini)/n;

    [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble(y_ini, t_ini, h, n);
    tita_array = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];
    t_array = linspace(0, n*h, n);

    e_ampli_b_array = erroramplitud_b(t_array, tita_array,g);

    #Grafico
    plot(t_array, e_ampli_b_array, 'r');
    pause(10);
endif

#Calculo el error a t = pi en función de h
plotear = false
if plotear == true
    y_ini = [pi/2, pi/2, 0, 0];
    t_ini = 0;
    t_max = pi;
    #n_array = [10,20,50,100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000];
    n_array = [10,20,50,100, 200, 500, 1000, 2000, 5000, 10000];
    h_array = zeros(length(n_array),1);
    e_ampli_b_array_h = zeros(length(n_array),1);
    for ii=1:length(n_array)
        n = n_array(ii);
        h_array(ii) = (t_max - t_ini)/n;
        h = h_array(ii);
        [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array] = Pendulo_doble(y_ini, t_ini, h, n);
        tita_array = [tita1_array, tita2_array, tita1_punto_array, tita2_punto_array];
        t_array = linspace(0, n*h, n);
        e_ampli_b_array_h(ii) = abs(amplitud_b(tita_array(1,:), g) - amplitud_b(tita_array(length(tita1_array),:), g));
    endfor

    #Grafico
    loglog(h_array, e_ampli_b_array_h, 'r');
    pause(10);
endif