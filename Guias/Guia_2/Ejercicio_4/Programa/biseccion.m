function x = biseccion(f, a, b, Nmax_bis, tol_bis)
    # Plantea el método de bisección para la ecuación f(x) = 0 considerando que f(a) y f(b) tienen signos opuestos
    # a y b son límites
    # x es el valor buscado
    # Nmax_bis es el nro de iteraciones
    # tol_bis es la tolerancia

    fa = f(a);
    fb = f(b);

    # Verifico si a o b ya son solución
    if (fa == 0)
        x = a;
        return
    endif
    if (fb == 0)
        x = b;
        return
    endif

    #Verifico que no tengan el mismo signo
    if (fa*fb>0)
        error("No se pudo ejecutar el metodo de biseccion")
    endif

    longitud = b-a; #longitud del intervalo
    c = (a+b)/2; #valor medio del intervalo

    fc = f(c);

    contador = 1; #útil para verificar que no se hagan más bisecciones que Nmax_bis
    while (contador <= Nmax_bis && longitud > tol_bis & fc != 0)
        aux = c;
        if fc*fa<0
            b = c;
            fb = fc;
        else
            a = c;
            fa = fc;
        endif
        c = (a + b) / 2;
        fc = f(c);
        longitud = abs(c-aux);
        contador = contador + 1;
    endwhile

    if (fc == 0) || (longitud <= tol_bis)
        x = c;
        return
    else
        error("en el metodo de biseccion se supero el numero maximo de iteraciones")
    endif

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