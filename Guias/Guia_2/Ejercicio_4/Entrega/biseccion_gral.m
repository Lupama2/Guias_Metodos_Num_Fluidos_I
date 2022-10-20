function [y1_new, y2_new] = biseccion_gral(metodo, y_old, t_old, h, a, b, Nmax_bis, tol_bis, alea)

    #Planteo el método de bisección para el problema en particular y retorna los valores de y1_new y y2_new
    #metodo es el método a utilizar. Se trata de un puntero a una función
    #y_old es el vector [y1_n, y2_n]
    #t_old es el tiempo t_n
    #h es el paso de tiempo


    #Defino la función a la que le aplico el método de bisección, es decir, f tal que f(y2_new) = 0
    function f_biseccion_= f_biseccion(y2_new)
        #Planteo la primer componente como una igualdad
        y1_new = metodo([alea,y2_new], y_old, t_old, h)(1); #independiente del valor de y1_new. Esto se cumple para EI y CN
        f_biseccion_ = y2_new - metodo([y1_new,y2_new], y_old, t_old, h)(2); #ec. a resolver con el método de biyección
    endfunction

    y2_new = biseccion(@f_biseccion, a, b, Nmax_bis, tol_bis);
    y1_new = metodo([alea,y2_new], y_old, t_old, h)(1);
endfunction

