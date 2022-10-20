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