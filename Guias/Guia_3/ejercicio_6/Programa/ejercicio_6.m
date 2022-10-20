% Ejercicio 6 de la guía 3 de Métodos Numéricos en Fluidos I
% Chehade Pablo

% El objetivo es resolver una ecuación en derivadas parciales numéricamente a través de RK4 y algún método espacial, analizando posteriormente las características de los mismos

a = 1; %Este statement es necesario para que octave no piense que este archivo es un archivo de funciones

function u0_ = u0(x)
    %Condición inicial: velocidad en x a t = 0 
    u0_ = exp(-200(x-0.25)^2);
endfunction

function c_ = c(x, inciso)
    %Velocidad del sonido en función de la posición para el inciso inciso
    switch (inciso)
        case "a"
            %Arenisca porosa
            c_ = 1;
        case "b"
            %Transición a arenisca impermeable
            c_ = 1.25 - 0.25*tanh(40*(0.75-x))
        case "c"
            %Arenisca impermeable
            c_ = 1.5;
        case "d"
            %Nave alienígena sepultada
            c_ = 1.5 - exp(-300*(x-1.75)^2);
        otherwise
            error("La velocidad c no fue asignada")
    endswitch
endfunction


% El método RK4 se encuentra definido en RK4.m




%Resuelvo el inciso a
inciso = "a";
t_ini = 0;
t_max = 8;

%Discretización
M = 5;
Delta_x = 1/(M+1);
N = 2*M;
Delta_t = 1/N;



%Creo el vector solución de tamaño 2N + 2
z = zeros(2*M + 2,1)
%Defino coeficientes para la resolución numérica
a = c(0,inciso)/(2*Delta_x);
b = c(4,inciso)/(2*Delta_x);
d = 1/Delta_x^2;

%Defino la matriz sparse gamma de tamaño (2M + 2) x (2M + 2)
size_ = 2*M + 2;

#1er cuadrante
gamma_1_1 = sparse([1],[1],-3*a, size_,size_,0);
gamma_1_2 = sparse([1],[2],4*a, size_,size_,0);
gamma_1_3 = sparse([1],[3],-a, size_,size_,0);
gamma_1_4 = sparse([M+2],[M],-b, size_,size_,0);
gamma_1_5 = sparse([M+2],[M+1],4*b, size_,size_,0);
gamma_1_6 = sparse([M+2],[M+2],-3*b, size_,size_,0);
%Sumo
gamma_1 = gamma_1_1 + gamma_1_2 + gamma_1_3 + gamma_1_4 + gamma_1_5 + gamma_1_6;

#2do cuadrante
gamma_2 = sparse([2:M+1],[1+ (M+2):M + (M+2)], ones(1,M), size_,size_,0)

#3er cuadrante
gamma_3 = sparse([1],[1],0,size_,size_,0); %Creo un sparse de ceros
for ii = 1:M
    cii = c(ii*Delta_x, inciso);
    gamma_3_1 = sparse([ii + (M+2)],[ii], d*cii^2, size_,size_, 0);
    gamma_3_2 = sparse([ii + (M+2)],[ii+1],-2*cii^2*d, size_, size_, 0);
    gamma_3_3 = sparse([ii + (M+2)],[ii+1],d*cii^2, size_, size_, 0);
    gamma_3 = gamma_3 + gamma_3_1 + gamma_3_2 + gamma_3_3;
endfor
gamma_3

#4to cuadrante: son nulos

gamma_matriz = gamma_1 + gamma_2 + gamma_3;
#VERIFICAR QUE ESTÉ BIEN DESCRIPTA

% #Nodos internos:
% P_S1 = sparse([2:Nn-1], [2:Nn-1], (10/12)*ones(1,Nn-2), N, Nn, 0);
% P_S2 = sparse([2:Nn-1], [1:Nn-2], (1/12)*ones(1,Nn-2), Nn, Nn, 0);
% P_S3 = sparse([2:Nn-1], [3:Nn],(1/12)*ones(1,Nn-2), Nn, Nn, 0);

% #Nodos de frontera:
% P_F1 = sparse([1], [1], (1), Nn, Nn, 0);
% P_F2 = sparse([Nn], [Nn], (1), Nn, Nn, 0);

% A2 = P_S1 + P_S2 + P_S3 + P_F1 + P_F2;



%A partir de la matriz gamma defino la función f(z) tal que dz/dt = f(z)
%Esta está definida como f(z) = gamma * z
function f(z)

endfunction

%Aplico el método RK4 y ploteo la solución u en cada paso de tiempo



% Grafico c(x)