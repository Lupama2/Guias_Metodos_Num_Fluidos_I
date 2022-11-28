#---------------------------------------------
#   SOLVER NAVIER-STOKES 2D CAVIDAD CUADRADA HIDRODINAMICA CON GRILLA UNIFORME STAGERRED - ADIMENSIONAL CON U Y L (Re=U*L/nu)
#       METOD0 SIMPLER
#       TIEMPO: EULER IMPLICITO          
#       ADVECCION: VELOCIDAD EN LAS CARAS CON DC2, UPWIND THE PRIMER ORDEN + CORRECCION DIFERIDA PARA VELOCIDAD CON DC2  
#       DIFUSION: DC2 
#---------------------------------------------

#---------------------------------------------
#   NOTA 1: Si se quiere continuar una simulacion comentar el clear (# clear)
#---------------------------------------------
% clear

function Chehade_NSdc2(Re, U0_top, n1, dt, tol_estacionario, Ndeltat, nsimpler, termino_advectivo, metodo_temporal, archivo_evolucion_variables, archivo_velocidades_centrales, archivo_parametros)
#Re: numero de Re basado en Utop L / nu
#n1: Numero de volumenes por direccion.
#dt: delta tiempo adimensional t*U/L
#nsimpler: nro de iteraciones internas en SIMPLER
#termino_advectivo: esquema para el termino advectivo. Posibilidades: "UP1", "DC2" o "QUI" haciendo referencia a "QUICK".
#U0_top: función velocidad U0(t) en el borde superior 
#metodo_temporal : "EI" (Euler Implícito) o "CN" (Crank-Nicolson). Técnicamente si no se usa "CN" se está usando "EI". En ningún lado metí un if con "EI". 
#tol_estacionario. Condición de estado estacionario (agregar una función que elija Ndeltat de modo de correr el código hasta que se llegue a una determinada condición)
#N_deltat. Nro máximo de pasos de tiempo. Condición adicional por las dudas. Estaría bueno imprimir en la terminal qué condición se cumplió para terminar la simulación
#archivo_evolucion_variables. Nombre del archivo donde se guarda la evolución de las variables.
#archivo_velocidades_centrales. Nombre del archivo donde se guarda la evolución de las velocidades en los puntos centrales.
#archivo_parametros. Nombre del archivo donde se guardan algunos parámetros

#---------------------------------------------
#Datos del problema
#---------------------------------------------

% Re = 100.0; #X numero de Re basado en Utop L / nu. DEFINIDO COMO INPUT
L = 1.0; #Dimension de la cabidad cuadrada adimensionalisada

function Utop1 = Utop1(k)
  #velocidad top normalizada = 1
  #k: índice de la iteración
  #Ndeltat: número de pasos de tiempo
  t = k*dt;
  Utop1= U0_top(t);
endfunction

% Utop1 = 1.0;
Utop2 = 0.0; #velocidad right normalizada = 1
Utop3 = 0.0; #velocidad bottom normalizada = 1
Utop4 = 0.0; #velocidad left normalizada = 1

#---------------------------------------------
#Datos para la discretizacion
#---------------------------------------------

% n1 = 20; #Numero de volumenes por direccion. DEFINIDO COMO INPUT.
% dt=0.5; #delta tiempo adimensional t*U/L. DEFINIDO COMO INPUT
% Ndeltat=80; #numero de pasos de tiempo

#---------------------------------------------
#Condición límite para detener la evolución temporal
#---------------------------------------------
function criterio = condicion_limite_temporal(derivadatuvel_k, derivadatvvel_k)
  % Evalúa si se llegó o no al estado estacionario empleando como criterio que las derivadas de las velocidades en el tiempo sean menores a una tolerancia
  % derivadatuvel(k)=norm(uvel-u0)/dt;
  % derivadatvvel(k)=norm(vvel-v0)/dt;
  derivadatuvel_k/tol_estacionario #Para imprimir
  if derivadatuvel_k < tol_estacionario && derivadatvvel_k < tol_estacionario
    criterio = true;
  else
    criterio = false;
  endif
endfunction

#---------------------------------------------
#Auxiliares
#---------------------------------------------

n2=n1*n1; #Numero de volumenes totales
dx = L / n1; #lado del volumen = distancia entre centros
dx2 = dx*dx; #volumen del volumen
re1=1/Re; #inversa numero de Reynolds
#---------------------------------------------
#   NOTA 1: Si se quiere continuar una simulacion comentar el tiempo (#tiempo=0)
#---------------------------------------------
tiempo=0;
t1=time();

#---------------------------------------------
# Datos para guardar evolucion temporal de variables
#---------------------------------------------

i1=2;
j1=2;

#---------------------------------------------
#--------- abro archivo de evolucion de variables
fid = fopen(archivo_evolucion_variables, "w");

#---------------------------------------------
#Reservar Vectores para solucion y definimos condicion inicial = todo quieto
#---------------------------------------------

# ------------------------------------------
# Variables para la ecuacion de momento en x
# ------------------------------------------
#---------------------------------------------
#   NOTA 1: Si se quiere continuar una simulacion comentar la velocidad u (# uvel = zeros(n1-1,n1))
#---------------------------------------------
uvel = zeros(n1-1,n1); # velocidad-x forma matricial (n1,n1)
u=zeros(n2-n1,1); # velocidad-x forma vectorial (n2)
u0 = zeros(n1-1,n1); # velocidad-x forma matricial paso anterior (n1,n1)
#  Coeficientes matriz momento en x 
eu = zeros(n1-1,n1); # coeficiente al este 
wu = zeros(n1-1,n1); # coeficiente al oeste
nu = zeros(n1-1,n1); # coeficiente al norte
su = zeros(n1-1,n1); # coeficiente al sur
au = zeros(n1-1,n1); # coeficiente diagonal
bu = zeros(n1-1,n1); # lado izquierdo 
# Velocidades en las caras
uwest = zeros(n1,n1); # velocidad en la cara oeste para el volumen de u
usouth = zeros(n1-1,n1+1); # velocidad en la cara sur para el volumen de u
# Factores de upwind en las caras
ufw = zeros(n1,n1); # factor en la cara oeste
ufs = zeros(n1-1,n1); # factor en la cara sur
# Fuentes diferidas de upwind en las caras
usourcew = zeros(n1,n1); # fuente en la cara oeste
usources = zeros(n1-1,n1); # fuente en la cara sur
# Matriz momento en x 
MMx = zeros(n2-n1,n2-n1); # matriz momento x
bux = zeros(n2-n1,1); # lado derecho MMx * u = bux


# ------------------------------------------
# Variables para la ecuacion de momento en y
# ------------------------------------------
#---------------------------------------------
#   NOTA 1: Si se quiere continuar una simulacion comentar la velocidad u (# vvel = zeros(n1-1,n1))
#---------------------------------------------
vvel = zeros(n1,n1-1); # velocidad-y forma matricial (n1,n1)
v=zeros(n2-n1,1); # velocidad-y forma vectorial (n2)
v0 = zeros(n1,n1-1); # velocidad-y forma matricial (n1,n1) paso anterior
#  Coeficientes matriz momento en y 
ev = zeros(n1,n1-1); # coeficiente al este 
wv = zeros(n1,n1-1); # coeficiente al oeste
nv = zeros(n1,n1-1); # coeficiente al norte
sv = zeros(n1,n1-1); # coeficiente al sur
av = zeros(n1,n1-1); # coeficiente diagonal
bv = zeros(n1,n1-1); # lado izquierdo 
# Velocidades en las caras
vwest = zeros(n1+1,n1-1); # velocidad en la cara oeste para el volumen de v
vsouth = zeros(n1,n1); # velocidad en la cara sur para el volumen de v
# Factores de upwind en las caras
vfw = zeros(n1,n1-1); # factor en la cara oeste
vfs = zeros(n1,n1); # factor en la cara sur
# Fuentes diferidas de upwind en las caras
vsourcew = zeros(n1,n1-1); # fuente en la cara oeste
vsources = zeros(n1,n1); # fuente en la cara sur
# Matriz momento en x 
MMy = zeros(n2-n1,n2-n1); # matriz momento y
buy = zeros(n2-n1,1); # lado derecho MMy * v = buy

# ------------------------------------------
#Inicializo las matrices nulas b_CNu y b_CNv para aplicar el método de Crank-Nicolson en función de la recomendación de la cátedra en el pdf del TP_final. También divido por la mitad el paso temporal en caso de ser CN
# ------------------------------------------
if metodo_temporal == "C"


  b_CNu_matriz = zeros(n1-1,n1); #matriz actual (n1,n1)
  b_CNu0_matriz = zeros(n1-1,n1); #matriz del paso anterior (n1,n1)
  b_CNv_matriz = zeros(n1,n1-1); # matriz actual (n1,n1)
  b_CNv0_matriz = zeros(n1,n1-1); # matriz del paso anterior (n1,n1)

  #Redefino dt y Ndeltat
  Ndeltat = 2*Ndeltat;
  dt = dt/2;

endif



# ------------------------------------------
# Variables para la presion
# ------------------------------------------
#  Coeficientes matriz presion
du = zeros(n1-1,n1); # coeficiente en x
dv = zeros(n1,n1-1); # coeficiente en y
Usom = zeros(n1-1,n1); # lado derecho en x
Vsom = zeros(n1,n1-1); # lado derecho en y
# Matriz presion
MP = zeros(n2,n2); # matriz presion
bp = zeros(n2,1); # lado derecho MP * pres = bp
pres=zeros(n2,1); # presion en forma vectorial (n2)
prescor=zeros(n2,1); # correcion de presion en forma vectorial (n2)
presion=zeros(n1,n1); # presion en forma matricial

#---------------------------------------------
#---------------------------------------------
#
#  COMIENZA LOOP DE ITERATION TEMPORAL - METODO SIMPLER (PATANKAR 1982)
#
#
for k=1:Ndeltat   # Este es el loop de tiempo!!!
k #Para imprimir
n1;
Re;
dt;
#
#
#---------------------------------------------
#---------------------------------------------

#---------------------------------------------
# Actualiza la velocidad u del tiempo anterior
#---------------------------------------------

for i=1:n1-1
  for j=1:n1-1
    u0(i,j)=uvel(i,j);
    v0(i,j)=vvel(i,j);
  endfor
  j=n1;
  u0(i,j)=uvel(i,j);
endfor
i=n1;
for j=1:n1-1
  v0(i,j)=vvel(i,j);
endfor

if metodo_temporal == "C"
  uvel = uvel + b_CNu_matriz*dt;
  vvel = vvel + b_CNv_matriz*dt;
endif


for lsimpler=1:nsimpler     #Este es el loop de iteraciones simpler!!!

#---------------------------------------------
# Calculo de velocidades en los bordes de los Volumenes, se usa atrasada para eliminar la parte no lineal
#---------------------------------------------

#-  Velocidad en x al oeste y al sur de u(i,j) Diferencias centradas
i=1;
j=1;
uwest(i,j)=uvel(i,j)/2;
usouth(i,j)=0;
for j=2:n1
  uwest(i,j)=uvel(i,j)/2;
  usouth(i,j)=(vvel(i,j-1)+vvel(i+1,j-1))/2;
endfor
j=n1+1;
usouth(i,j)=0;
for i=2:n1-1
  j=1;
  uwest(i,j)=(uvel(i,j)+uvel(i-1,j))/2;
  usouth(i,j)=0;
  for j=2:n1
    uwest(i,j)=(uvel(i,j)+uvel(i-1,j))/2;
    usouth(i,j)=(vvel(i,j-1)+vvel(i+1,j-1))/2;
  endfor
  j=n1+1;
  usouth(i,j)=0;
endfor
i=n1;
for j=1:n1
  uwest(i,j)=uvel(i-1,j)/2;
endfor


#-  Velocidad en y al sur y al oeste de v(i,j) Diferencias centradas

i=1;
j=1;
vwest(i,j)=0;
vsouth(i,j)=vvel(i,j)/2;
for j=2:n1-1
  vwest(i,j)=0;
  vsouth(i,j)=(vvel(i,j)+vvel(i,j-1))/2;
endfor
j=n1;
vsouth(i,j)=vvel(i,j-1)/2;
for i=2:n1
  j=1;
  vwest(i,j)=(uvel(i-1,j)+uvel(i-1,j+1))/2;
  vsouth(i,j)=vvel(i,j)/2;
  for j=2:n1-1
    vwest(i,j)=(uvel(i-1,j)+uvel(i-1,j+1))/2;
    vsouth(i,j)=(vvel(i,j)+vvel(i,j-1))/2;
  endfor
  j=n1;
  vsouth(i,j)=vvel(i,j-1)/2;
endfor
i=n1+1;
for j=1:n1-1
  vwest(i,j)=0;
endfor

#---------------------------------------------
# Calculo de los factores de upwind basados en la velocidad calculada en las caras
#---------------------------------------------
#La siguiente cuenta es para usar el esquema UP1 para el término advectivo. Como DC2 se escribe como UP1 + fuente, es necesario correrlo tmb si queremos aplicar DC2. Lo mismo ocurre para QUICK que tmb se escribe como UP1 + fuente
if termino_advectivo == "U" || termino_advectivo == "D"  || termino_advectivo == "Q"
  #   PARA VOLUMENES U
  for i=1:n1-1
    j=1;
    if (uwest(i,j)>0)
      ufw(i,j)=1;
    else
      ufw(i,j)=0;
    endif
    ufs(i,j)=0; 
    for j=2:n1
    if (uwest(i,j)>0)
            ufw(i,j)=1;
          else
            ufw(i,j)=0;
          endif 
    if (usouth(i,j)>0)
            ufs(i,j)=1;
          else
            ufs(i,j)=0;
          endif 
    endfor
  endfor

  i=n1;
  j=1;
  if (uwest(i,j)>0)
    ufw(i,j)=1;
  else
    ufw(i,j)=0;
  endif 
  for j=2:n1
    if (uwest(i,j)>0)
      ufw(i,j)=1;
    else
      ufw(i,j)=0;
    endif 
  endfor

  #   PARA VOLUMENES V
  i=1;
  for j=1:n1-1
    if (vsouth(i,j)>0)
    vfs(i,j)=1;
    else
      vfs(i,j)=0;
    endif 
    vfw(i,j)=0;
  endfor

  j=n1;
  if (vsouth(i,j)>0)
    vfs(i,j)=1;
  else
    vfs(i,j)=0;
  endif

  for i=2:n1
    for j=1:n1-1
    if (vwest(i,j)>0)
            vfw(i,j)=1;
          else
            vfw(i,j)=0;
          endif 
    if (vsouth(i,j)>0)
            vfs(i,j)=1;
          else
            vfs(i,j)=0;
          endif 
    endfor
    j=n1;
    if (vsouth(i,j)>0)
      vfs(i,j)=1;
    else
      vfs(i,j)=0;
    endif 
  endfor

endif

#---------------------------------------------
# Calculo de las fuentes diferidas debido al upwind para DC2
#---------------------------------------------
#La siguiente cuenta corresponde a la fuente en el esquema DC2 diferido usando UPWIND
if termino_advectivo == "D"
  #   PARA VOLUMENES U
  i=1;
  j=1;
  if (uwest(i,j)>0)
    usourcew(i,j)=(uvel(i,j))/2*dx*uwest(i,j);
  else
    usourcew(i,j)=(-uvel(i,j))/2*dx*uwest(i,j);
  endif 
  usources(i,j)=0;
  for j=2:n1
    if (uwest(i,j)>0)
      usourcew(i,j)=(uvel(i,j))/2*dx*uwest(i,j);
    else
      usourcew(i,j)=(-uvel(i,j))/2*dx*uwest(i,j);
    endif 
    if (usouth(i,j)>0)
      usources(i,j)=(uvel(i,j)-uvel(i,j-1))/2*dx*usouth(i,j);
    else
      usources(i,j)=(-uvel(i,j)+uvel(i,j-1))/2*dx*usouth(i,j);
    endif  
  endfor

  for i=2:n1-1
    j=1;
    if (uwest(i,j)>0)
      usourcew(i,j)=(uvel(i,j)-uvel(i-1,j))/2*dx*uwest(i,j);
    else
      usourcew(i,j)=(-uvel(i,j)+uvel(i-1,j))/2*dx*uwest(i,j);
    endif 
    usources(i,j)=0;
    for j=2:n1
      if (uwest(i,j)>0)
        usourcew(i,j)=(uvel(i,j)-uvel(i-1,j))/2*dx*uwest(i,j);
      else
        usourcew(i,j)=(-uvel(i,j)+uvel(i-1,j))/2*dx*uwest(i,j);
      endif 
      if (usouth(i,j)>0)
        usources(i,j)=(uvel(i,j)-uvel(i,j-1))/2*dx*usouth(i,j);
      else
        usources(i,j)=(-uvel(i,j)+uvel(i,j-1))/2*dx*usouth(i,j);
      endif  
    endfor
  endfor

  i=n1;
  j=1;
  if (uwest(i,j)>0)
    usourcew(i,j)=(-uvel(i-1,j))/2*dx*uwest(i,j);
  else
    usourcew(i,j)=(+uvel(i-1,j))/2*dx*uwest(i,j);
  endif 
  for j=2:n1
    if (uwest(i,j)>0)
      usourcew(i,j)=(-uvel(i-1,j))/2*dx*uwest(i,j);
    else
      usourcew(i,j)=(+uvel(i-1,j))/2*dx*uwest(i,j);
    endif 
  endfor


  #   PARA VOLUMENES V
  i=1;
  j=1;
  if (vsouth(i,j)>0)
    vsources(i,j)=(vvel(i,j))/2*dx*vsouth(i,j);
  else
    vsources(i,j)=(-vvel(i,j))/2*dx*vsouth(i,j);     
  endif  
  vsourcew(i,j)=0;   
  for j=2:n1-1
    if (vsouth(i,j)>0)
      vsources(i,j)=(vvel(i,j)-vvel(i,j-1))/2*dx*vsouth(i,j);         
    else
      vsources(i,j)=(-vvel(i,j)+vvel(i,j-1))/2*dx*vsouth(i,j);     
    endif
    vsourcew(i,j)=0;     
  endfor
  j=n1;
  if (vsouth(i,j)>0)
    vsources(i,j)=(-vvel(i,j-1))/2*dx*vsouth(i,j);          
  else
    vsources(i,j)=(+vvel(i,j-1))/2*dx*vsouth(i,j);     
  endif 

  for i=2:n1
    j=1;
    if (vwest(i,j)>0)
      vsourcew(i,j)=(vvel(i,j)-vvel(i-1,j))/2*dx*vwest(i,j);
    else
      vsourcew(i,j)=(-vvel(i,j)+vvel(i-1,j))/2*dx*vwest(i,j);
    endif 
    if (vsouth(i,j)>0)
      vsources(i,j)=(vvel(i,j))/2*dx*vsouth(i,j);          
    else
      vsources(i,j)=(-vvel(i,j))/2*dx*vsouth(i,j);    
    endif     
    for j=2:n1-1
      if (vwest(i,j)>0)
        vsourcew(i,j)=(vvel(i,j)-vvel(i-1,j))/2*dx*vwest(i,j);     
      else
        vsourcew(i,j)=(-vvel(i,j)+vvel(i-1,j))/2*dx*vwest(i,j);
      endif 
      if (vsouth(i,j)>0)
        vsources(i,j)=(vvel(i,j)-vvel(i,j-1))/2*dx*vsouth(i,j);          
      else
        vsources(i,j)=(-vvel(i,j)+vvel(i,j-1))/2*dx*vsouth(i,j);     
      endif     
    endfor
    j=n1;
    if (vsouth(i,j)>0)
      vsources(i,j)=(-vvel(i,j-1))/2*dx*vsouth(i,j);          
    else
      vsources(i,j)=(+vvel(i,j-1))/2*dx*vsouth(i,j);     
    endif 
  endfor

endif

#---------------------------------------------
# Calculo con QUICK
#---------------------------------------------
if termino_advectivo == "Q"
  % error("Aun no implemente QUICK para termino advectivo")
  #   PARA VOLUMENES U
  #Esquina inferior izquierda:
  i=1;
  j=1;
  if (uwest(i,j)>0)
    usourcew(i,j)=(uvel(i,j))/2*dx*uwest(i,j);
  else
    usourcew(i,j)=(-uvel(i,j))/2*dx*uwest(i,j);
  endif 
  usources(i,j)=0;

  #Borde izquierdo:
  for j=2:n1
    if (uwest(i,j)>0)
      usourcew(i,j)=(uvel(i,j))/2*dx*uwest(i,j);
    else
      usourcew(i,j)=(-uvel(i,j))/2*dx*uwest(i,j);
    endif 
    if (usouth(i,j)>0)
      usources(i,j)=(uvel(i,j)-uvel(i,j-1))/2*dx*usouth(i,j);
    else
      usources(i,j)=(-uvel(i,j)+uvel(i,j-1))/2*dx*usouth(i,j);
    endif  
  endfor

  for i=2:n1-1
    #Borde inferior:
    j=1;
    if (uwest(i,j)>0)
      usourcew(i,j)=(uvel(i,j)-uvel(i-1,j))/2*dx*uwest(i,j);
    else
      usourcew(i,j)=(-uvel(i,j)+uvel(i-1,j))/2*dx*uwest(i,j);
    endif 
    usources(i,j)=0;

    #Interior desde i=2 hasta n1-1 y desde j=2 hasta n1:
    #Ejecuto como si fuera uno de DC-2. Luego voy a cambiar algunos términos de fuente para QUICK
    for j=2:n1
      if (uwest(i,j)>0)
        usourcew(i,j)=(uvel(i,j)-uvel(i-1,j))/2*dx*uwest(i,j);
      else
        usourcew(i,j)=(-uvel(i,j)+uvel(i-1,j))/2*dx*uwest(i,j);
      endif 
      if (usouth(i,j)>0)
        usources(i,j)=(uvel(i,j)-uvel(i,j-1))/2*dx*usouth(i,j);
      else
        usources(i,j)=(-uvel(i,j)+uvel(i,j-1))/2*dx*usouth(i,j);
      endif  
    endfor
  endfor

  #Interior - QUICK
  #Corrijo los valores de fuente para usar QUICK en lugar de DC2
  for i=3:n1-2
      for j=3:n1-1
      if (uwest(i,j)>0)
        usourcew(i,j)=  1/8*(-uvel(i-2,j)-2*uvel(i-1,j) + 3*uvel(i,j)) *dx*uwest(i,j);
      else
        usourcew(i,j)= 1/8*(-uvel(i+1,j)-2*uvel(i,j) + 3*uvel(i-1,j))*dx*uwest(i,j);
      endif 
      if (usouth(i,j)>0)
        usources(i,j)=   1/8*(-uvel(i,j-2)-2*uvel(i,j-1) + 3*uvel(i,j))*dx*usouth(i,j);
      else
        usources(i,j)= 1/8*(-uvel(i,j+1)-2*uvel(i,j) + 3*uvel(i,j-1))*dx*usouth(i,j);
      endif  
    endfor
  endfor

  #Esquina inferior derecha:
  i=n1;
  j=1;
  if (uwest(i,j)>0)
    usourcew(i,j)=(-uvel(i-1,j))/2*dx*uwest(i,j);
  else
    usourcew(i,j)=(+uvel(i-1,j))/2*dx*uwest(i,j);
  endif 

  #Borde derecho:
  for j=2:n1
    if (uwest(i,j)>0)
      usourcew(i,j)=(-uvel(i-1,j))/2*dx*uwest(i,j);
    else
      usourcew(i,j)=(+uvel(i-1,j))/2*dx*uwest(i,j);
    endif 
  endfor


  #   PARA VOLUMENES V
  #Esquina inferior izquierda:
  i=1;
  j=1;
  if (vsouth(i,j)>0)
    vsources(i,j)=(vvel(i,j))/2*dx*vsouth(i,j);
  else
    vsources(i,j)=(-vvel(i,j))/2*dx*vsouth(i,j);     
  endif  
  vsourcew(i,j)=0;   

  #Borde izquierdo
  for j=2:n1-1
    if (vsouth(i,j)>0)
      vsources(i,j)=(vvel(i,j)-vvel(i,j-1))/2*dx*vsouth(i,j);         
    else
      vsources(i,j)=(-vvel(i,j)+vvel(i,j-1))/2*dx*vsouth(i,j);     
    endif
    vsourcew(i,j)=0;     
  endfor

  #Esquina superior izquierda
  j=n1;
  if (vsouth(i,j)>0)
    vsources(i,j)=(-vvel(i,j-1))/2*dx*vsouth(i,j);          
  else
    vsources(i,j)=(+vvel(i,j-1))/2*dx*vsouth(i,j);     
  endif 


  for i=2:n1
    #Borde inferior:
    j=1;
    if (vwest(i,j)>0)
      vsourcew(i,j)=(vvel(i,j)-vvel(i-1,j))/2*dx*vwest(i,j);
    else
      vsourcew(i,j)=(-vvel(i,j)+vvel(i-1,j))/2*dx*vwest(i,j);
    endif 
    if (vsouth(i,j)>0)
      vsources(i,j)=(vvel(i,j))/2*dx*vsouth(i,j);          
    else
      vsources(i,j)=(-vvel(i,j))/2*dx*vsouth(i,j);    
    endif     

    #Interior y creo borde derecho
    
    #Interior desde i=2 hasta n1 y desde j=2 hasta n1-1:
    #Ejecuto como si fuera uno de DC-2. Luego voy a cambiar algunos términos de fuente para QUICK 
    for j=2:n1-1
      if (vwest(i,j)>0)
        vsourcew(i,j)=(vvel(i,j)-vvel(i-1,j))/2*dx*vwest(i,j);     
      else
        vsourcew(i,j)=(-vvel(i,j)+vvel(i-1,j))/2*dx*vwest(i,j);
      endif 
      if (vsouth(i,j)>0)
        vsources(i,j)=(vvel(i,j)-vvel(i,j-1))/2*dx*vsouth(i,j);          
      else
        vsources(i,j)=(-vvel(i,j)+vvel(i,j-1))/2*dx*vsouth(i,j);     
      endif     
    endfor



    #Borde superior:
    j=n1;
    if (vsouth(i,j)>0)
      vsources(i,j)=(-vvel(i,j-1))/2*dx*vsouth(i,j);          
    else
      vsources(i,j)=(+vvel(i,j-1))/2*dx*vsouth(i,j);     
    endif 
  endfor


  #Interior-corrijo QUICK:
  for i=3:n1-1
    for j=3:n1-2
      if (vwest(i,j)>0)
        vsourcew(i,j)=1/8*(-vvel(i-2,j)-2*vvel(i-1,j) + 3*vvel(i,j))*dx*vwest(i,j);     
      else
        vsourcew(i,j)= 1/8*(-vvel(i+1,j)-2*vvel(i,j) + 3*vvel(i-1,j))*dx*vwest(i,j);
      endif 
      if (vsouth(i,j)>0)
        vsources(i,j)=1/8*(-vvel(i,j-2)-2*vvel(i,j-1) + 3*vvel(i,j))*dx*vsouth(i,j);          
      else
        vsources(i,j)= 1/8*(-vvel(i,j+1)-2*vvel(i,j) + 3*vvel(i,j-1))*dx*vsouth(i,j);     
      endif     
    endfor
  endfor


endif

# -------------------------------------------------------------------
#
#   PASO 1 SIMPLER: CALCULO EL CAMPO DE PRESIONES TENIENDO LAS VELOCIDADES COMO DATO
#
# -------------------------------------------------------------------

#---------------------------------------------
# Calculo de los coeficientes de la ecuacion de momento en x, todo del lado izquierdo salvo bu(i,j)
#---------------------------------------------

for i=1:n1-1
  j=1;
  eu(i,j)=dx*uwest(i+1,j)*(1-ufw(i+1,j))-re1;
  wu(i,j)=-dx*uwest(i,j)*ufw(i,j)-re1;
  nu(i,j)=dx*usouth(i,j+1)*(1-ufs(i,j+1))-re1;
  su(i,j)=0; # Definido como cero
  au(i,j)=dx2/dt+2*re1-eu(i,j)-wu(i,j)-nu(i,j)-su(i,j); # Tension de corte en la pared con primer orden
  bu(i,j)=2*Utop3*re1+dx2/dt*u0(i,j)+usourcew(i,j)-usourcew(i+1,j)-usources(i,j+1);
  #Si estoy en CN, agrego el término b de la propuesta de la catedra en el pdf TP_Final
  % if metodo_temporal == "C"
  %   bu(i,j) = bu(i,j) + b_CNu_matriz(i,j);
  % endif

  for j=2:n1-1
    eu(i,j)=dx*uwest(i+1,j)*(1-ufw(i+1,j))-re1;
    wu(i,j)=-dx*uwest(i,j)*ufw(i,j)-re1;
    nu(i,j)=dx*usouth(i,j+1)*(1-ufs(i,j+1))-re1;
    su(i,j)=-dx*usouth(i,j)*ufs(i,j)-re1;
    au(i,j)=dx2/dt-eu(i,j)-wu(i,j)-nu(i,j)-su(i,j);
    bu(i,j)=dx2/dt*u0(i,j)+usourcew(i,j)-usourcew(i+1,j)+usources(i,j)-usources(i,j+1);
    #Si estoy en CN, agrego el término b de la propuesta de la catedra en el pdf TP_Final
    % if metodo_temporal == "C"
    %   bu(i,j) = bu(i,j) + b_CNu_matriz(i,j);
    % endif
  endfor

  j=n1;
  eu(i,j)=dx*uwest(i+1,j)*(1-ufw(i+1,j))-re1;
  wu(i,j)=-dx*uwest(i,j)*ufw(i,j)-re1;
  nu(i,j)=0; # Definido como cero
  su(i,j)=-dx*usouth(i,j)*ufs(i,j)-re1;
  au(i,j)=dx2/dt+2*re1-eu(i,j)-wu(i,j)-nu(i,j)-su(i,j); # Tension de corte en la pared con primer orden
  bu(i,j)=2*Utop1(k)*re1+dx2/dt*u0(i,j)+usourcew(i,j)-usourcew(i+1,j)+usources(i,j);
  #Si estoy en CN, agrego el término b de la propuesta de la catedra en el pdf TP_Final
  % if metodo_temporal == "C"
  %   bu(i,j) = bu(i,j) + b_CNu_matriz(i,j);
  % endif

endfor

% if metodo_temporal == "C"
%   bu = bu + b_CNu_matriz;
% endif

#---------------------------------------------
# Calculo de los coeficientes de la ecuacion de momento en y, todo del lado izquierdo salvo bv(i,j)
#---------------------------------------------

i=1;
for j=1:n1-1
  ev(i,j)=dx*vwest(i+1,j)*(1-vfw(i+1,j))-re1;
  wv(i,j)=0;  # Definido como cero
  nv(i,j)=dx*vsouth(i,j+1)*(1-vfs(i,j+1))-re1;
  sv(i,j)=-dx*vsouth(i,j)*vfs(i,j)-re1;
  av(i,j)=dx2/dt+2*re1-ev(i,j)-wv(i,j)-nv(i,j)-sv(i,j); # Tension de corte en la pared con primer orden
  bv(i,j)=2*Utop4*re1+dx2/dt*v0(i,j)-vsourcew(i+1,j)+vsources(i,j)-vsources(i,j+1);
  #Si estoy en CN, agrego el término b de la propuesta de la catedra en el pdf TP_Final
  % if metodo_temporal == "C"
  %   bv(i,j) = bv(i,j) + b_CNv_matriz(i,j);
  % endif

endfor

for i=2:n1-1
  for j=1:n1-1
    ev(i,j)=dx*vwest(i+1,j)*(1-vfw(i+1,j))-re1;
    wv(i,j)=-dx*vwest(i,j)*vfw(i,j)-re1;
    nv(i,j)=dx*vsouth(i,j+1)*(1-vfs(i,j+1))-re1;
    sv(i,j)=-dx*vsouth(i,j)*vfs(i,j)-re1;
    av(i,j)=dx2/dt-ev(i,j)-wv(i,j)-nv(i,j)-sv(i,j);
    bv(i,j)=dx2/dt*v0(i,j)+vsourcew(i,j)-vsourcew(i+1,j)+vsources(i,j)-vsources(i,j+1);
    #Si estoy en CN, agrego el término b de la propuesta de la catedra en el pdf TP_Final
    % if metodo_temporal == "C"
    %   bv(i,j) = bv(i,j) + b_CNv_matriz(i,j);
    % endif

  endfor
endfor
i=n1;
for j=1:n1-1
  ev(i,j)=0;  # Definido como cero
  wv(i,j)=-dx*vwest(i,j)*vfw(i,j)-re1;
  nv(i,j)=dx*vsouth(i,j+1)*(1-vfs(i,j+1))-re1;
  sv(i,j)=-dx*vsouth(i,j)*vfs(i,j)-re1;
  av(i,j)=dx2/dt+2*re1-ev(i,j)-wv(i,j)-nv(i,j)-sv(i,j); # Tension de corte en la pared con primer orden
  bv(i,j)=2*Utop2*re1+dx2/dt*v0(i,j)+vsourcew(i,j)+vsources(i,j)-vsources(i,j+1);
  #Si estoy en CN, agrego el término b de la propuesta de la catedra en el pdf TP_Final
  % if metodo_temporal == "C"
  %   bv(i,j) = bv(i,j) + b_CNv_matriz(i,j);
  % endif

endfor

% if metodo_temporal == "C"
%   bv = bv + b_CNv_matriz;
% endif

#---------------------------------------------
# Calculo de los coeficientes para la ecuacion de presión
#---------------------------------------------

i=1;
j=1;
du(i,j)=dx/au(i,j);   
Usom(i,j)=(bu(i,j)-(eu(i,j)*uvel(i+1,j)+nu(i,j)*uvel(i,j+1)))/au(i,j);
dv(i,j)=dx/av(i,j);   
Vsom(i,j)=(bv(i,j)-(ev(i,j)*vvel(i+1,j)+nv(i,j)*vvel(i,j+1)))/av(i,j);
for j=2:n1-2
  du(i,j)=dx/au(i,j);   
  Usom(i,j)=(bu(i,j)-(eu(i,j)*uvel(i+1,j)+nu(i,j)*uvel(i,j+1)+su(i,j)*uvel(i,j-1)))/au(i,j);
  dv(i,j)=dx/av(i,j);   
  Vsom(i,j)=(bv(i,j)-(ev(i,j)*vvel(i+1,j)+nv(i,j)*vvel(i,j+1)+sv(i,j)*vvel(i,j-1)))/av(i,j);
endfor
j=n1-1;
du(i,j)=dx/au(i,j);   
Usom(i,j)=(bu(i,j)-(eu(i,j)*uvel(i+1,j)+nu(i,j)*uvel(i,j+1)+su(i,j)*uvel(i,j-1)))/au(i,j);
dv(i,j)=dx/av(i,j);   
Vsom(i,j)=(bv(i,j)-(ev(i,j)*vvel(i+1,j)+sv(i,j)*vvel(i,j-1)))/av(i,j);
j=n1;
du(i,j)=dx/au(i,j);   
Usom(i,j)=(bu(i,j)-(eu(i,j)*uvel(i+1,j)+su(i,j)*uvel(i,j-1)))/au(i,j);
dv(i,j)=0;   
Vsom(i,j)=0;

for i=2:n1-2
  j=1;
  du(i,j)=dx/au(i,j);   
  Usom(i,j)=(bu(i,j)-(eu(i,j)*uvel(i+1,j)+wu(i,j)*uvel(i-1,j)+nu(i,j)*uvel(i,j+1)))/au(i,j);
  dv(i,j)=dx/av(i,j);   
  Vsom(i,j)=(bv(i,j)-(ev(i,j)*vvel(i+1,j)+wv(i,j)*vvel(i-1,j)+nv(i,j)*vvel(i,j+1)))/av(i,j);
  for j=2:n1-2
    du(i,j)=dx/au(i,j);   
    Usom(i,j)=(bu(i,j)-(eu(i,j)*uvel(i+1,j)+wu(i,j)*uvel(i-1,j)+nu(i,j)*uvel(i,j+1)+su(i,j)*uvel(i,j-1)))/au(i,j);
    dv(i,j)=dx/av(i,j);   
    Vsom(i,j)=(bv(i,j)-(ev(i,j)*vvel(i+1,j)+wv(i,j)*vvel(i-1,j)+nv(i,j)*vvel(i,j+1)+sv(i,j)*vvel(i,j-1)))/av(i,j);
  endfor
  j=n1-1;
  du(i,j)=dx/au(i,j);   
  Usom(i,j)=(bu(i,j)-(eu(i,j)*uvel(i+1,j)+wu(i,j)*uvel(i-1,j)+nu(i,j)*uvel(i,j+1)+su(i,j)*uvel(i,j-1)))/au(i,j);
  dv(i,j)=dx/av(i,j);   
  Vsom(i,j)=(bv(i,j)-(ev(i,j)*vvel(i+1,j)+wv(i,j)*vvel(i-1,j)+sv(i,j)*vvel(i,j-1)))/av(i,j);
  j=n1;
  du(i,j)=dx/au(i,j);   
  Usom(i,j)=(bu(i,j)-(eu(i,j)*uvel(i+1,j)+wu(i,j)*uvel(i-1,j)+su(i,j)*uvel(i,j-1)))/au(i,j);
  dv(i,j)=0;   
  Vsom(i,j)=0;
endfor

i=n1-1;
j=1;
du(i,j)=dx/au(i,j);   
Usom(i,j)=(bu(i,j)-(wu(i,j)*uvel(i-1,j)+nu(i,j)*uvel(i,j+1)))/au(i,j);
dv(i,j)=dx/av(i,j);   
Vsom(i,j)=(bv(i,j)-(ev(i,j)*vvel(i+1,j)+wv(i,j)*vvel(i-1,j)+nv(i,j)*vvel(i,j+1)))/av(i,j);
for j=2:n1-2
  du(i,j)=dx/au(i,j);   
  Usom(i,j)=(bu(i,j)-(wu(i,j)*uvel(i-1,j)+nu(i,j)*uvel(i,j+1)+su(i,j)*uvel(i,j-1)))/au(i,j);
  dv(i,j)=dx/av(i,j);   
  Vsom(i,j)=(bv(i,j)-(ev(i,j)*vvel(i+1,j)+wv(i,j)*vvel(i-1,j)+nv(i,j)*vvel(i,j+1)+sv(i,j)*vvel(i,j-1)))/av(i,j);
endfor
j=n1-1;
du(i,j)=dx/au(i,j);   
Usom(i,j)=(bu(i,j)-(wu(i,j)*uvel(i-1,j)+nu(i,j)*uvel(i,j+1)+su(i,j)*uvel(i,j-1)))/au(i,j);
dv(i,j)=dx/av(i,j);   
Vsom(i,j)=(bv(i,j)-(ev(i,j)*vvel(i+1,j)+wv(i,j)*vvel(i-1,j)+sv(i,j)*vvel(i,j-1)))/av(i,j);
j=n1;
du(i,j)=dx/au(i,j);   
Usom(i,j)=(bu(i,j)-(wu(i,j)*uvel(i-1,j)+su(i,j)*uvel(i,j-1)))/au(i,j);
dv(i,j)=0;   
Vsom(i,j)=0;


i=n1;
j=1;
du(i,j)=0;   
Usom(i,j)=0;
dv(i,j)=dx/av(i,j);   
Vsom(i,j)=(bv(i,j)-(wv(i,j)*vvel(i-1,j)+nv(i,j)*vvel(i,j+1)))/av(i,j);
for j=2:n1-2
  du(i,j)=0;   
  Usom(i,j)=0;
  dv(i,j)=dx/av(i,j);   
  Vsom(i,j)=(bv(i,j)-(wv(i,j)*vvel(i-1,j)+nv(i,j)*vvel(i,j+1)+sv(i,j)*vvel(i,j-1)))/av(i,j);
endfor
j=n1-1;
du(i,j)=0;   
Usom(i,j)=0;
dv(i,j)=dx/av(i,j);   
Vsom(i,j)=(bv(i,j)-(wv(i,j)*vvel(i-1,j)+sv(i,j)*vvel(i,j-1)))/av(i,j);
j=n1;
du(i,j)=0;   
Usom(i,j)=0;
dv(i,j)=0;   
Vsom(i,j)=0;


#---------------------------------------------
# Matriz para la presion
#---------------------------------------------
mindex=1;
i=1;
j=1;
MP(mindex,mindex+1)=-dv(i,j);
MP(mindex,mindex+n1)=-du(i,j);
MP(mindex,mindex)=du(i,j)+dv(i,j);
bp(mindex)=-Usom(i,j)-Vsom(i,j);
mindex=mindex+1;
for j=2:n1-1
  MP(mindex,mindex+1)=-dv(i,j);
  MP(mindex,mindex-1)=-dv(i,j-1);
  MP(mindex,mindex+n1)=-du(i,j);
  MP(mindex,mindex)=du(i,j)+dv(i,j)+dv(i,j-1);
  bp(mindex)=-Usom(i,j)+Vsom(i,j-1)-Vsom(i,j);
  mindex=mindex+1;
endfor
j=n1;
MP(mindex,mindex-1)=-dv(i,j-1);
MP(mindex,mindex+n1)=-du(i,j);
MP(mindex,mindex)=du(i,j)+dv(i,j-1);
bp(mindex)=-Usom(i,j)+Vsom(i,j-1);
mindex=mindex+1;

for i=2:n1-1
  j=1;
  MP(mindex,mindex+1)=-dv(i,j);
  MP(mindex,mindex+n1)=-du(i,j);
  MP(mindex,mindex-n1)=-du(i-1,j);
  MP(mindex,mindex)=du(i,j)+du(i-1,j)+dv(i,j);
  bp(mindex)= Usom(i-1,j)-Usom(i,j)-Vsom(i,j);
  mindex=mindex+1;
  for j=2:n1-1
    MP(mindex,mindex+1)=-dv(i,j);
    MP(mindex,mindex-1)=-dv(i,j-1);
    MP(mindex,mindex+n1)=-du(i,j);
    MP(mindex,mindex-n1)=-du(i-1,j);
    MP(mindex,mindex)=du(i,j)+du(i-1,j)+dv(i,j)+dv(i,j-1);
    bp(mindex)= Usom(i-1,j)-Usom(i,j)+Vsom(i,j-1)-Vsom(i,j);
    mindex=mindex+1;
  endfor
  j=n1;
  MP(mindex,mindex-1)=-dv(i,j-1);
  MP(mindex,mindex+n1)=-du(i,j);
  MP(mindex,mindex-n1)=-du(i-1,j);
  MP(mindex,mindex)=du(i,j)+du(i-1,j)+dv(i,j-1);
  bp(mindex)= Usom(i-1,j)-Usom(i,j)+Vsom(i,j-1);
  mindex=mindex+1;
endfor

i=n1;
j=1;
MP(mindex,mindex+1)=-dv(i,j);
MP(mindex,mindex-n1)=-du(i-1,j);
MP(mindex,mindex)=du(i-1,j)+dv(i,j);
bp(mindex)= Usom(i-1,j)-Vsom(i,j);
mindex=mindex+1;
for j=2:n1-1
  MP(mindex,mindex+1)=-dv(i,j);
  MP(mindex,mindex-1)=-dv(i,j-1);
  MP(mindex,mindex-n1)=-du(i-1,j);
  MP(mindex,mindex)=du(i-1,j)+dv(i,j)+dv(i,j-1);
  bp(mindex)= Usom(i-1,j)+Vsom(i,j-1)-Vsom(i,j);
  mindex=mindex+1;
endfor
# Forzamos el volumen final para que tenga presion cero (eliminamos la singularidad de la matrix)
MP(mindex,mindex)=1;
bp(mindex)=0;

# Utilizamos la matriz rala para calcular
MPS=sparse(MP);

#-------------------------------------------------------
# Resolvemos la presion (MPS * pres = bp) obtenemos el vector de presiones pres
#-------------------------------------------------------

pres=MPS\bp;

mindex=1;
for i=1:n1
  for j=1:n1
    presion(i,j)=pres(mindex);
    mindex=mindex+1;
  endfor
endfor

% if metodo_temporal == "C"
%   bu = bu + b_CNu_matriz;
%   bv = bv + b_CNv_matriz;
% endif


# -------------------------------------------------------------------
#
#   PASO 2 SIMPLER: RESUELVO ECUACIONES DE MOMENTO DADO EL CAMPO DE PRESIONES - LAS VELOCIDADES OBTENIDAS NO SON DIVERGENCIA LIBRE
#
# -------------------------------------------------------------------

#-------------------------------------------------------
# Resolvemos la ecuacion de momento en x
#-------------------------------------------------------
mindex=1;

i=1;
j=1;
MMx(mindex,mindex+1)=nu(i,j);
MMx(mindex,mindex+n1)=eu(i,j);
MMx(mindex,mindex)=au(i,j);
bux(mindex)= bu(i,j)+dx*(pres(mindex)-pres(mindex+n1));
mindex=mindex+1;
for j=2:n1-1
  MMx(mindex,mindex+1)=nu(i,j);
  MMx(mindex,mindex-1)=su(i,j);
  MMx(mindex,mindex+n1)=eu(i,j);
  MMx(mindex,mindex)=au(i,j);
  bux(mindex)= bu(i,j)+dx*(pres(mindex)-pres(mindex+n1));
  mindex=mindex+1;
endfor
j=n1;
MMx(mindex,mindex-1)=su(i,j);
MMx(mindex,mindex+n1)=eu(i,j);
MMx(mindex,mindex)=au(i,j);
bux(mindex)= bu(i,j)+dx*(pres(mindex)-pres(mindex+n1));
mindex=mindex+1;

for i=2:n1-2
  j=1;
  MMx(mindex,mindex+1)=nu(i,j);
  MMx(mindex,mindex+n1)=eu(i,j);
  MMx(mindex,mindex-n1)=wu(i,j);
  MMx(mindex,mindex)=au(i,j);
  bux(mindex)= bu(i,j)+dx*(pres(mindex)-pres(mindex+n1));
  mindex=mindex+1;
  for j=2:n1-1
    MMx(mindex,mindex+1)=nu(i,j);
    MMx(mindex,mindex-1)=su(i,j);
    MMx(mindex,mindex+n1)=eu(i,j);
    MMx(mindex,mindex-n1)=wu(i,j);
    MMx(mindex,mindex)=au(i,j);
    bux(mindex)= bu(i,j)+dx*(pres(mindex)-pres(mindex+n1));
    mindex=mindex+1;
  endfor
  j=n1;
  MMx(mindex,mindex-1)=su(i,j);
  MMx(mindex,mindex+n1)=eu(i,j);
  MMx(mindex,mindex-n1)=wu(i,j);
  MMx(mindex,mindex)=au(i,j);
  bux(mindex)= bu(i,j)+dx*(pres(mindex)-pres(mindex+n1));
  mindex=mindex+1;
endfor

i=n1-1;
j=1;
MMx(mindex,mindex+1)=nu(i,j);
MMx(mindex,mindex-n1)=wu(i,j);
MMx(mindex,mindex)=au(i,j);
bux(mindex)= bu(i,j)+dx*(pres(mindex)-pres(mindex+n1));
mindex=mindex+1;
for j=2:n1-1
  MMx(mindex,mindex+1)=nu(i,j);
  MMx(mindex,mindex-1)=su(i,j);
  MMx(mindex,mindex-n1)=wu(i,j);
  MMx(mindex,mindex)=au(i,j);
  bux(mindex)= bu(i,j)+dx*(pres(mindex)-pres(mindex+n1));
  mindex=mindex+1;
endfor
j=n1;
MMx(mindex,mindex-1)=su(i,j);
MMx(mindex,mindex-n1)=wu(i,j);
MMx(mindex,mindex)=au(i,j);
bux(mindex)= bu(i,j)+dx*(pres(mindex)-pres(mindex+n1));
mindex=mindex+1;

MMxS=sparse(MMx);
u=MMxS\bux;

mindex=1;
for i=1:n1-1
  for j=1:n1
    uvel(i,j)=u(mindex);
    mindex=mindex+1;
  endfor
endfor

#-------------------------------------------------------
# Resolvemos la ecuacion de momento en y
#-------------------------------------------------------
mindex=1;

i=1;
j=1;
MMy(mindex,mindex+1)=nv(i,j);
MMy(mindex,mindex+n1-1)=ev(i,j);
MMy(mindex,mindex)=av(i,j);
buy(mindex)= bv(i,j)+dx*(pres((i-1)*n1+j)-pres((i-1)*n1+j+1));
mindex=mindex+1;
for j=2:n1-2
  MMy(mindex,mindex+1)=nv(i,j);
  MMy(mindex,mindex-1)=sv(i,j);
  MMy(mindex,mindex+n1-1)=ev(i,j);
  MMy(mindex,mindex)=av(i,j);
  buy(mindex)= bv(i,j)+dx*(pres((i-1)*n1+j)-pres((i-1)*n1+j+1));
  mindex=mindex+1;
endfor
j=n1-1;
MMy(mindex,mindex-1)=sv(i,j);
MMy(mindex,mindex+n1-1)=ev(i,j);
MMy(mindex,mindex)=av(i,j);
buy(mindex)= bv(i,j)+dx*(pres((i-1)*n1+j)-pres((i-1)*n1+j+1));
mindex=mindex+1;



for i=2:n1-1
  j=1;
  MMy(mindex,mindex+1)=nv(i,j);
  MMy(mindex,mindex+n1-1)=ev(i,j);
  MMy(mindex,mindex-n1+1)=wv(i,j);
  MMy(mindex,mindex)=av(i,j);
  buy(mindex)= bv(i,j)+dx*(pres((i-1)*n1+j)-pres((i-1)*n1+j+1));
  mindex=mindex+1;
  for j=2:n1-2
    MMy(mindex,mindex+1)=nv(i,j);
    MMy(mindex,mindex-1)=sv(i,j);
    MMy(mindex,mindex+n1-1)=ev(i,j);
    MMy(mindex,mindex-n1+1)=wv(i,j);
    MMy(mindex,mindex)=av(i,j);
    buy(mindex)= bv(i,j)+dx*(pres((i-1)*n1+j)-pres((i-1)*n1+j+1));
    mindex=mindex+1;
  endfor
  j=n1-1;
  MMy(mindex,mindex-1)=sv(i,j);
  MMy(mindex,mindex+n1-1)=ev(i,j);
  MMy(mindex,mindex-n1+1)=wv(i,j);
  MMy(mindex,mindex)=av(i,j);
  buy(mindex)= bv(i,j)+dx*(pres((i-1)*n1+j)-pres((i-1)*n1+j+1));
  mindex=mindex+1;
endfor

i=n1;
j=1;
MMy(mindex,mindex+1)=nv(i,j);
MMy(mindex,mindex-n1+1)=wv(i,j);
MMy(mindex,mindex)=av(i,j);
buy(mindex)= bv(i,j)+dx*(pres((i-1)*n1+j)-pres((i-1)*n1+j+1));
mindex=mindex+1;
for j=2:n1-2
  MMy(mindex,mindex+1)=nv(i,j);
  MMy(mindex,mindex-1)=sv(i,j);
  MMy(mindex,mindex-n1+1)=wv(i,j);
  MMy(mindex,mindex)=av(i,j);
  buy(mindex)= bv(i,j)+dx*(pres((i-1)*n1+j)-pres((i-1)*n1+j+1));
  mindex=mindex+1;
endfor
j=n1-1;
MMy(mindex,mindex-1)=sv(i,j);
MMy(mindex,mindex-n1+1)=wv(i,j);
MMy(mindex,mindex)=av(i,j);
buy(mindex)= bv(i,j)+dx*(pres((i-1)*n1+j)-pres((i-1)*n1+j+1));

MMyS=sparse(MMy);
v=MMyS\buy;

mindex=1;
for i=1:n1
  for j=1:n1-1
    vvel(i,j)=v(mindex);
    mindex=mindex+1;
  endfor
endfor

# -------------------------------------------------------------------
#
#   PASO 3 SIMPLER: ENCONTRAMOS LA CORRECION DE LA PRESION QUE HACE EL CAMPO DE VELOCIDADES DIVERGENCIA LIBRE
#
# -------------------------------------------------------------------

#-------------------------------------------------------
# Resolvemos para la correccion en la presion
#-------------------------------------------------------

mindex=1;
i=1;
j=1;
bp(mindex)=-uvel(i,j)-vvel(i,j);
mindex=mindex+1;
for j=2:n1-1
  bp(mindex)=-uvel(i,j)+vvel(i,j-1)-vvel(i,j);
  mindex=mindex+1;
endfor
j=n1;
bp(mindex)=-uvel(i,j)+vvel(i,j-1);
mindex=mindex+1;

for i=2:n1-1
  j=1;
  bp(mindex)= uvel(i-1,j)-uvel(i,j)-vvel(i,j);
  mindex=mindex+1;
  for j=2:n1-1
    bp(mindex)= uvel(i-1,j)-uvel(i,j)+vvel(i,j-1)-vvel(i,j);
    mindex=mindex+1;
  endfor
  j=n1;
  bp(mindex)= uvel(i-1,j)-uvel(i,j)+vvel(i,j-1);
  mindex=mindex+1;
endfor

i=n1;
j=1;
bp(mindex)= uvel(i-1,j)-vvel(i,j);
mindex=mindex+1;
for j=2:n1-1
  bp(mindex)= uvel(i-1,j)+vvel(i,j-1)-vvel(i,j);
  mindex=mindex+1;
endfor
# Forzamos el volumen final para que tenga presion cero (eliminamos la singularidad de la matrix)
bp(mindex)=0;

#-------------------------------------------------------
# Resolvemos la presion (MPS * pres = bp) obtenemos el vector de presiones pres
#-------------------------------------------------------

prescor=MPS\bp;

# -------------------------------------------------------------------
#
#   PASO 4 SIMPLER: CORREGIMOS EL CAMPO DE VELOCIDADES CON LA CORRECION DE LA PRESION
#
# -------------------------------------------------------------------

#-------------------------------------------------------
# Corregimos la velocidad para que sea divergencia libre
#-------------------------------------------------------
mindex=1;
for i=1:n1-1
  for j=1:n1-1
    uvel(i,j)=uvel(i,j)+du(i,j)*(prescor(mindex)-prescor(mindex+n1));
    vvel(i,j)=vvel(i,j)+dv(i,j)*(prescor(mindex)-prescor(mindex+1));
    mindex=mindex+1;
  endfor
  j=n1;
  uvel(i,j)=uvel(i,j)+du(i,j)*(prescor(mindex)-prescor(mindex+n1));
  mindex=mindex+1;
endfor

i=n1;
for j=1:n1-1
  vvel(i,j)=vvel(i,j)+dv(i,j)*(prescor(mindex)-prescor(mindex+1));
  mindex=mindex+1;
endfor


# -------------------------------------------------------------------
#
#
endfor       #termina iteracion simpler, indice lsimpler
#
#
# -------------------------------------------------------------------


# -------------------------------------------------------------------
#
#
#Al final del loop, actualizo la matriz b_CN siempre y cuando quiera aplicar CN y no EI. Para EI tengo que dejarla nula
#Asumo que los objetos con la info. de interés son uvel, u0, vvel y v0
#
# -------------------------------------------------------------------
if metodo_temporal == "C"
  % b_CNu_matriz = zeros(n1-1,n1); #matriz actual (n1,n1)
  % b_CNu_vector =zeros(n2-n1,1); # misma matriz en forma vectorial (n2)
  % b_CNu0_matriz = zeros(n1-1,n1); #matriz del paso anterior (n1,n1)
  % b_CNv_matriz = zeros(n1,n1-1); # matriz actual (n1,n1)
  % b_CNv_vector =zeros(n2-n1,1); # misma matriz en forma vectorial (n2)
  % b_CNv0_matriz = zeros(n1,n1-1); # matriz del paso anterior (n1,n1)

  % b_u = -b_u + (uvel-u0)/dt;
  % b_v = -b_v + (vvel-v0)/dt;
  b_CNu_matriz = - b_CNu_matriz + (uvel-u0)/dt;
  b_CNv_matriz = - b_CNv_matriz + (vvel-v0)/dt;

  % b_CNu_matriz =  - b_CNu0_matriz + (uvel - u0)/dt*dx2;
  % b_CNv_matriz =  - b_CNv0_matriz + (vvel - v0)/dt*dx2;
  % b_CNu_matriz =  - b_CNu0_matriz + bu;
  % b_CNv_matriz =  - b_CNv0_matriz + bv;


  #Redefino b_CNu0_matriz y b_CNv0_matriz
  % b_CNu0_matriz = b_CNu_matriz;
  % b_CNv0_matriz = b_CNv_matriz;

  #2*Utop2*re1 + dx2/dt*v0(i,j) +vsourcew(i,j)+vsources(i,j)-vsources(i,j+1);
endif

# -------------------------------------------------------------------
#
#   PASO 5 SIMPLER: TENEMOS EL NUEVO CAMPO DE VELOCIDADES EN n+1, verificamos
#
# -------------------------------------------------------------------

#-------------------------------------------------------
# Verificacion de la evolucion temporal de la solucion numerica
#-------------------------------------------------------

#---------- conservacion de la masa
masa(k)=0;
i=1;
j=1;
masa(k)=masa(k)+uvel(i,j)+vvel(i,j);
for j=2:n1-1
  masa(k)=masa(k)+uvel(i,j)+vvel(i,j)-vvel(i,j-1);
endfor
j=n1;
masa(k)=masa(k)+uvel(i,j)-vvel(i,j-1);
for i=2:n1-1
  j=1;
  masa(k)=masa(k)+uvel(i,j)-uvel(i-1,j)+vvel(i,j);
  for j=2:n1-1
    masa(k)=masa(k)+uvel(i,j)-uvel(i-1,j)+vvel(i,j)-vvel(i,j-1);
  endfor
  j=n1;
  masa(k)=masa(k)+uvel(i,j)-uvel(i-1,j)-vvel(i,j-1);
endfor
i=n1;
j=1;
masa(k)=masa(k)-uvel(i-1,j)+vvel(i,j);
for j=2:n1-1
  masa(k)=masa(k)-uvel(i-1,j)+vvel(i,j)-vvel(i,j-1);
endfor
j=n1;
masa(k)=masa(k)-uvel(i-1,j)-vvel(i,j-1);

#---------- variacion temporal del campo de velocidades = ver estado estacionario
derivadatuvel(k)=norm(uvel-u0)/dt;
derivadatvvel(k)=norm(vvel-v0)/dt;


#---------- variacion temporal = ver evaluacion temporal de alguna variable
uvel1(k)=uvel(i1,j1);
vvel1(k)=vvel(i1,j1);
pre1(k)=presion(i1,j1);



#Calculo las velocidades en el centro a cada tiempo
% u(0.5, i) = (uvel(int8(n1/2),i)+uvel(int8((n1-1)/2),i))/2

% v(i, 0.5) = (vvel(i,int8(n1/2))+vvel(i,int8((n1-1)/2)))/2

posicion_centro_izq = n1/2;
posicion_centro_der = n1/2 + 1;

iii = posicion_centro_izq;
u_05_izq = (uvel(int8(n1/2),iii)+uvel(int8((n1-1)/2),iii))/2;
iii = posicion_centro_der;
u_05_der = (uvel(int8(n1/2),iii)+uvel(int8((n1-1)/2),iii))/2;

u_05_05 = (u_05_izq + u_05_der)/2;

iii = posicion_centro_izq;
v_izq_05 = (vvel(iii,int8(n1/2))+vvel(iii,int8((n1-1)/2)))/2;
iii = posicion_centro_der;
v_der_05 = (vvel(iii,int8(n1/2))+vvel(iii,int8((n1-1)/2)))/2;
v_05_05 = (v_izq_05 + v_der_05)/2;

#termino el cálculo de las velocidades en el centro

tiempo =tiempo +dt;

#---------------------------------------------
#--------- imprimo en archivo variables para controlar simulacion
fprintf(fid, "%e %e %e %e %e %e\n", tiempo, masa(k),derivadatuvel(k),derivadatvvel(k), u_05_05, v_05_05);

#-------------------------------------------------------
# Evalúo si se llegó al estado estacionario
#-------------------------------------------------------
if condicion_limite_temporal(derivadatuvel(k), derivadatvvel(k))
  strcat("Se llego al estado estacionario con Ndeltat =  ", num2str(k))
  Ndeltat = k; #Test para ver si se guarda el nro de pasos de tiempo ejecutados
  break;
endif

if k == Ndeltat - 1
  strcat("Se llego al nro maximo de pasos de tiempo Ndeltat =  ", num2str(k))
endif

#---------------------------------------------
#---------------------------------------------
#
#  FINALIZA LOOP DE ITERATION TEMPORAL  
#
#
endfor
#
#
#---------------------------------------------
#---------------------------------------------

t2=time();
"Tiempo de Calculo en segundos para";
Ndeltat;
"iteraciones";
t2-t1;

#---------------------------------------------
#--------- cierro archivo
fclose(fid);


#---------------------------------------------
#
# Archivo de parámetros
#
#---------------------------------------------
#abro archivo
fid = fopen(archivo_parametros, "w");
fprintf(fid,"%e %e %e %e %e\n", n1, Re, Ndeltat, dt, t2-t1);
#cierro archivo
fclose(fid);
#---------------------------------------------
#
#              REALIZACION DE GRAFICOS
#
#---------------------------------------------

#---------------------------------------------
#--------- abro archivo para imprimir velocidades en el centro - Nota: por conservacion masa la integral de cada una de estas debe ser = 0
#---------  Ideal si
fid = fopen(archivo_velocidades_centrales, "w");
xtot=dx/2;
for i=1:n1
  fprintf(fid, "%e %e %e\n", xtot, (uvel(int8(n1/2),i)+uvel(int8((n1-1)/2),i))/2,(vvel(i,int8(n1/2))+vvel(i,int8((n1-1)/2)))/2);
  xtot=xtot+dx;
endfor
fclose(fid);

#---------- verificacion simulacion
% plot(masa); pause(3);
% semilogy(derivadatuvel);pause(3);
% semilogy(derivadatvvel);pause(3);
% plot(uvel1);pause(3);
% plot(vvel1);pause(3);
% plot(pre1);pause(3);


#---------- ploteo de resultados

#---------- grafico 3D de velocidad u
% txu=linspace(0+1/n1,1-1/n1,n1-1);
% tyu=linspace(0+1/n1/2,1-1/n1/2,n1);
% [xx, yy] = meshgrid (txu, tyu);
% mesh(xx,yy,uvel');pause(10);

#---------- grafico 3D de velocidad v
% mesh(xx,yy,vvel);pause(10);

#---------- grafico 3D de presion
% txp=linspace(0+1/n1/2,1-1/n1/2,n1);
% typ=linspace(0+1/n1/2,1-1/n1/2,n1);
% [xx, yy] = meshgrid (txp, typ);
% mesh(xx,yy,presion);pause(10);

#---------- 2D coloreado
% contour(txp,typ,presion);pause(10);

#---------- lineas de velocidad 
% plot(uvel(int8((n1-1)/2),:));pause(10);
% plot(vvel(:,int8((n1-1)/2)));pause(10);
% plot(tyu,vvel(:,int8((n1-1)/2)),tyu,uvel(int8((n1-1)/2),:));pause(10);

#---------- Vector velocidad - simil streamlines
i=1;
j=1;
usl(i,j)=(uvel(i,j))/2;
vsl(i,j)=(vvel(i,j))/2;
for j=2:n1-1
  usl(i,j)=(uvel(i,j))/2;
  vsl(i,j)=(vvel(i,j)+vvel(i,j-1))/2;
endfor
j=n1;
usl(i,j)=(uvel(i,j))/2;
vsl(i,j)=(vvel(i,j-1))/2;
for i=2:n1-1
  j=1;
  usl(i,j)=(uvel(i,j)+uvel(i-1,j))/2;
  vsl(i,j)=(vvel(i,j))/2;
  for j=2:n1-1
    usl(i,j)=(uvel(i,j)+uvel(i-1,j))/2;
    vsl(i,j)=(vvel(i,j)+vvel(i,j-1))/2;
  endfor
  j=n1;
  usl(i,j)=(uvel(i,j)+uvel(i-1,j))/2;
  vsl(i,j)=(vvel(i,j-1))/2;
endfor
i=n1;
j=1;
usl(i,j)=(uvel(i-1,j))/2;
vsl(i,j)=(vvel(i,j))/2;
for j=2:n1-1
  usl(i,j)=(uvel(i-1,j))/2;
  vsl(i,j)=(vvel(i,j)+vvel(i,j-1))/2;
endfor
j=n1;
usl(i,j)=(uvel(i-1,j))/2;
vsl(i,j)=(vvel(i,j-1))/2;
# COmando para graficar velocidad proporcional a su módulo
% quiver(usl',vsl',1);
for i=1:n1
for j=1:n1
  vsln(i,j)=vsl(i,j)/(usl(i,j)^2+vsl(i,j)^2)^.5;
  usln(i,j)=usl(i,j)/(usl(i,j)^2+vsl(i,j)^2)^.5;
endfor
endfor
# COmando para graficar velocidad normalizada a Módulo 1 (.4 es la longitud de la flecha en el grafico)
#quiver(usln',vsln',.4);

#---------- Vorticidad - wz = dv/dy-du/dx 
for i=1:n1-1
  for j=1:n1-1
    wz(i,j)=(-uvel(i,j+1)+uvel(i,j)+vvel(i+1,j)-vvel(i,j))/dx;
  endfor
endfor
#---------- grafico de vorticidad
txwz=linspace(dx,1-dx,n1-1);
tywz=linspace(dx,1-dx,n1-1);
#---------- 2D coloreado
#contour(txwz,tywz,wz',-20:.2:20);


endfunction