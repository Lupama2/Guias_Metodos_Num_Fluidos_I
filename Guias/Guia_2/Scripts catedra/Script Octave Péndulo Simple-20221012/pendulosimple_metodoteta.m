# Pendulo Simple Metodo Tita

clear uee;
clear uei;
clear ucn;
clear uexacta;

#Parámettros Fisicos y datos
w=1; #Frecuencia 
cv=0; #Coeficiente friccion
tita0=pi/4;
titapunto0=0;

#### PASOS DE TIEMPO
n=100;
####

#Tiempos de la Solucion
Tini=0;
nperiodos=2;
Tfin=nperiodos*2*pi/w; 
t=linspace(Tini, Tfin, n)'; # n=pasos de tiempo
h=(Tfin-Tini)/n;

# Solucion exacta
uexacta=zeros(2,n);
for i=1:n
   uexacta(1,i)=tita0*cos(w*h*i)+titapunto0/w*sin(w*h*i);
   uexacta(2,i)=w*tita0*sin(w*h*i)-titapunto0*cos(w*h*i);
endfor

#Definicion Matrices cada Metodo
# Euler Explicito
AEE=[1 h; -h*w^2  1];
# Euler Implicito
BEI=[1 -h; h*w^2  1];
# Tita=1/2 - Crank Nicholson
ACN=[1 h/2; -h/2*w^2  1];
BCN=[1 -h/2; h/2*w^2  1];
uee=zeros(2,n);
uei=zeros(2,n);
ucn=zeros(2,n);

# Solucion numérica metodos EE, EI y tita=1/2
cinicial=[tita0 titapunto0]';
uee(:,1)=AEE*cinicial;
uei(:,1)=BEI\cinicial;
ucn(:,1)=BCN\ACN*cinicial;
for i = 2:n
   uee(:,i)=AEE*uee(:,i-1);
   uei(:,i)=BEI\uei(:,i-1);
   ucn(:,i)=BCN\ACN*ucn(:,i-1);
endfor

# Ploteo de la solucion Trayectoria
% title ("Trayectoria EI");
% hold on;
% pause(20)
% comet(sin(uei(1,:)),-cos(uei(1,:)));
% hold on;
% title ("Trayectoria Crank-Nicolson");
% comet(sin(ucn(1,:)),-cos(ucn(1,:)));
% hold on;
% title ("Trayectoria EE");
% comet(sin(uee(1,:)),-cos(uee(1,:)));

% plot(t,uexacta(1,:),t,uee(1,:),t,uei(1,:),t,ucn(1,:))
% legend("Exacto","EE","EI","CN");
% ylabel ("Theta");
% xlabel ("tiempo");
% print -djpg foo.jpg

# Ploteo trayectoria tipo video
title ("Trayectorias");
#pause(20)
for l=1:n
plot(sin(uexacta(1,l)),-cos(uexacta(1,l)),sin(uee(1,l)),-cos(uee(1,l)),sin(uei(1,l)),-cos(uei(1,l)),sin(ucn(1,l)),-cos(ucn(1,l)),'linewidth',1.5);
legend("Exacto","EE","EI","CN");
axis([-1.2 1.2 -1 -0.5]); 
pause(0.1);
endfor

#Error en funcion de h
#n=2;
#for j=1:12
#n=n*2;
#inter(j)=1/n;
#pendulosimple_metodoteta;
#erroree(j)=norm(uexacta(:,n)-uee(:,n));
#errorei(j)=norm(uexacta(:,n)-uei(:,n));
#errorcn(j)=norm(uexacta(:,n)-ucn(:,n));
#endfor
#loglog(inter,erroree,inter,errorei,inter,errorcn,inter,y1,'linewidth',2,inter,y2,'linewidth', 2)
#legend("EE","EI","CN","1","2")
#y1=inter.^1
#y2=inter.^2


