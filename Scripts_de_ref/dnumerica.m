fre = 5;
#---------------------------------------------
#Datos para la discretizacion
#---------------------------------------------
nx = 19; #Numero de nodos internos
dx = 1 / (nx+1); #Distancia entre nodos

#---------------------------------------------
#Datos de la funcion en los nodos internos - f(x)=sin(fre*pi*x), donde fre es la frecuencia
#---------------------------------------------
for j=1:nx+1
  f(j)=sin(j*dx*fre*pi());
endfor

for j=3:nx-2
  fpexacta(j-2)=fre*pi()*cos(j*dx*fre*pi());	
endfor

#---------------------------------------------
# f'calculada con DC-2
#---------------------------------------------
for j=3:nx-2
  x(j-2)=dx*j; 
  fdc2(j-2)=(f(j+1)-f(j-1))/2/dx;
endfor

#---------------------------------------------
# f'calculada con DC-4
#---------------------------------------------
for j=3:nx-2
  fdc4(j-2)=(-f(j+2)+8*f(j+1)-8*f(j-1)+f(j-2))/12/dx;
endfor

#---------------------------------------------
# f'calculada con pade
#---------------------------------------------
A=zeros(nx-4,nx-4);
j=3;
#A(j-2,j-2)=1;
#A(j-2,j-1)=2;
#B(j-2)=(-5/2*f(j)+2*f(j+1)+1/2*f(j+2))/dx;
A(j-2,j-2)=1;
B(j-2)=(-f(j+2)+8*f(j+1)-8*f(j-1)+f(j-2))/12/dx;
for j=4:nx-3
  A(j-2,j-3)=1;
  A(j-2,j-2)=4;
  A(j-2,j-1)=1;
  B(j-2)=(3*(f(j+1)-f(j-1)))/dx;
endfor
j=nx-2;
#A(j-2,j-2)=1;
#A(j-2,j-3)=2;
#B(j-2)=(5/2*f(j)-2*f(j-1)-1/2*f(j-2))/dx;
A(j-2,j-2)=1;
B(j-2)=(-f(j+2)+8*f(j+1)-8*f(j-1)+f(j-2))/12/dx;
fpade=A\B';

plot(x,fpexacta,".k", "markersize", 10,x,fdc2,"linewidth", 2,";DC-2;",x,fdc4,";DC-4;","linewidth", 2,x,fpade,";PADE;","linewidth", 2)

pause(10)