Comentarios varios


Puede ser que no llegue al mismo comportamiento que en la tabla para los errores de fase y amplitud. Eso se puede deber a que en la tabla valen las ecuaciones linealizadas. Debería llegarse al mismo resultado cuando tit-> 0, ej 10^-4

Las tolerancias que se eligen en bisección y en la determinación de Tau tienen que ser tales que cuando plotee los errores finales tienen que ser mayores a esos errores

Se podría resolver el problema sin el método de bisección. En Moodle están los códigos de los métodos numéricos. En ese caso actúan iterativamente:
y^{n+1} = y^n - sin(y^{n+1})Deltat
y^* = y^n
do while (Error < tol y k<10000)
	y^{n+1} = y^n - sin(y^*)
	Error = |y^{n+1} - y^*|
	y^* = y^{n+1}
enddo
Esto converge más rápido que bisección

¿Qué ventaja tiene elegir 2CN+LeapFrog ?
Ver el error de fase de cada método


El error de fase me da del orden de 1. Puede ser que esté haciendo mal la cuenta o estoy eligiendo mal el tiempo tal que se cumple un período.
Solución:
*Ya que ejecuto hasta 2 veces el período, elijo calcular el error de fase hasta la mitad de los pasos
*Ejecuto hasta una vez el período y calculo el error de fase en el último valor


En el gráfico de comet del péndulo doble NO HAY que hacer un video. Estoy haciendo mal esto. Debería hacer la trayectoria graficada completamente (plot). Al hacer comet aparece una linda animación nomás, pero realmente tengo que hacer un plot de la trayectoria
