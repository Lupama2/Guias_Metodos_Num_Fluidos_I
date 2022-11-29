



* Dudas:

* 1. Identifiqué bien cómo hacer solo un UP1, un DC2 y un QUICK? 
#RTA: Esquemáticamente, sí

* 2. En QUICK solo tengo que modificar el loop que dice "#Interior desde i=2 hasta n1-1 y desde j=2 hasta n1:", no?
#RTA: Sólo cambiar el for del interior. Sí.

* 3. En ese loop viejo se multiplica por dx*uwest(i,j) para directamente calcular un término en el cálculo de la integral 2 (ver mis notas)? O sea que para implementar QUICK eso no lo tengo que cambiar, no?
#RTA: El factor dx*uwest(i,j) es del término advectivo. Lo tengo que dejar al aplicar QUICK

* 4. Está bien dónde cambié la velocidad Utop1? Está cerca de la línea 620
#RTA: No. Modifiqué mal la función de velocidad. No depende de x sino de t!. Entonces tengo que  cambiar la constante al inicio de cada paso de tiempo y listo


* 5. Cuando ejecuto el código aparece el resultado
* ans = 39.309.
* Pero no siempre aparece el mismo valor. El código no es determinista?
#RTA: Este resultado es el tiempo de ejecución del programa. El valor de ans (14.796 o 42.115) da distinto si conecto o no la notebook a lacorriente. Es el tiempo que le toma al programa hacer la cuenta

* 6. Euler implícito solo se usa en la integral 1, no?
#RTA: Sí. Hay que usar la ayuda. Hay que aplicarlo afuera del loop de simpler al final del loop de tiempo

* 7. No entiendo cómo aplicar CN a partir de EI. Por qué tengo que usar la mitad del paso de tiempo?
#RTA: Después me di cuenta cómo hacerlo

* 2. La cuenta para CN la puse cuando terminó la iteración simples, está bien?
#RTA: Sí.

* * 9. La matriz b de CN solo afecta a las matrices bu y bv en el PASO 1 SIMPLER, no?
#RTA: Creo que sí. No lo pregunté para después verlo yo solo


#LO ESTOY AMPLICANDO MAL! Hay que tener cuidado al agregar C-N. Porque el problema no es y' = f, sino y' = f \Deltax^2 por la integral espacial. Hay que dijarse que cuando sumemos las velocidades tengan las mismas unidades

* 10. Están bien los resultados que encontré?
#RTA: en base a lo que me dijo tengo que reveer los resultados

* 11. Puedo usar números impares para n1 de modo de no tener que interpolar el valor en 0.5?
#RTA: Sí o sí voy a tener que hacer una interpolación porque las grillas de u y v están desplazadas. La aproximación (interpolación) es de segundo orden, del mismo que el esquema numérico para el dominio espacial.

* 12. En el inciso c me diverge usando UP1, puede ser que lo haya implementado mal? O será el valor de Deltat?
#RTA: entiendo que está bien implementado, tengo que disminuir dt.

* 13. Está bien el algoritmo que apliqué en el inciso e?
* Tengo que ir cambiando dt?
# RTA: e. Variar dt tmb. Como es estado estacionario no depende de dt. Entonces, tengo que encontrar el n1 mínimo que cumpla la condición. Luego, para ese n1 tengo que encontrar dt. Ninguna corrida con dt y n1 particular me tiene que tomar más de una hora. No es necesario llegar al n1 exacto, se busca responder las siguientes preguntas: ¿Qué pasa si aumento o disminuyo el dt? ¿Cómo afecta al costo computacional? ¿Qué sentido tiene? Podés usar un transitorio no físico y sin sentido para llegar al estado estacionario

* 14. Me tardar mucho tiempo en correr el programa. ¿Qué se puede hacer?
# RTA: Como estoy en el estado estacionario, la solución no debería depender de t ni de dt (paso de tiempo). Así que puedo usar cualquiera! Como dt está adimensionalizado entonces puedo usar dt alto a medida que subo Re.

* 15. ¿Qué significa el warning "matrix singular to machine precision"?
# RTA: El error indica que está resolviendo un sistema lineal y la matriz que resuelve está cerca de ser singular.

#No sé a qué hace referencia el siguiente comentario
#En algún momento tengo que sumar su efecto en la ecuación diferencial.

Tuve que eliminar el clear del inicio del código. No sé en qué afectará
Terminar de programar QUICK


DUDAS:


Dejé corriendo el inciso b 

* El e me toma mucho tiempo porque tengo que evolucionar hasta el estado final para evaluar el error. Hay otra forma de estimarlo?
* Está mal poner una tolerancia de 1e-3 en lugar de 1e-6 para que corra más rápido?
No está bueno. Pero está bien. Para algunos casos puede no tener sentido

* ¿Qué puedo hacer para acelerar la convergencia? Dejé el b y c corriendo toda la noche y no terminaron aún. Comentar cómo estimo el dt.
* u ya está en la posición vertical 0.5 y cuando extraigo los valores, el correspondiente a 0.5 es u(0.5,0.5)? (Al revés con v?)
Sí. Esto vale siempre que el n1 sea PAR
Para n1 impar NO QUEDA u() en x = 0.5

* d. No sé en qué zona calcular el error. No me queda lineal

Creo que en el inciso e no estoy guardando u(0.5,0.5), sino en otra posición. Además, no lo estoy guardando para todo tiempo

Dejé corriendo varios incisos con una tolerancia de e-3


Estuvimos corriendo el código en una clase de consulta. María y yo no obtenemos los mismos resultados que Federico Teruel. Se ve que el problema de la lentitud del programa está en un problema con el código de la cátedra (o con Octave). Federico va a revisar el código y luego nos avisa.


Terminó el a

Para Re = 100, n1 = 80, dt = 0.5 NO CONVERGE. No estoy en el caso del profesor "Para Re=100 empecé a subir deltat y vi que para deltat=0.5 llego a un criterio de derivada de aprox. 10^-5 en 100 iteraciones o aprox. 4 minutos"

d. Con R1 = 1 el término advectivo tiene poca importancia. Para Re alto van a valores distinos a altos n1. 

Hay que tener cuidado la comparación entre métodos numéricos. Si comparo DC2 a distintos n1 con DC2 a n1 = 80, no debería sorprenderme que a n1 = 80 tengo un error cero y cerca tengo una convergencia bastante rápida.





























Cambié el rate de cambio en best_dt de 2 a 1.2. El mejor dt es
best_dt =0.1331. Con ese dt llega a 1.121e5 en 100 pasos de tiempo.

Cambié el criterio de va_a_converger. Ahora la pendiente se calcula en escala log para el eje vertical. Obtengo que el mejor dt es 5.1034. Con ese dt llega a 41.455 en 100 pasos de tiempo

El problema es que puede llegar a tomarle mucho tiempo encontrar el mejor dt si arranco de un dt demasiado chico.

Cambié el criterio de va_a_converger. Ahora calcula la pendiente ignorando el primer dato. Obtengo que el mejor dt es 0.1917. Con ese dt llega a 4.5e4. en 100 pasos de tiempo. Me quedo con el ajuste de todos los datos, a pesar de que pueda llegar a ser problemático.

Testeé un caso con n1 alto y le toma mucho tiempo hacer la búsqueda del mejor dt. Cambié el rate a 2 y llegué a 5.12 como mejor dt.


Corro el inciso b para n1 = 81 y Re = 5000. Mido cuánto tiempo le toma usando el mecanismo de búsqueda nuev
o y el anterior.

Con el nuevo (rate 2 y escala log en el eje vertical para el ajuste lineal) le toma . considerando también el tiempo necesario para encontrar tal dt

Con el viejo (rate 2 y escala normal para el ajuste lineal) le toma .




Cuando pruebo el caso n1 = 80, Re = 5000 no obtengo el mismo dt que el profesor, sino uno menor. Voy a analizar este caso en particular. Evoluciono 20 pasos de tiempo con distintos dt y grafico cada uno de ellos en Origin
dt = 0.005
dt = 0.02
dt = 0.04

En base a los gráficos podría cambiar la condición de va_a_converger
Cuando encuentre un punto con derivada menor va a calcular la difere

Lo dejé buscando el best_dt para el caso de n1 = 80 y Re 5000. La única diferencia es que puse un Ndeltat_max = 20




































En el inciso f se pide otra cosa. Lo voy a cambiar si me queda tiempo

En el e resolví un caso simplificado por las dudas


Mencionar en algún lado que el código fue escrito por la cátedra

Y si el error está en el promedio para calcular u(0.5,0.5)?

Asumí que QUICK es de orden 2 (podría ser de orden 3)

El e con Re 1000 llegó a un error de 0.7642











