



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
