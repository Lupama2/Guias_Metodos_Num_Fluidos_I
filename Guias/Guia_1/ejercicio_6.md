Readme.md

DUDAS:

*El método de Padé me da que no necesito en ningún caso usar b2 y me queda una diferencia centrada de orden 2 en ambos extremos. Es así? Sí.


*La solución tiene que ser para K=6? o para K grande? Cómo obtener la solución?

Puedo suponer que la solución exacta es aquella con N suficientemente grande? Tal que el "error" definido como la diferencia entre considerar N y N+1 (o 2N) es menor a un error relativo determinado?

NO. SE SUPONE QUE SABEMOS RESOLVERLA EXACTAMENTE. Tenemos que encontrar una solución homogénea y una particular. 



*e. Cómo determino el orden de truncamiento en ese punto?

Hago un gráfico de log error en función de log h?? Sí
Error_j = |y_j - y^~_j|
y grafico Error_j en función de DeltaX en escala log-log. Padé es orden 4 y DFC es orden 2.

NO ENTIENDO EL COMPORTAMIENTO DEL ERROR DE TRUNCAMIENTO EN FUNCIÓN DEL h

Arreglar en el código:
*X Creo que hay un error en el código de Padé porque la solución cambia bruscamente con Nn. Encontré el problema pero las soluciones convergen a funciones distintas. Para saber dónde está el error tengo que graficar la función solución.

*X Ver si ambas aproximaciones convergen a la solución exacta para N grande. No lo hacen, creo que hay un error en alguna de ellas.

Encontré el error en ambas. Ahora van a converger a la misma función. Falta ver que esa sea la función solución
