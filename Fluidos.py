#Importamos paquetes:
import numpy as np                 #Para usar la función exponencial y crear una matriz.
import matplotlib.pyplot as plt    #Para graficar la solución.
plt.style.use('bmh')               #Cambiamos el estilo del gráfico por gusto personal.

#Definimos las constantes:
p0 = 101.325 #N/m^2
a = 1.5 #m
b = 3.5 #m
g = 9.8 #m/s^2
h = 10 #m
T0 = 293.15 #K
R = 8.314 #J/(K mol)

def fuerza(H):                     #Definimos una función que calcula F_p.
    return (b * p0 / 2) * (-a + H + (R * T0 / g) * (np.exp(g * (h + a - H) / (R * T0)) - np.exp(g * h / (R * T0))))

def biseccion(F_p, tol=1e-6, n=100):#Definimos un método de bisección que resuelve numéricamente la ecuación trascendental para H.
    #Extraemos los valores máximo y mínimo de H para definir un intervalo en el cual buscar la solución.
    H_min = 0
    H_max = a
    for _ in range(n): #Iteramos el proceso hasta 100 veces con índice "_" ya que no necesitamos registrar el número de la iteración.
        H_med = (H_min + H_max)/2 #Calculamos el valor intermedio de los valores posibles de H.
        F_p_med = fuerza(H_med) #Calculamos F_p para ese valor de H.
        
        if abs(F_p_med - F_p) < tol: #Comparamos el valor entregado de F_p con el calculado en la línea anterior.
            return H_med #Si la diferencia es menor que la tolerancia, H_med es la solución de H para el valor de F_p dado.
        #Calculamos si el valor está en la mitad superior o inferior del intervalo.
        if (fuerza(H_min)-F_p)*(F_p_med-F_p)<0:#Revisamos esto planteando un posible cambio de signo.
            H_max = H_med #Si hay cambio, acotamos el invervalo a su mitad superior.
        else:
            H_min = H_med #Si no hay cambio, acotamos el intervalo a su mitad inferior.
    #De esta forma se repite el proceso acotando el intervalo hasta encontrar el valor de H para el P_b dado.
    #Si después de las 100 iteraciones no se encuentra H se dice que diverge. 
    return np.nan #Se regresa un nan para evitar errores al trabajar con arreglos, se agregó por completitud, ya que el código no debería divergir.

H_l = np.linspace(0, a, 500)       #Limitamos los valores de H por el sentido físico del sistema.
F_p_l = fuerza(H_l)                #Limitando nuestra variable a sus valores que entreguen resultados con sentido.

H_nun = [biseccion(F_p) for F_p in F_p_l]
#Evaluamos término a término para poder hacer las comparaciones de los if y después agrupamos los resultados en un arreglo.

#Graficamos H en función de F_p para todos los valores con sentido físico.
plt.figure(figsize=(8,6))
plt.plot(F_p_l, H_nun)
plt.xlabel("$F_p$ [N]")
plt.ylabel("$H$ [m]")
plt.savefig("H en función de F_p.pdf")
plt.show()