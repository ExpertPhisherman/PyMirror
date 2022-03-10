# Imports
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# Initialize variables
x = sp.Symbol('x', real = True)
mir = x + 1
obj = x ** 2
xBounds = (-5, 5)
xStep = 0.1
decPlaces = 4

# Flatten list
def flatten(myList):
    try:
        return [i for subList in myList for i in subList]
    except:
        return myList

# Convert list items to float
def floatList(myList):
    return [float(i) for i in myList]

# Approximate function zeros using Householder's method
def householderZeros(func, var, guess, order, epochs):
    for epoch in range(epochs):
        guessOld = guess
        guessOld += order * ((1 / func).diff(var, order - 1).subs(var, guessOld) / (1 / func).diff(var, order).subs(var, guessOld))
        guessOld = guessOld.evalf()
        if guessOld == sp.nan:
            break
        else:
            guess = guessOld
    return guess

# Normal line to function at a point
def norm(func, var, pt):
    return -1 / func.diff(var).subs(var, pt) * (var - pt) + func.subs(var, pt)

# Reflect point about center
def reflect(pt, center):
    return 2 * center - pt

mir = sp.sympify(mir)
obj = sp.sympify(obj)

ptsMir_x = np.arange(xBounds[0], xBounds[1] + xStep, xStep)
print([norm(mir, x, ptMir_x) for ptMir_x in ptsMir_x])
ptsObj_x = [sp.solve(norm(mir, x, ptMir_x) - obj) for ptMir_x in ptsMir_x]
ptsRefl_x = [reflect(ptObj_x, ptMir_x) for ptObj_x, ptMir_x in zip(ptsObj_x, ptsMir_x)]

ptsMir_x = floatList(flatten(ptsMir_x))
ptsObj_x = floatList(flatten(ptsObj_x))
ptsRefl_x = floatList(flatten(ptsRefl_x))

ptsMir_y = [mir.subs(x, i) for i in ptsMir_x]
ptsObj_y = [obj.subs(x, i) for i in ptsObj_x]
ptsRefl_y = [reflect(ptObj_y, ptMir_y) for ptObj_y, ptMir_y in zip(ptsObj_y, ptsMir_y)]

ptsMir_y = floatList(flatten(ptsMir_y))
ptsObj_y = floatList(flatten(ptsObj_y))
ptsRefl_y = floatList(flatten(ptsRefl_y))

ptsRefl_x = ptsRefl_x[:len(ptsRefl_y)]

print('ptsObj_x:', ptsObj_x, '\n')
print('ptsMir_x:', ptsMir_x, '\n')
print('ptsRefl_x:', ptsRefl_x, '\n')

print(len(ptsObj_x), '\n')
print(len(ptsMir_x), '\n')
print(len(ptsRefl_x), '\n')

print('ptsObj_y:', ptsObj_y, '\n')
print('ptsMir_y:', ptsMir_y, '\n')
print('ptsRefl_y:', ptsRefl_y, '\n')

print(len(ptsObj_y), '\n')
print(len(ptsMir_y), '\n')
print(len(ptsRefl_y), '\n')

plt.plot(ptsMir_x, ptsMir_y, label = 'Mirror')
plt.plot(ptsObj_x, ptsObj_y, label = 'Object')
plt.plot(ptsRefl_x, ptsRefl_y, label = 'Reflection')
plt.legend()
plt.show()