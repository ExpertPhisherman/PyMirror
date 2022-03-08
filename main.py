# Imports
import sympy as sp

# Initialize variables
x = sp.Symbol('x')
mir = 2 * x + 1
obj = x - 1
xBounds = (-5, 5)
xStep = 0.1
decPlaces = 4

# Custom range method to support floats
def fRange(start, stop, increment):
    out = []
    while start <= stop:
        out += [start]
        start += increment
    return out

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
        #print(guess)
    return guess

# Normal line to function at a point
def norm(func, var, pt):
    return -1 / func.diff(var).subs(var, pt) * (var - pt) + func.subs(var, pt)

# Reflect point about center
def reflect(pt, center):
    return 2 * center - pt

ptsMir = fRange(xBounds[0], xBounds[1], xStep)
ptsObj = [householderZeros(norm(mir, x, ptMir) - obj, x, 4, 1, 8) for ptMir in ptsMir]
ptsRefl = [reflect(ptObj, ptMir) for ptObj, ptMir in zip(ptsObj, ptsMir)]

print('ptsMir:', ptsMir)
print('ptsObj:', ptsObj)
print('ptsRefl:', ptsRefl)