"""
TODO: module/program docstring
"""

import turtle
import random


def fractal(x0, y0, x1, y1, scale):
    width = x1 - x0
    if width < 1:
        turtle.goto(x1, y1)
    else:
        x_mid = (x0 + x1) / 2
        y_mid = (y0 + y1) / 2 + random.gauss(0, width*scale)
        fractal(x0, y0, x_mid, y_mid, scale)
        fractal(x_mid, y_mid, x1, y1, scale)

def main():
    """The entry point of the program."""
    turtle.tracer(0, 0)
    
    for x in range(10, 0, -1):
        turtle.goto(-400, 0)
        c = x / 10
        turtle.color((c,c,c))
        turtle.begin_fill()
        fractal(-400, 0, 400, 0, scale=0.02*x)
        turtle.end_fill()
    
    turtle.mainloop()

if __name__ == "__main__":
    main()
