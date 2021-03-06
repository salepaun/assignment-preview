#! /usr/bin/python3


from random import random
from sys import argv


def random_selection(*args):
    rnd = []
    for a in args:
        mn, mx = a
        rnd.append(random() * (mx -  mn) + mn)
    return tuple(rnd)



def generate_balls(num, base):
    """
    """
    mass_range = tuple(((0.2, 10),))
    size_range = tuple(((0.03, 0.30),))
    color_range = tuple(((0.1, 1.0), (0.1, 1.0), (0.1, 1.0)))
    x_range = tuple(((-10, 10),))
    y_range = tuple(((-10, 10),))
    vxy_range = tuple(((-1,1),(-1,1)))


    for i in range(num):
        m, = random_selection(*mass_range)
        ri, = random_selection(*size_range)
        r, g, b = random_selection(*color_range)
        x, = random_selection(*x_range)
        y, = random_selection(*y_range)
        vx,vy = random_selection(*vxy_range)
        print('  <particle m="%5.2f" px="%5.2f" py="%5.2f" vx="%5.2f" vy="%5.2f" fixed="0" radius="%5.2f"/>' \
                % (m, x, y, vx, vy, ri))
        print('  <particlecolor i="%i" r="%5.3f" g="%5.2f" b="%5.2f"/>' \
                % (i+base, r, g, b))
    print('</scene>')



def main(name, args):
    """
    """
    num = 100
    base = 4

    if args:
        for i, arg in enumerate(args):
            if i==0:
                num = int(arg)
            elif i==1:
                base = int(arg)
            else:
                break

    generate_balls(num, base)



if __name__ == "__main__":
    main(argv[0], argv[1:])

