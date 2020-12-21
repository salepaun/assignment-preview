#! /usr/bin/python3


from random import random
from sys import argv
from math import sqrt

import numpy as np


def random_selection(*args):
    rnd = []
    for a in args:
        mn, mx = a
        rnd.append(random() * (mx -  mn) + mn)
    return tuple(rnd)



def generate_particles(num, base, x, y, size, psize, mass, cmin, cmax, vrange):
    """
    """
    gap = 1 * psize
    step = 2*psize+gap
    a = int(sqrt(num));
    xmin = x - step * a/2
    ymin = y - step * a/2
    xmax = xmin + a * step
    ymax = ymin + a * step
    X = np.asarray((x,y))
    Xmin = np.asarray((xmin, ymin))
    Xmax = np.asarray((xmax, ymax))
    dmax = np.linalg.norm(Xmax - Xmin) / 2
    vxy_range = tuple((vrange, vrange))

    xc = xmin
    i = 0
    while xc < xmax:
        yc = ymin
        while yc < ymax:
            d = np.linalg.norm(X - np.asarray((xc,yc)))
            colors= tuple([cmin[i] + d*(cmax[i]-cmin[i])/dmax for i in range(len(cmin))])
            r,g,b = tuple([i if i<=1.0 else 1 for i in colors])
            vx,vy = random_selection(*vxy_range)

            print('  <particle m="%5.3f" px="%5.3f" py="%5.3f" vx="%5.3f" vy="%5.3f" fixed="0" radius="%5.3f"/>' \
                    % (mass, xc, yc, vx, vy, psize))
            print('  <particlecolor i="%i" r="%5.3f" g="%5.3f" b="%5.3f"/>' \
                    % (i+base, r, g, b))

            i += 1
            yc = yc + step

        xc = xc + step

    print('</scene>')



def help():
    """
    """

    print("""
This program is a helper for creative assignment for t2m2.
It generates objects at certain location, with "num" particles of "size".

Usage: {0} num base x,y size psize mass vrange
""".format(argv[0]))



def main(name, args):
    """
    """
    num = 100
    base = 24
    x,y = ((0,0))
    size = 2
    psize = 0.1
    mass = 1
    cmin = ((0.1,0.1,0.1))
    cmax = ((0.1,0.1,0.9))
    vrange = ((-10,10))


    if args:
        if args[0] == '-h':
            help()
            return

        for i, arg in enumerate(args):
            if i==0:
                num = int(arg)
            elif i==1:
                base = int(arg)
            elif i==2:
                x,y = tuple([float(i) for i in arg.split(',')])
            elif i==3:
                size = float(arg)
            elif i==4:
                psize = float(arg)
            elif i==5:
                mass = float(arg)
            elif i==6:
                cmin = tuple([float(i) for i in arg.split(',')])
            elif i==7:
                cmax = tuple([float(i) for i in arg.split(',')])
            elif i==8:
                vrange = tuple([float(i) for i in arg.split(',')])
            else:
                break

    generate_particles(num, base, x, y, size, psize, mass, cmin, cmax, vrange)



if __name__ == "__main__":
    main(argv[0], argv[1:])

