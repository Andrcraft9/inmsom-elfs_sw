from hilbert_curve import d2xy
import numpy as np
import math as mth
import matplotlib.pyplot as plt
import sys

def get_random_color(pastel_factor = 0.5):
    return [(x+pastel_factor)/(1.0+pastel_factor) for x in [np.random.uniform(0,1.0) for i in [1,2,3]]]

def color_distance(c1,c2):
    return sum([abs(x[0]-x[1]) for x in zip(c1,c2)])

def generate_new_color(existing_colors,pastel_factor = 0.1):
    max_distance = None
    best_color = None
    for i in range(0, 100):
        color = get_random_color(pastel_factor = pastel_factor)
        if not existing_colors:
            return color
        best_distance = min([color_distance(color,c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    
    return best_color

def read_decomposition(namef):
    data = {}
    with open(namef) as f:
        headIsRead = False
        for line in f: # read lines
            if headIsRead:
                x = int(line.split()[0])
                y = int(line.split()[1])
                val = int(line.split()[2])
                #weight = int(float(line.split()[3]))
                data[x, y] = val
            else:
                bx = int(line.split()[0])
                by = int(line.split()[1])
                px = int(line.split()[2])
                py = int(line.split()[3])
                headIsRead = True
    
    return bx, by, px, py, data

if ( __name__ == '__main__' ):
    print("Decomposition program!")
    
    namef = sys.argv[1]
    if len(sys.argv) > 2:
        msize = int(sys.argv[2])
    else:
        msize = 14.0
    print("Msize: ", msize)

    bx, by, px, py, data = read_decomposition(namef)
    blocks = bx*by
    print("File %s was read! bx = %i, by = %i Total blocks: %i" % (namef, bx, by, blocks))
    print("Proc grid: px = %i, py = %i" % (px, py))

    xdata = {i : [] for i in range(px*py)}
    ydata = {i : [] for i in range(px*py)}

    exist_colors = []
    pltcolors = {}
    for val in range(px*py):
        color = generate_new_color(exist_colors)
        exist_colors.append(color)
        pltcolors.update({val: exist_colors[len(exist_colors) - 1]})
        #pltcolors.update({val: np.random.rand(3,)})

    m = int(mth.log(blocks, 2))
    print("Hilbert Curve index: %i" % m)
    n = 2**m
    print("Total points: %i" % n)

    xplt = []
    yplt = []
    xmax, ymax = 0, 0
    for d in range(n):
        x, y = d2xy(m, d)
        
        if x > xmax:
            xmax = x
        if y > ymax:
            ymax = y

        xplt.append(x)
        yplt.append(y)

        val = data[x+1, y+1]
        if val != -1:
            xdata[val].append(x)
            ydata[val].append(y)
    
    print("xmax = %i and ymax = %i" % (xmax, ymax))
    plt.plot(xplt, yplt, 'ks', markersize=msize)
    plt.plot(xplt, yplt, 'k')

    for val in range(px*py):
        plt.plot(xdata[val], ydata[val], marker="s", markersize=msize, linestyle='None', color=pltcolors[val])
        
    plt.title("Hilbert curve, blocks = %i, cores = %i" % (blocks, px*py))
    plt.xticks(np.linspace(0.5, xmax+0.5, 16))
    plt.yticks(np.linspace(0.5, ymax+0.5, 16))
    plt.xlim((-0.5, xmax+0.5))
    plt.ylim((-0.5, ymax+0.5))
    plt.grid(True, linewidth=2)
    plt.show()