import numpy as np
import sys

def get_mesh_data(namef):
    data = {}
    with open(namef) as f:
        headIsRead = False
        for line in f: # read lines
            if headIsRead:
                x = int(line.split()[0])
                y = int(line.split()[1])
                #val = int(line.split()[2])
                weight = int(float(line.split()[3]))
                data[x, y] = weight
            else:
                bx = int(line.split()[0])
                by = int(line.split()[1])
                headIsRead = True
    
    return bx, by, data

def mesh2graph(data, m, n):
    edges_count = 0
    nodes_count = m*n

    xrange = range(1, m+1)
    yrange = range(1, n+1)

    total_edges = []
    for i in xrange:
        for j in yrange:
            
            edges = []
            # 4 point scheme:
            x = i
            y = j-1
            if x in xrange and y in yrange:   
                edges.append(y + (x-1)*n)
                
            x = i
            y = j+1
            if x in xrange and y in yrange:
                edges.append(y + (x-1)*n)
            x = i-1
            y = j
            if x in xrange and y in yrange:
                edges.append(y + (x-1)*n)
            x = i+1
            y = j
            if x in xrange and y in yrange:
                edges.append(y + (x-1)*n)

            edges_count = edges_count + len(edges)

            # vertex weight 
            w = data[i, j]
            total_edges.append([w, *edges])

    edges_count = int(edges_count / 2)
    #print("Nodes: ", nodes_count, " Edges: ", edges_count)
    print(nodes_count, edges_count, 10)

    for e in total_edges:
        print(*e)

if ( __name__ == '__main__' ):

    # Read mesh in model format
    mesh_name = sys.argv[1]
    bx, by, data = get_mesh_data(mesh_name)

    if len(sys.argv) > 2:
        # Read decomposition from METIS and print decomposition in model format
        part_name = sys.argv[2]

        part = []
        with open(part_name) as f:
            for line in f:
                part.append(int(line.split()[0]))

        xrange = range(1, bx+1)
        yrange = range(1, by+1)
        for i in xrange:
            for j in yrange:
                node = j + (i - 1)*by
                if data[i, j] == 0:
                    print(i, j, -1, data[i, j])
                else:
                    print(i, j, part[node-1], data[i, j])
    else:
        # Print mesh as graph (METIS format)
        mesh2graph(data, bx, by)
