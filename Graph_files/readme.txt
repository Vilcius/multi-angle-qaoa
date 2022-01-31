The graph files contain lists of the adjacency matrix upper-triangles for each graph.  These contain all edges in the graphs. 

For example, the upper triangle of an adjacency matrix for a 3-node graph could read

01
1

The first row shows edges connected to qubit=0.  The first "0" entry indicates there is not an edge (0,1).  The second "1" entry indicates there is an edge (0,2).  

The next row shows the connections to qubit=1, apart from edges in previous rows.  The "1"  in this row shows there is an edge (1,2). 

So for this example, the total graph has two edges, (0,2) and (1,2).  

The same type of description holds for each graph in these files.  For n=50, n=100, there are also separate files listing the optimal cut value C_max. 
