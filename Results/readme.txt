
*****FOR n=8 DATA

columns are:
 
--graph number,	C_{max},	C(p=0 QAOA),	C(new modified QAOA),	P(C_max),	p,	beta_(qubit=0)/pi,	beta_(qubit=1)/pi,	..., beta_(qubit=n-1)/pi, gamma_edges/pi from  here on, 


The gamma for each edge are listed in left-to-right and top-to-bottom correspondance with the adjacency matrix upper-triangles in the graph files. For example, for a three-qubit graph the adjacency matrix could read 

01
1

This means there is an edge between qubits 0,2 (the 1 in the top right) and between qubits 1,2 (the one at the bottom).  For this example the gamma angles at the end of the file would be listed as gamma(edge = 0,2), gamma(edge = 1,2).


******For n=50,100 

columns are:

--graph number, C_max, min observed C during optimization, max observed C during optimization, beta angles, gamma angles

The beta and gamma angles are stored with the same ordering as the n=8 data

These have results for both regular QAOA and for ma-QAOA.

