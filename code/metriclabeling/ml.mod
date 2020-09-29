param N integer > 0;  # number of nodes
param M integer > 0;  # number of edges
param LabelNum integer > 0;  # number of labels

set V := 0 .. N-1;  #the set of vertices 
set L := 0 .. LabelNum-1;  #the set of labels

param EDGES{u in V,v in V}; #The edges in the graph

param COST{l in L,u in V}; # Costs

param DISTANCE{l1 in L,l2 in L}; #Distances

var X1{u in V , l in L} >= 0 ; #if X1[u,l]=1 then u is assigned label l.

var X2{u in V , l1 in L, v in V, l2 in L} >=0; #u is assigned label l1 and v is assigned label l2.




# OBJECTIVE FUNCTION
minimize totalcost: (sum{u in V,l in L }(COST[l,u]*X1[u,l]) + sum{u in V,l1 in L,v in V,l2 in L}(EDGES[u,v]*DISTANCE[l1,l2]*X2[u,l1,v,l2]));


# CONSTRAINTS

# Each vertex is assigned only one label.
 subject to uniqueness{u in V}: sum{l in L}(X1[u,l]) = 1;

# Assigned labels do not violate each other
  subject to nonviolation{u in V,v in V,l1 in L}: sum{l2 in L}(X2[u,l1,v,l2]) = X1[u,l1];

#Symmetry for label assignments
  subject to symmetry {u in V,v in V,l1 in L,l2 in L}: X2[u,l1,v,l2] = X2[v,l2,u,l1]; 


solve;

# OUTPUT


printf "#Mincost = %.4f \n", (sum{u in V,l in L }(COST[l,u]*X1[u,l]) + sum{u in V,l1 in L,v in V,l2 in L}(EDGES[u,v]*DISTANCE[l1,l2]*X2[u,l1,v,l2]));   
    printf {i in V, j in L  }: "%d %d %f \n", i, j , X1[i, j] ; 

end;