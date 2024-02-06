# This program demonstrates that the Coxeter group $H_3$
# does not admit a symmetric chain decomposition.

from collections import deque
from SCD import scd_instance, get_solution
from sys import *
from sage.all import *
from pysat.solvers import Cadical153, Glucose4


if __name__ == '__main__':
    solver = Glucose4()

    G = CoxeterGroup(['H', 3])
    min_el = G.from_reduced_word([])

    # vertices and edges
    V = set()
    E = set()
    level = {}

    d = deque()
    d.append(min_el)

    print("Started creating graph")
    while len(d) > 0:
        v0 = d.popleft()
        v0_red = tuple(v0.reduced_word())
        V.add(v0_red)
        level[v0_red] = v0.length()
        for v1 in v0.upper_covers():
            v1_red = tuple(v1.reduced_word())
            if not v1_red in V:
                d.append(v1)
            E.add((v0_red, v1_red))
    print("Finished creating graph")

    sol_obj = scd_instance(V, E, level, solver)
    sol = get_solution(1, solver, sol_obj)

    # The output will be empty because no symmetric chain decomposition of $H_3$ will be found
    for s in sol:
        for c in s:
            print(c)