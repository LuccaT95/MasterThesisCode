# This program counts the total number of symmetric chain decompositions
# which are permitted by $P_4$

from collections import deque
from SCD import scd_instance, get_solution
from sys import *
from sage.all import *
from pysat.solvers import Cadical153, Glucose4

if __name__ == '__main__':
    solver = Cadical153()

    G = CoxeterGroup(['A', 3], implementation="permutation")
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
    # We pass the option -1 to obtain all possible models (satisfiable assignments) for the created SAT problem
    sol = get_solution(-1, solver, sol_obj)

    # The output will be 192 meaning that there are 192 satisfiable assignments to our posed SAT problem.
    # By our SAT encoding the chains of a symmetric chain decomposition are ordered. As a result, swapping chains of
    # equal length in the symmetric chain decomposition will be counted as different solutions to the SAT problem.
    # However, we are not interested in this distinction. Since there are 2 chains of length 3 and 2 chains of length 5
    # in any symmetric chain decomposition of $P_4$, 192 / 4 = 48 gives the actual number of different symmetric chain
    # decompositions.
    print(len(sol))
