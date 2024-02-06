from itertools import *

# Create a SAT instance stating that the graph given by
# the vertices 'V' and edges 'E' has a symmetric chain decomposition
# using the solver 'solver' from the PySAT framework.
# The rank of each vertex is given by the dictionary 'level'.
# We call the returned tuple from this function a 'solver_obj'
# to be passed to the function 'get_solution' below.
def scd_instance(V, E, level, solver, output=False):
    level_to_vtx = {}
    for v in V:
        l = level[v]
        if l in level_to_vtx:
            level_to_vtx[l].append(v)
        else:
            level_to_vtx[l] = [v]

    height = len(level_to_vtx)
    width = len(level_to_vtx[height // 2])
    if output:
        print("Height: " + str(height))
        print("Width: " + str(width))

    # if number of levels is odd, we need loops on middle level
    if height % 2 == 1:
        for v in level_to_vtx[height // 2]:
            E.add((v, v))

    v_id = {}
    counter = 0
    for v in V:
        v_id[v] = counter
        counter += 1

    N = len(V)
    X = lambda i,v: 1 + v_id[v] + N*i

    chain_start = chain_start_positions(level_to_vtx, height)

    # add all assumptions to solver
    if output: print("Assuming possible vertices for each chain")
    assume_vtx_in_chain(chain_start, level_to_vtx, width, height, X, solver)
    if output: print("Assuming that all vertices are visited")
    assume_vtx_used(chain_start, level_to_vtx, width, height, X, solver)
    if output: print("Assuming that vertices are only used once")
    assume_vtx_once(chain_start, level_to_vtx, width, height, X, solver)
    if output: print("Assuming that only valid edges are used")
    assume_valid_edges(chain_start, level_to_vtx, E, width, height, X, solver)

    return X, chain_start, level_to_vtx, width, height


# Compute 'n' possible solutions to the SAT instance given by 'solver_obj'
# using the solver 'solver'
def get_solution(n, solver, solver_obj, output=False):
    X = solver_obj[0]
    chain_start = solver_obj[1]
    level_to_vtx = solver_obj[2]
    width = solver_obj[3]
    height = solver_obj[4]

    # ----- Set up solver -----
    solution = list()
    solution_iter = solver.enum_models()

    counter = 0
    for sol in solution_iter:
        if output: print("Solution " + str(counter + 1))
        solution.append([])
        sol = set(sol)
        for i in range(width):
            l0 = chain_start[i]
            chain = set()
            for l in range(l0, height-l0):
                for v in level_to_vtx[l]:
                    if X(i, v) in sol:
                        chain.add(v)
            chain = sorted(chain)
            solution[counter].append(chain)
        if output: print("Solution end")
        counter += 1
        if n != -1 and counter >= n:
            break
    return solution

# ------------- Helper functions -------------------------------------------------------

# Compute for each chain on which rank the chain starts. The total number of chains
# is given by the width of the poset.
def chain_start_positions(level_to_nodes, height):
    start_positions = list(repeat(0, len(level_to_nodes[0])))
    for l in range(1, (height // 2) + 1):
        for n in range(len(level_to_nodes[l]) - len(level_to_nodes[l - 1])):
            start_positions.append(l)
    return start_positions


def assume_vtx_in_chain(chain_start, level_to_vtx, width, height, X, solver):
    for i in range(width):
        l0 = chain_start[i]
        # reachable vertices
        for l in range(l0, height-l0):
            solver.add_clause([X(i,v) for v in level_to_vtx[l]])
        # not reachable vertices
        for l in chain(range(l0), range(height-l0, height)):
            for v in level_to_vtx[l]:
                solver.add_clause([-X(i,v)])


def assume_vtx_used(chain_start, level_to_vtx, width, height, X, solver):
    for l in range(height):
        for v in level_to_vtx[l]:
            solver.add_clause([X(i,v) for i in range(width) if chain_start[i] <= l < height - chain_start[i]])


def assume_vtx_once(chain_start, level_to_vtx, width, height, X, solver):
    for i1, i2 in combinations(range(width),2):
        l0 = chain_start[i2]
        assert(chain_start[i1] <=l0)
        for l in range(l0, height-l0):
            for v in level_to_vtx[l]:
                solver.add_clause([-X(i1,v), -X(i2, v)])


def assume_valid_edges(chain_start, level_to_vtx, E, width, height, X, solver):
    for l0 in range(height-1):
        l1 = l0 + 1
        for v0 in level_to_vtx[l0]:
            for v1 in level_to_vtx[l1]:
                for i in range(width):
                    if chain_start[i] <= l0 and l1 <= height - chain_start[i]:
                        if not (v0, v1) in E:
                            solver.add_clause([-X(i,v0), -X(i, v1)])
# -----------------------------------------------------------------------
