# This Python script holds the functionality to create a braid cover graph from a maximal chain
# of the Permutahedron

from collections import deque
from pyrsistent import v, pvector
import pygraphviz as pgv


# ------------- Draw Graph using PyGraphviz -----------------------
# Draw graph given by edges into name.jpg
def draw_graph(edges, name):
    graph = pgv.AGraph()
    graph.add_edges_from(edges)
    graph.graph_attr["overlap"] = "scale"
    graph.layout(prog="dot")
    graph.draw(name + '.jpg')


# Draw Hesse diagram of poset encoded by braid cover graph into name.jpg
# dim is the order of the Permutahedron in which the maximal chains can be found
def draw_poset(edges_BCG, name, dim):
    graph = pgv.AGraph()
    nodes = set()
    for edge in edges_BCG:
        for el in edge:
            chain = cox_to_nums(el, dim)
            for i in range(len(chain) - 1, 0, -1):
                nodes.add(chain[i])
                nodes.add(chain[i-1])
                graph.add_edge(chain[i], chain[i-1])

    graph.layout()
    graph.layout(prog="dot")
    graph.draw(name + '.jpg')
# -----------------------------------------------------------------


# ------------- Simple Persistent Union Find Data Structure -------
def uf_find(w_memo, x):
    memo = w_memo[0]
    if memo[x] != x:
        w_memo[0] = memo.set(x, uf_find(w_memo, memo[x]))
        return w_memo[0][x]
    else:
        return x


def uf_union(memo, x, y):
    x = uf_find(memo, x)
    y = uf_find(memo, y)
    if x == y:
        return memo
    elif x < y:
        return [memo[0].set(y, x)]
    else:
        return [memo[0].set(x, y)]


def uf_equal(memo, x, y):
    x = uf_find(memo, x)
    y = uf_find(memo, y)
    return x == y
# ---------------------------------------------------------------

# -------------- Helper functions -------------------------------
def is_braid_move(i, w):
    if w[i] == w[i+2]:
        if w[i+1] == w[i] + 1 or w[i+1] == w[i] - 1:
            return True
    return False


def do_braid_move(i, w):
    w1 = list(w)
    w1[i] = w1[i+1]
    w1[i+1] = w1[i+2]
    w1[i+2] = w1[i]
    return tuple(w1)


def is_comm_move(i, w):
    if abs(w[i] - w[i+1]) > 1:
        return True
    else:
        return False


def do_comm_move(i, w):
    w1 = list(w)
    temp = w1[i]
    w1[i] = w1[i+1]
    w1[i+1] = temp
    return tuple(w1)


# Convert a sequence of inversions to a Coxeter word
def seq_to_coxeter(s, dim):
    cox = list()
    n = list(range(1,dim+1))
    for i in range(len(s)):
        j = n.index(s[i][0])
        cox.append(j+1)
        tmp = n[j]
        n[j] = n[j+1]
        n[j+1] = tmp
    return cox


# Convert a Coxeter word to a sequence of inversions
def coxeter_to_seq(c, dim):
    s = list()
    n = list(range(1,dim+1))
    for i in c:
        i -= 1
        j1 = n[i]
        j2 = n[i+1]
        s.append((j1,j2))
        tmp = n[i]
        n[i] = n[i+1]
        n[i+1] = tmp
    return tuple(s)


# Convert a Coxeter word to the denoted permutation
def cox_to_nums(c, dim):
    n = list(range(1, dim+1))
    chain = [tuple(n)]
    for i in c:
        tmp = n[i-1]
        n[i-1] = n[i]
        n[i] = tmp
        chain.append(tuple(n))
    return chain


def subchain(w, dim):
    r = list()
    for inv in w:
        if inv[1] == dim:
            continue
        else:
            r.append(inv)
    return tuple(r)
# ------------------------------------------------------


# Return the braid covers of w been disjoint from p to len(w)-p
def get_braid_cover(w, p):
    if p > len(w) // 2:
        return set()

    path = [pvector(range(0, len(w)))]
    visited = set()
    visited.add(w)
    d = deque()
    d.append((w, path))
    res = set()

    while len(d) > 0:
        u, u_path = d.popleft()
        for i in range(p, len(w)-p-2):
            if is_braid_move(i, u):
                if not uf_equal(u_path, i, i+2) and not uf_equal(u_path, i, i+1) and not uf_equal(u_path, i+1, i+2):
                    u1 = do_braid_move(i, u)
                    if u1 not in visited:
                        u1_path = uf_union(u_path, i, i+1)
                        u1_path = uf_union(u1_path,i+1, i+2)
                        if uf_equal(u1_path, p, len(w)-p-1):
                            res.add(u1)
                        else:
                            d.append((u1, u1_path))
                        visited.add(u1)
        for i in range(p, len(w)-p-1):
            if is_comm_move(i, u):
                if not uf_equal(u_path, i, i+1):
                    u2 = do_comm_move(i, u)
                    if u2 not in visited:
                        u2_path = uf_union(u_path, i, i+1)
                        if uf_equal(u2_path, p, len(w)-p-1):
                            res.add(u2)
                        else:
                            d.append((u2, u2_path))
                        visited.add(u2)
    return res


# Return the braid cover graph obtained from w
def braid_cover_graph(w, restricted=False, K = set()):
    d = deque()
    d.append(w)
    visited = set()
    n = len(w)
    h = n // 2 + 1
    graph = []

    while len(d) > 0:
        w1 = d.popleft()
        for i in range(h):
            succs = get_braid_cover(w1, i)
            for s in succs:
                if not restricted or subchain(coxeter_to_seq(s, 5), 5) in K:
                    if s not in visited:
                        graph.append((w1, s))
                        visited.add(s)
                        d.append(s)
                    else:
                        graph.append((w1, s))
    return graph


# Example creating a chain cover poset for P4
def example_P4():
    w = (1,2,1,3,2,1)
    bcg = braid_cover_graph(w)
    draw_graph(bcg, "BCG(121321)")
    draw_poset(bcg, "P(121321)", 4)


# Example creating a restricted braid cover graph for P5.
# As explained in the thesis, this graph does not readily
# give rise to a chain cover poset of P5.
def example_P5():
    w1 = (1,2,1,3,2,1,4,3,2,1)
    w2 = (3,2,1,2,3,4,3,2,1,3)
    # Subset of maximal chains of chain cover poset of P4
    K1 = {
        ((1, 2), (1, 3), (2, 3), (1, 4), (2, 4), (3, 4)),
        ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)),
        ((2, 3), (2, 4), (3, 4), (1, 4), (1, 3), (1, 2)),
        ((2, 3), (1, 3), (2, 4), (1, 4), (3, 4), (1, 2)),
    }
    # Remaining two chains from chain cover poset of P4
    K2 = {
        ((3, 4), (1, 2), (1, 4), (1, 3), (2, 4), (2, 3)),
        ((3, 4), (2, 4), (1, 4), (1, 2), (1, 3), (2, 3))
    }

    restricted_bcg1 = braid_cover_graph(w1, True, K1)
    restricted_bcg2 = braid_cover_graph(w2, True, K2)

    bcg_P5 = restricted_bcg1 + restricted_bcg2
    draw_graph(bcg_P5, "Restricted_BCG_P5")
    # It can be seen that the resulting poset is not a chain cover poset yet;
    # There are to many chains present
    draw_poset(bcg_P5, "P_P5", 5)


if __name__ == '__main__':
    example_P5()
