import pyparsing as pp
import networkx as nx


def copyAttributes(A, B):
    for u in B:
        if "latent" in A.nodes[u]:
            B.nodes[u]["latent"] = True


def vertexFlowGraph(G):
    vfG = nx.DiGraph()
    for n in G.nodes():
        vfG.add_edge(n, n + "#", capacity=1)
        for desc in G.successors(n):
            vfG.add_edge(n + "#", desc)

    return vfG


def flowSet(G, S, T):
    if len(T) == 0 or len(S) == 0:
        return set(), set()
    G = vertexFlowGraph(G)

    # Create the incoming edges
    for s in S:
        G.add_edge("STARTFLOW", s)
    for t in T:
        G.add_edge(t + "#", "ENDFLOW")

    f, v = nx.algorithms.flow.maximum_flow(G, "STARTFLOW", "ENDFLOW")

    Sm = set()
    for n in v["STARTFLOW"]:
        if v["STARTFLOW"][n] > 0:
            Sm.add(n)

    Tm = set()
    for n in G.predecessors("ENDFLOW"):
        if v[n]["ENDFLOW"] > 0:
            Tm.add(n[:-1])

    return Sm, Tm


def matchBlock(G, S, T):
    tlen = 99999999
    while len(T) < tlen:
        # print("LOOP")
        tlen = len(T)
        Sp, Tp = flowSet(G, S, T)
        # Remove the ancestors of T that had no flow to them from S
        tNoFlow = T - Tp
        tAn = set()
        for v in tNoFlow:
            tAn = tAn.union(nx.ancestors(G, v))

        S = Sp - tAn
        T = Tp
        if len(S) == 0 or len(T) == 0:
            return set(), set()

    return S, T


def closestMinCut(G, S, T):
    if len(T) == 0 or len(S) == 0:
        return set()
    # Find the min-cut between S and T, closest to T
    G = G.copy()

    # Create the incoming edges
    for s in S:
        G.add_edge("STARTFLOW", s)
    for t in T:
        G.add_edge(t, "ENDFLOW")

    # G = G.reverse()
    # NOTE: We are exploiting the underlying flow algorithm of networkx to return the closest cutset.
    return nx.algorithms.connectivity.minimum_st_node_cut(G, "STARTFLOW", "ENDFLOW")


def flowGraph(G):
    fG = nx.DiGraph()
    for n in G:
        if "latent" in G.nodes[n]:
            # Create the bidirected edge between the nodes
            conn = list(G.successors(n))
            fG.add_edge(conn[0], conn[1] + "'")
            fG.add_edge(conn[1], conn[0] + "'")
        else:
            # Create the node, and links to its successors
            fG.add_edge(n, n + "'")
            for de in G.successors(n):
                fG.add_edge(n + "'", de + "'")
                fG.add_edge(de, n)
    return fG


def auxGraph(G, known=set()):
    aG = nx.DiGraph()
    for n in G:
        if "latent" in G.nodes[n]:
            if not n in aG:
                aG.add_node(n)
            aG.nodes[n]["latent"] = True
        else:
            aG.add_edge(n + "*", n)
            for pa in G.predecessors(n):
                if (pa, n) in known:
                    aG.add_edge(pa, n)
                else:
                    aG.add_edge(pa, n + "*")
    return aG


def auxFlowGraph(G, known=set()):
    aG = auxGraph(G, known)
    faG = flowGraph(aG)

    # The auxiliary flow graph has the "top" epsilon in its AVs only.
    for n in G:
        if not "latent" in G.nodes[n]:
            faG.remove_edge(n, n + "'")
    return faG


def ICvar(G, y, known=set()):
    faG = auxFlowGraph(G, known)

    yfaG = faG.copy()

    T = set()

    # Remove all parents of y* in yfaG, and find ancestors
    pred = set(yfaG.predecessors(y + "*'"))
    for pa in pred:
        if pa[-1] == "'":
            yfaG.remove_edge(pa, y + "*'")
            T.add(pa)
    anY = nx.ancestors(yfaG, y + "*'")

    S = set(
        [
            x + "*"
            for x in G
            if not x == y
            and not x in nx.descendants(G, y)
            and not "latent" in G.nodes[x]
        ]
    )

    S = S - anY

    C = closestMinCut(faG, S, T)

    S, _ = flowSet(faG, S, C)

    # Remove incoming edges to C
    cfaG = faG.copy()
    ied = set(cfaG.in_edges(C))
    #print(C, ied)
    cfaG.remove_edges_from(ied)

    Sm, Tm = matchBlock(cfaG, C, T)
    _, T = flowSet(cfaG, C - Sm, T - Tm)
    T = T | Tm
    # Returns 3 sets. S is all the source nodes. Next, is T\Tm, sink nodes that are not part of the match-block.
    # Finally, Tm are the sink nodes that *are* part of the match-block, and we can identify t->y for t in Tm
    return S, {t[:-1] for t in T}, {t[:-1] for t in Tm}


def ICID(G, known=set()):
    known = known.copy()
    lk = len(known) - 1
    while len(known) > lk:
        lk = len(known)
        for n in G:
            if not "latent" in G.nodes[n]:
                _, _, v = ICvar(G, n, known)
                for vi in v:
                    known.add((vi, n))
                # print(n, v)
    return known


# Manually inputting the graph is a task I heavily dislike,
# so I am including a simple parser for a simplified version of
# the text graph input used in fusion.
# That is, we simply draw the arrows:
#
# z->x
# x->y
# x--y
#
# The above 3 lines represent the instrumental variable.

# Set up the variable names - var represents a node name
varname = pp.Combine(
    pp.Word(pp.alphanums + "_", exact=1) +
    pp.Optional(pp.Word(pp.alphanums + "_"))
)
arrow = pp.Or(["--", "->"])
edge = pp.Group(varname + arrow + varname)
graphparser = pp.OneOrMore(edge)


def generateGraph(txt):
    parseresult = graphparser.parseString(txt)

    G = nx.DiGraph()
    for edge in parseresult:
        if edge[1] == "->":
            G.add_edge(edge[0], edge[2])
        else:
            # Uh oh, latent alert!
            latentName = "U({},{})".format(edge[0], edge[2])
            G.add_edges_from([(latentName, edge[0]), (latentName, edge[2])])
            G.nodes[latentName]["latent"] = True
    return G


def pg(G):
    gstring = []
    for node in G:
        if "latent" in G.nodes[node]:
            gstring.append("--".join(list(G.successors(node))))
        else:
            for succ in G.successors(node):
                gstring.append(node + "->" + succ)
    print(" ".join(gstring))


if __name__ == "__main__":
    G = generateGraph(
        # "z1->x1 z1->w z2->w w->x3 w->x2 w--y x1--y x2--y x3--y x1->y x2->y x3->y"
        # "w->z1 w->z2 z2->z1 z1->x1 z2->x2 x1->y x2->y x1--y x2--y w--y z1--w"
        # "x1->y x2->y x3->y x1--y x2--y x3--y w--y w->x2 w->x3 x1->w z1->x1 z2->w"
        "1->2 1->6 1--6 1--4 1->3 2->6 2->3 2->4 2->5 2--6 2--5 2--3 3->4 4->5"
    )

    # pg(vertexFlowGraph(G))
    # print(matchBlock(G, {"1", "2", "3"}, {"4", "5"}))

    # pg(auxGraph(G))

    # print(ICvar(G, "y"))
    print(ICvar(G, '5'))

    # print(closestMinCut(G, {"1", "2"}, {"5", "6"}))
