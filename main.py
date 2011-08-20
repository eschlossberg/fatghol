#! /usr/bin/env python
"""Command-line front-end for `rg.py`.
"""
__docformat__ = 'reStructuredText'


from rg import all_graphs,vertex_valences_for_given_g_and_n,Graph
import sys
from utils import positive_int


def graph_to_xypic(graph):
    r"""Print XY-Pic code snippet to render graph `graph`.

    Examples::
      #>>> print graph_to_xypic([[3,2,1],[3,1,2]])
      \xy 0;<2cm,0cm>:%
      (1,1)="v1",%
      (-1,1)="v2",%
      "v1",{\xypolygon4"v1l"{~:{(1.20,0):(0,-1)::}~={90}~>{}}},%
      "v2",{\xypolygon4"v2l"{~:{(1.20,0):}~={270}~>{}}},%
      "v1"*\txt{[321]},%
      "v2"*\txt{[312]},%
      "v1";"v2"**\crv{"v1l3"&"v2l3"},%
      "v1";"v2"**\crv{"v1l2"&"v2l1"},%
      "v1";"v2"**\crv{"v1l1"&"v2l2"},%
      \endxy
    #>>> print graph_to_xypic([[1,2,3,1],[1,2,3,2],[1,2,3,3]])
      \xy 0;<2cm,0cm>:%
      {\xypolygon3"v"{~={0}~>{}}},% mark v1,v2,v3
      "v1",{\xypolygon8"v1l"{~:{(1.20,0):}~={90}~>{}}},%
      "v2",{\xypolygon8"v2l"{~:{(1.20,0):}~={210}~>{}}},%
      "v3",{\xypolygon8"v3l"{~:{(1.20,0):}~={330}~>{}}},%
      "v1"*\txt{[1231]},
      "v2"*\txt{[1232]},
      "v3"*\txt{[1233]},
      "v1";"v2"**\crv{"v1l2"&"v2l1"},%
      "v1";"v3"**\crv{"v1l3"&"v3l1"},%
      "v2";"v3"**\crv{"v2l3"&"v3l2"},%
      "v1";"v1"**\crv{"v1l1"&"v1l4"},%
      "v2";"v2"**\crv{"v2l2"&"v2l4"},%
      "v3";"v3"**\crv{"v3l3"&"v3l4"},%
      \endxy
    """
    def vertex_label(v):
        return '['+("".join(map(str,v)))+']'
    label = map(vertex_label, graph)
    K = len(graph) # number of vertices
    result = r'\xy 0;<2cm,0cm>:%'+'\n'
    # put graph vertices on a regular polygon
    if K < 3:
        result += '(1,1)="v1",%\n(-1,1)="v2",%\n'
    else:
        result += r'{\xypolygon%d"v"{~={0}~>{}}},%% mark vertices' \
                  % (max(K,3)) \
                  + '\n'
    # mark invisible "control points" for bezier curves connecting vertices
    def rotation_angle(K,k):
        return 90+(k-1)*360/K
    for k in range(1,K+1):
        # "vK",{\xypolygonL"vKl"{~{1.20,0):~={...}~>{}}},%
        result += r'"v%d",{\xypolygon%d"v%dl"{~:{(1.20,0):q}~={%d}~>{}}},%%' \
                  % (k, 2*len(graph[k-1])-2, k, rotation_angle(K,k)) \
                  + '\n'
    for k in range(1,K+1):
        # "vK"*\txt{[...]},%
        result += r'"v%d"*\txt{%s},%%' \
                  % (k, label[k-1]) \
                  + '\n'
    for l in range(1,num_edges(graph)+1):
        (v1,v2) = endpoints(l,graph)
        if v1 != v2:
            result += r'"v%d";"v%d"**\crv{"v%dl%d"&"v%dl%d"},%%' \
                      % (v1+1, v2+1, v1+1, 1+graph[v1].index(l),
                         v2+1, 1+graph[v2].index(l)) \
                      + '\n'
        else:
            h = graph[v1].index(l)
            result += r'"v%d";"v%d"**\crv{"v%dl%d"&"v%dl%d"},%%' \
                      % (v1+1, v2+1, v1+1, h+1, v2+1, 1+graph[v1].index(l,h+1)) \
                      + '\n'
    (g,n) = graph.classify()
    result += r'0*\txt{g=%d,n=%d}' % (g,n) + '\n'
    result += r'\endxy'
    return result


## main

# parse command-line options
from optparse import OptionParser
parser = OptionParser(usage="""Usage: %prog [options] action [arg ...]

    Actions:

      vertices G N
        Print the vertex valences occurring in M_{g,n} graphs

      graphs V1,V2,...
        Print the graphs having only vertices of the specified valences.

      test
        Run internal code tests and report results.
        
    """)
parser.add_option("-n", "--silent",
                  action='store_true', dest='silent', default=False,
                  help="No output. (Use mainly for timing the algorithm.)")
parser.add_option("-L", "--latex",
                  action='store_true', dest='latex', default=False,
                  help="Print Xy-Pic code to draw graphs.")
parser.add_option("-o", "--output", dest="outfile", default=None,
                  help="Output file for `vertices` action.")
(options, args) = parser.parse_args()

if 0 == len(args) or 'help' == args[0]:
    parser.print_help()
    sys.exit(1)

elif 'test' == args[0]:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

elif 'vertices' == args[0]:
    del args[0]
    if len(args) < 2:
        parser.print_help()
        sys.exit(1)
    try:
        g = int(args[0])
        if g < 0:
            raise ValueError
    except ValueError:
        sys.stderr.write("Bad value '%s' for argument G: " \
                         "should be positive integer.\n" \
                         % (args[0],))
        sys.exit(1)
    try:
        n = positive_int(args[1])
    except ValueError, msg:
        sys.stderr.write("Bad value '%s' for argument N: " \
                         "should be non-negative integer.\n" \
                         % (args[1],))
        sys.exit(1)
    vvs = vertex_valences_for_given_g_and_n(g,n)
    for vv in vvs:
        print vv

elif 'graphs' == args[0]:
    # parse command line
    del args[0]
    if len(args) == 0:
        parser.print_help()
        sys.exit(1)
    valences_list = []
    for pattern in args:
        try:
            valences = [ positive_int(v) for v in pattern.split(',') ]
        except ValueError, msg:
            sys.stderr.write("Bad argument '%s': " \
                             "vertex valences must be a comma-separated " \
                             "list of positive integers: %s" \
                             % (pattern, msg))
            sys.exit(1)
        if sum(valences) % 2 != 0:
            sys.stderr.write("Invalid argument '%s': " \
                             "sum of vertex valences must be an even number. " \
                             "Aborting." \
                             % (pattern,))
            sys.exit(1)
        valences_list.append(valences)

    # open output file
    if options.outfile is None:
        outfile = sys.stdout
    else:
        outfile = open(options.outfile, 'w')
    # compute graphs matching given vertex sequences
    graphs = []
    for vertex_pattern in args:
        vertex_list = map(int, vertex_pattern.split(','))
        graphs += all_graphs(vertex_list)

    # output results
    if not options.silent:
        if options.latex:
            outfile.write(r"""
    \documentclass[a4paper,twocolumn]{article}
    \usepackage[curve,poly,xdvi]{xy}
    \begin{document}
    """)
        for graph in graphs:
            if options.latex:
                outfile.write(graph_to_xypic(graph)+'\n')
            else:
                outfile.write("%s\n" % ((graph,
                                         graph.genus(),
                                         graph.num_boundary_components(),
                                         graph.has_orientation_reversing_automorphism()),))
        outfile.write("\n")
        outfile.write("Found %d graphs.\n" % len(graphs))
        outfile.write("\n")
        if options.latex:
            outfile.write(r"\end{document}")
            outfile.write("\n")
