#! /usr/bin/env python
#
"""Classes and routines for writing Fatgraphs to files.
"""
__docformat__ = 'reStructuredText'


class LaTeXFile(file):
    """
    
    """

    def __new__(cls, path, buffering=1):
        return file.__new__(cls, path, 'w', buffering)


    def __init__(self, path, buffering=1):
        file.__init__(self, path, 'w', buffering)
        self.write(r"""
\documentclass[a4paper]{article}
\raggedbottom

\usepackage{amsmath}
\usepackage{colortbl}
\usepackage{multirow}
\usepackage[usenames,dvipsnames]{xcolor}
\usepackage[xdvi,color,curve,frame,line,poly]{xy}%
  \UseCrayolaColors

\begin{document}
""")
        self._section_level = -1



    def close(self, final_text=""):
        self._close_pending_sections()
        self.write(final_text)
        self.write("\n")
        self.write(r"""
\end{document}
""")
        file.close(self)


    _SECTIONING_COMMANDS = [
        'section*',       # level=0
        'subsection*',
        'subsubsection*',
        'paragraph',
        'subparagraph',
        ]

    def section(self, title, level=0, intro=""):
        """
        Start a new section with the given title text.

        The `intro` text is printed between the title text and the
        ensuing list content.
        """
        assert level in range(len(LaTeXFile._SECTIONING_COMMANDS)), \
               ("Invalid value '%s' for `level`:"
                " must be an integer in range 0 to %d"
                % (level, len(LaTeXFile._SECTIONING_COMMANDS)))
        self._close_pending_sections()
        if level < 2:
            self.write(r"\vspace{2em}")
            self.write("\n")
        self.write(r"""

  \%(sectioncmd)s{%(title)s}%%
  %(intro)s%%

""" % {
                       'sectioncmd': LaTeXFile._SECTIONING_COMMANDS[level],
                       'title':title,
                       'intro':intro,
                       })
        self._section_level = level

    def _close_pending_sections(self):
        if self._section_level > -1:
            # close the section/multicols
            self.write(r"""
""")
        

    def write_graph(self, graph, name):
        """
        Output the given graph to this file.
        """
        self.write(r"""

\vspace{-1em}
\begin{tabular}{lr}
  \begin{minipage}{0.6\textwidth}
  {% left column: Xy-Pic diagram
""")
        self.write(LaTeXFile.render_to_xypic(graph))
        self.write(r"""
  }%
  \end{minipage}
  &% right column: Python code
  \begin{minipage}{0.3\textwidth}
""")
        self.write_repr(graph)
        self.write(r"""
  \end{minipage}
\end{tabular}
\vspace{-3em}
""")


    def write_automorphisms(self, graph, Aut=None):
        """
        Output a table of graph automorphisms (if any).
        """
        if Aut is None:
            Aut = list(graph.automorphisms())
        vfmt = 'c' * graph.num_vertices
        efmt = 'c' * graph.num_edges
        self.write(r"""
\begin{center}
  \newcommand\C{\cellcolor{gray!20}}
  \begin{tabular}{|c|%s|%s|}
  \hline
""" % (vfmt, efmt))
        nr = 0
        for a in Aut:
            if LaTeXFile._is_identity(a):
                continue # with next `a`
            # use uppercase latin letter for automorphism
            nr += 1
            self.write(r"\multirow{2}{1em}{%s}" % (chr(64 + nr)))
            self.write(" & ")
            self.write(str.join("&", [LaTeXFile._vertex_label(v)
                                      for v in xrange(graph.num_vertices)]))
            self.write(" & ")
            self.write(str.join("&", [str(graph.edge_numbering[e])
                                      for e in xrange(graph.num_edges)]))
            self.write(r"\\")
            self.write("\n")
            #self.write(r"\rowcolor{gray!20}")
            self.write(" & ")
            self.write(str.join("&", [(r"\C{%s}" % LaTeXFile._vertex_label(a.pv[v]))
                                      for v in xrange(graph.num_vertices)]))
            self.write(" & ")
            self.write(str.join("&", [(r"\C{%s}" % graph.edge_numbering[a.pe[e]])
                                      for e in xrange(graph.num_edges)]))
            self.write(r"\\")
            self.write(r"\hline")
            self.write("\n")
        self.write(r"""
  \end{tabular}
\end{center}
""")

    @staticmethod
    def _is_identity(a):
        for v in xrange(len(a.pv)):
            if v != a.pv[v]:
                return False
        for e in xrange(len(a.pe)):
            if e != a.pe[e]:
                return False
        return True


    def write_boundary_cycles(self, bcys):
        """
        Output a table of the given boundary cycles.
        """
        self.write(r"""
\begin{tabular}{rl}
""")
        for (bcy, nr) in sorted(bcys.iteritems(),
                                key=LaTeXFile._sort_bcy_marking):
            self.write(r"""\textsl{%s} & (%s) \\ """ % (
                LaTeXFile._fmt_bcy_marking(nr),
                str.join(",", [("(%s,%d,%d)"
                                % (LaTeXFile._vertex_label(corner[0]),
                                   corner[1], corner[2]))
                                for corner in bcy]),
                ))
            self.write('\n')
        self.write(r"""
\end{tabular}
""")


    def write_differential(self, D, j, name, labelfn=None):
        """
        Output the expansion of the `j`-th column of differential
        operator `D` as a relation between graphs.

        The string `name` is used to name the fatgraph at LHS;
        function `labelfn` is used to generate names for the graphs
        appearing at RHS.  Argument `labelfn` should return a string
        name of the graph at the `i`-th row of `D`; it is passed `D`
        and `i`.

        :param SimpleMatrix D:
        :param int i:
        :param str name:
        :param function labelfn:
        """
        self.write(r"""
\begin{equation*}
  D(%s) = """ % name)
        cnt = 0
        for i in xrange(D.num_rows):
            coeff = D.getEntry(i, j)
            if coeff == 0:
                continue # with next graph
            elif coeff == +1:
                coeff = "+"
            elif coeff == -1:
                coeff = "-"
            else:
                coeff = "%+d" % coeff
            self.write(" %s%s" % (coeff, labelfn(D, i)))
            cnt += 1
        if cnt == 0:
            self.write("0")
        self.write(r"""
\end{equation*}
        """)


    def write_repr(self, graph):
        """
        Pretty-print the representation of `graph`.
        """
        self.write(r"""
\begin{verbatim}
""")
        assert isinstance(graph, Fatgraph)
        self.write("""Fatgraph([\n""")
        for V in graph.vertices:
            self.write("  %s,\n" % V)
        self.write("""])\n""")
        self.write(r"""
\end{verbatim}
""")
            
    ##
    ## helper methods
    ##
    
    @staticmethod
    def render_to_xypic(graph, cross_out=True):
        """Return XY-Pic code snippet to render graph `graph`."""

        g = graph.genus
        n = graph.num_boundary_cycles

        result = [ ]

        K = graph.num_vertices
        result.append(r'\xy 0;<0.25\linewidth,0cm>:0,')

        # put graph vertices on a regular polygon
        if K == 1:
            result.append('(-1,0)="v1",')
        if K == 2:
            result.append('(-1,0)="v1",%\n(1,0)="v2",')
        else:
            result.append((r'{\xypolygon%d"v"{~={0}~>{}}},'
                           % (max(K,3))))

        # mark invisible "control points" for bezier curves connecting vertices
        for k in range(1,K+1):
            result.append(r'"v%d",{\xypolygon%d"v%dl"{~:{(1.20,0):}~={%d}~>{}}},'
                          % (k, 2*len(graph.vertices[k-1])-2, k,
                             LaTeXFile._rotation_angle(K,k)))

        for l in xrange(graph.num_edges):
            ((v1, i1), (v2, i2)) = graph.edges[l].endpoints
            if v1 != v2:
                result.append((r'"v%d"*+[o]\hbox{{%s}}*\frm{o};"v%d"*+\txt{{%s}}**\crv{"v%dl%d"&"v%dl%d"}?(.75)*\txt{\colorbox{gray!25}{\small\bfseries %d}},%%?(.3)*\dir{>},'
                               % (v1+1, LaTeXFile._vertex_label(v1),
                                  v2+1, LaTeXFile._vertex_label(v2),
                                  v1+1, 1+graph.vertices[v1].index(l),
                                  v2+1, 1+graph.vertices[v2].index(l),
                                  graph.edge_numbering[l])))
            else:
                h = graph.vertices[v1].index(l)
                result.append((r'"v%d"*+[o]\hbox{{%s}}*\frm{o};"v%d"*+\txt{%s}**\crv{"v%dl%d"&"v%dl%d"}?(.75)*\txt{\colorbox{gray!25}{\small\bfseries %d}},%%?(.3)*\dir{>},'
                               % (v1+1, LaTeXFile._vertex_label(v1),
                                  v2+1, LaTeXFile._vertex_label(v2),
                                  v1+1, h+1, v2+1, 1+graph.vertices[v1].index(l,h+1),
                                  graph.edge_numbering[l])))

        # cross-out graph, if not orientable
        if cross_out and not graph.is_oriented():
            result.append(r"0,(-1.00,+1.00);(+1.00,-1.00)**[red][|(10)]@{..},")
            result.append(r"0,(-1.00,-1.00);(+1.00,+1.00)**[red][|(10)]@{..},")

        result.append(r'\endxy')

        return str.join('%\n', result)


    @staticmethod
    def _vertex_label(v):
        # vertices are labeled with lowercase latin letters
        return chr(97 + v)

    @staticmethod
    def _rotation_angle(K,k):
        return 90+(k-1)*360/K

    @staticmethod
    def _fmt_bcy_marking(nr):
        try:
            return str.join(",", [str(elt) for elt in nr])
        except TypeError: # `nr` not iterable
            return str(nr)

    @staticmethod
    def _sort_bcy_marking(nr):
        try:
            return min(x)
        except TypeError: # `nr` not iterable
            return x


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="output",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
