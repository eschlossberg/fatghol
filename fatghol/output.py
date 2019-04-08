#! /usr/bin/env python
#
"""Classes and routines for writing Fatgraphs to files.
"""
#
#   Copyright (C) 2008-2012 Riccardo Murri <riccardo.murri@gmail.com>
#   All rights reserved.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
__docformat__ = 'reStructuredText'


import datetime
import os.path
import tempita

from fatghol.combinatorics import factorial
import iterators


class LaTeXOutput(object):

    def __init__(self, stream, **kw):
        self._output = stream
        self.data = tempita.bunch(
            filename=self._output.name,
            today=datetime.date.today(),
            **kw
            )
        try:
            import pkg_resources
            template_txt = pkg_resources.resource_string('fatghol', 'template.tex')
        except ImportError:
            # try to open in current directory
            if os.path.exists('template.tex') and os.path.isfile('template.tex'):
                with open('template.tex', 'r') as template_file:
                    template_txt = template_file.read()
            else:
                raise
        self._template = tempita.Template(
            template_txt,
            delimeters=('<<', '>>'),
            name='fatghol/template.tex',
            namespace=self.data,
            )
        self._total_graphs = 0
        self._total_marked_graphs = 0
        # document appendices
        self._appendix_graph_markings = None


    def close(self, **kw):
        are_there_appendices = (self._appendix_graph_markings is not None)
        self._output.write(self._template.substitute(
            appendices=are_there_appendices,
            appendix_graph_markings=self._appendix_graph_markings,
            total_num_graphs=self._total_graphs,
            total_num_marked_graphs=self._total_marked_graphs,
            **kw))
        self._output.close()


    def start_section(self, num_edges, num_vertices, intro=""):
        """
        Start a new section with the given title text.

        The `intro` text is printed between the title text and the
        ensuing list content.
        """
        if 'sections' not in self.data:
            self.data.sections = list()
        title = ("Fatgraphs with %(edges)s / %(vertices)s" % dict(
            edges = ("%d edges" % num_edges) if num_edges > 1 else "1 edge",
            vertices = ("%d vertices" % num_vertices) if num_vertices > 1 else "1 vertex",
            ))
        self._current_section = tempita.bunch(
            title=title,
            intro=intro,
            graphs=list(),
            )
        self._section_graphs = 0
        self._section_marked_graphs = 0
        self._section_orientable_marked_graphs = 0

    def end_section(self):
        self._current_section.num_graphs = self._section_graphs
        self._current_section.num_marked_graphs = self._section_marked_graphs
        self._current_section.num_orientable_marked_graphs = self._section_orientable_marked_graphs
        self._current_section.num_nonorientable_marked_graphs = self._section_marked_graphs - self._section_orientable_marked_graphs
        self.data.sections.append(self._current_section)


    def start_graph(self, name, graph, pool):
        """
        Start constructing the section relative to the given graph in the output.
        """
        assert num_orientable_markings >= 0
        num_orientable_markings = sum((1 if NG.is_oriented() else 0) for NG in pool)
        self._current_graph = graph
        self._current_marked_graph = pool[0]
        self._current_marked_graph_automorphisms = list(self._current_marked_graph.automorphisms())
        self._total_graphs += 1
        self._total_marked_graphs += len(pool)
        self._section_graphs += 1
        self._section_marked_graphs += len(pool)
        self._section_orientable_marked_graphs += num_orientable_markings
        self._current_graph_data = tempita.bunch(
            name=name,
            orientable=graph.is_oriented(),
            num_orientable_markings=num_orientable_markings,
            latex_repr=LaTeXOutput.render_to_xypic(graph, cross_out=False),
            python_repr=LaTeXOutput.render_to_python(graph),
            boundary_cycles=LaTeXOutput._fmt_boundary_cycles(graph.boundary_cycles)
            )


    def end_graph(self):
        """
        Output the last constructed graph.

        Any call to `add_*` intervening before the next `start_graph`
        will be discarded.
        """
        self._current_section.graphs.append(self._current_graph_data)


    def add_automorphisms(self, automorphisms):
        """
        Output a table of graph automorphisms (if any).
        """
        if len(automorphisms) == 0:
            return
        assert self._current_graph == automorphisms[0].source
        graph = self._current_graph
        vfmt = 'c' * graph.num_vertices
        efmt = 'c' * graph.num_edges
        bfmt = 'c' * graph.num_boundary_cycles
        parts = [ ]
        parts.append(r"""
\begin{center}
  \newcommand\Reven{\rowcolor{MidnightBlue!5}}
  \newcommand\Rodd{\rowcolor{MidnightBlue!25}}
  \begin{longtable}{l|%s|%s|%s}
""" % (vfmt, efmt, bfmt))
        # header
        parts.append(r"$A_0$")
        parts.append(" & ")
        parts.append(str.join("&", [(r"{%s}" % LaTeXOutput._vertex_label(v))
                                  for v in xrange(graph.num_vertices)]))
        parts.append(" & ")
        parts.append(str.join("&", [(r"{%s}" % str(graph.edge_numbering[e]))
                                  for e in xrange(graph.num_edges)]))
        parts.append(" & ")
        parts.append(str.join("&", [(r"$%s$" % LaTeXOutput._BOUNDARY_CYCLES_LABELS[b])
                                  for b in xrange(graph.num_boundary_cycles)]))
        parts.append(r"\\")
        parts.append("\n")
        parts.append(r"\hline")
        # one line per automorphism
        nr = 0
        for a in automorphisms:
            if a.is_identity():
                continue # with next `a`
            nr += 1
            # mark rows with alternating colors
            if nr % 2 == 0:
                parts.append(r"\Reven")
            else:
                parts.append(r"\Rodd")
            # automorphisms are named A_n
            parts.append(r"$A_{%d}$" % nr)
            if a.is_orientation_reversing():
                parts.append(r"${}^\dag$")
            if a not in self._current_marked_graph_automorphisms:
                parts.append(r"${}^\ddag$")
            # write vertices permutation
            parts.append(" & ")
            parts.append(str.join("&", [(r"{%s}" % LaTeXOutput._vertex_label(a.pv[v]))
                                      for v in xrange(graph.num_vertices)]))
            # edges permutation
            parts.append(" & ")
            parts.append(str.join("&", [(r"{%s}" % graph.edge_numbering[a.pe[e]])
                                      for e in xrange(graph.num_edges)]))
            # boundary cycles permutation
            parts.append(" & ")
            parts.append(str.join("&", [
                (r"$%s$" % LaTeXOutput._BOUNDARY_CYCLES_LABELS[
                    graph.boundary_cycles.index(
                        a.transform_boundary_cycle(bcy))])
                 for bcy in graph.boundary_cycles ]))
            parts.append(r"\\")
            parts.append("\n")
        parts.append(r"""
  \end{longtable}
\end{center}
""")
        # add it to the output
        self._current_graph_data.automorphisms = str.join('', parts)
        

    def add_differential_start(self, omit_null=True):
        """
        Start writing the "Differentials" section.
        """
        self._d_omit_null = omit_null
        self._d_parts = [ ]
        self._d_maxterms = 0
        
    def add_differential(self, D, j, name, labelfn=None):
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
        parts = [ r"$", (r"D(%s) = " % name) ]
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
            parts.append(" %s%s" % (coeff, labelfn(D, i)))
            cnt += 1
        if cnt == 0:
            parts.append("0")
        parts.append(r"$")
        if cnt != 0 or (cnt == 0 and not self._d_omit_null):
            self._d_parts.append(str.join("", parts))
            self._d_maxterms = max(cnt, self._d_maxterms)

    def add_differential_end(self):
        """
        Close the "Differentials" section.
        """
        if len(self._d_parts) > 0:
            parts = [ ]
            if self._d_maxterms < 4:
                # use two-column layout
                parts.append(r"""
\begin{multicols}{2}
""")
                parts.append(str.join("\n\n", self._d_parts))
                parts.append(r"""
\end{multicols}
""")
            else:
                # use 1-column layout
                parts.append(r"""
\begin{longtable}[c]{l}
""")
                parts.append(str.join(r"\\", self._d_parts))
                parts.append(r"""
\end{longtable}
""")
            # add it to the current graph
            self._current_graph_data.differentials = str.join('', parts)


    def add_markings(self, name, markings, per_row=8):
        """
        Output a table showing how the different markings of the same
        underlying fatgraph do number the boundary components.
        """
        graph = markings.graph
        n = graph.num_boundary_cycles
        N = len(markings.numberings)
        self._section_marked_graphs += N
        if N == factorial(n):
            self._current_graph_data.markings = (r"""
Fatgraph $%(name)s$ only has the identity automorphism, so the
marked fatgraphs $%(name)s^{(0)}$ to $%(name)s^{(%(num_markings)d)}$
are formed by decorating boundary cycles of $%(name)s$ with
all permutations of $(%(basetuple)s)$ in lexicographic order.
See Section ``Markings of fatgraphs with trivial automorphisms''
for a complete table.
""" % dict(name=name, num_markings=N,
           basetuple=(str.join(',', [str(x) for x in xrange(n)]))))
            if self._appendix_graph_markings is None:
                # only need to do this once per document
                self._appendix_graph_markings = self._fmt_numberings(
                        markings.numberings, 'G', n, per_row)
        else:
            self._current_graph_data.markings = self._fmt_numberings(
                markings.numberings, name, n, per_row)


    ##
    ## helper methods
    ##
    
    _BOUNDARY_CYCLES_LABELS = [
        r'\alpha',
        r'\beta',
        r'\gamma',
        r'\delta',
        r'\epsilon',
        r'\zeta',
        r'\eta',
        r'\theta',
        r'\iota',
        r'\kappa',
        r'\lambda',
        r'\mu',
        r'\nu',
        r'\xi',
        #r'$o$' # not visually distinct from 'o' or '0'
        r'\pi',
        r'\rho',
        r'\sigma',
        r'\tau',
        r'\upsilon',
        r'\phi',
        r'\chi',
        r'\psi',
        r'\omega'
        ]

    @staticmethod
    def _fmt_boundary_cycles(bcys):
        """
        Output a table of the given boundary cycles.
        """
        parts = [ ]
        parts.append(r"""
\begin{align*}
""")
        for nr, bcy in enumerate(sorted(bcys, key=LaTeXOutput._sort_bcy_marking)):
            parts.append(r"""%s &= (%s) \\ """ % (
                LaTeXOutput._BOUNDARY_CYCLES_LABELS[nr],
                str.join(r'\cornerjoin',
                         [(r"\corner{%s}{%d}{%d}"
                           % (
                               LaTeXOutput._vertex_label(corner[0]),
                               corner[1], # incoming
                               corner[2]  # outgoing
                               ))
                          for corner in bcy]),
                ))
        parts.append(r"""
\end{align*}
""")
        return str.join('', parts)


    @staticmethod
    def _fmt_numberings(numberings, name, num_boundary_cycles, per_row):
        N = len(numberings)
        W = 1 + (per_row if N>per_row else N)
        cfmt = '|c' * (W-1)
        parts = [ ]
        # make one table per each group of `per_row` markings
        parts.append(r"""
\begin{center}
  \newcommand\Rhead{\rowcolor{Tan!25}}
  \newcommand\Reven{\rowcolor{Tan!10}}
  \newcommand\Rodd{\rowcolor{white}}
  \begin{longtable}{l%(cfmt)s}
    \multicolumn{%(width)d}{l}{} \endfirsthead
    \multicolumn{%(width)d}{l}{\em (continued.)} \endhead
""" % dict(cfmt=cfmt, width=W))
        done = 0
        for ms in iterators.chunks([per_row] * (N / per_row) + [N % per_row], numberings):
            # if N % per_row == 0, we have an extra cycle; skip it
            if len(ms) == 0:
                break
            # table header lists graph/marking names
            if done > 0:
                parts.append(r"""\pagebreak[0]""")
            parts.append(r""" \Rhead""")
            for j in xrange(done, done + len(ms)):
                parts.append(r""" & $%s^{(%d)}$""" % (name, j))
            parts.append(r"""\\""")
            parts.append(r"""\nopagebreak""")
            parts.append('\n')
            # each row lists the marking of a certain boundary cycle
            for b in range(num_boundary_cycles):
                if b == num_boundary_cycles - 1:
                    parts.append(r"""\nopagebreak""")
                if b % 2 == 0:
                    parts.append(r"""\Reven""")
                else:
                    parts.append(r"""\Rodd""")
                parts.append("$%s$" % LaTeXOutput._BOUNDARY_CYCLES_LABELS[b])
                parts.append(r""" & """)
                parts.append(str.join(" & ",
                                    (str(ms[j][b]) for j in range(len(ms)))))
                parts.append(r""" \\""")
                parts.append('\n')
            done += per_row
        parts.append(r"""
  \end{longtable}
\end{center}
""")
        return str.join('', parts)


    @staticmethod
    def render_to_python(graph):
        """
        Return pretty-printed Python representation of `graph`.
        """
        assert isinstance(graph, Fatgraph)
        vertex_reprs = [ str(V) for V in graph.vertices ]
        m = max(len(vtxt) for vtxt in vertex_reprs)
        return str.join('', [
            """Fatgraph([\n""",
            str.join("\n", [
                str.join("", [
                    # initial indent
                    "  ",
                    # vertex representation
                    vtxt, ",",
                    # padding to align `#`'s
                    (' ' * (m - len(vtxt))),
                    # vertex label
                    "# ",
                    LaTeXOutput._vertex_label(v),
                    ])
                for v,vtxt in enumerate(vertex_reprs) ]),
            """\n])\n""",
            ])
    
            
    @staticmethod
    def render_to_xypic(graph, cross_out=True):
        """Return XY-Pic code snippet to render graph `graph`."""
        assert isinstance(graph, Fatgraph)

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
                          % (k,
                             LaTeXOutput._ensure_large_odd(len(graph.vertices[k-1])),
                             k,
                             LaTeXOutput._rotation_angle(K,k)))

        for l in xrange(graph.num_edges):
            ((v1, i1), (v2, i2)) = graph.edges[l].endpoints
            if v1 != v2:
                result.append((r'"v%d"*+[o]\txt{{%s}}*\frm{o};"v%d"*+[o]\txt{{%s}}*\frm{o}**\crv{"v%dl%d"&"v%dl%d"}?(.75)*\txt{\colorbox{CadetBlue!25}{\scriptsize\bfseries %d}},'
                               % (v1+1, LaTeXOutput._vertex_label(v1),
                                  v2+1, LaTeXOutput._vertex_label(v2),
                                  v1+1, 1+graph.vertices[v1].index(l),
                                  v2+1, 1+graph.vertices[v2].index(l),
                                  graph.edge_numbering[l])))
            else:
                h = graph.vertices[v1].index(l)
                result.append((r'"v%d"*+[o]\txt{{%s}}*\frm{o};"v%d"*+[o]\txt{{%s}}*\frm{o}**\crv{"v%dl%d"&"v%dl%d"}?(.75)*\txt{\colorbox{CadetBlue!25}{\scriptsize\bfseries %d}},'
                               % (v1+1, LaTeXOutput._vertex_label(v1),
                                  v2+1, LaTeXOutput._vertex_label(v2),
                                  v1+1, h+1, v2+1, 1+graph.vertices[v1].index(l,h+1),
                                  graph.edge_numbering[l])))

        # cross-out graph, if not orientable
        if cross_out and not graph.is_oriented():
            result.append(r"0,(-1.00,+1.00);(+1.00,-1.00)**[red]@{..},")
            result.append(r"0,(-1.00,-1.00);(+1.00,+1.00)**[red]@{..},")

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
            return min(nr)
        except TypeError: # `nr` not iterable
            return nr

    @staticmethod
    def _ensure_large_odd(nr):
        if nr % 2 == 0:
            return 1+nr
        else:
            return 2+nr


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="output",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
