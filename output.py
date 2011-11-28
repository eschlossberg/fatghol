#! /usr/bin/env python
#
"""Classes and routines for writing Fatgraphs to files.
"""
__docformat__ = 'reStructuredText'


import iterators


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

% use Knuth's "concrete" fonts, which blend better with the
% typewriter font used for code listings
\usepackage[T1]{fontenc}
\usepackage[boldsans]{ccfonts}

\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{colortbl}
\usepackage{hyperref}
\usepackage{longtable}%
  \setcounter{LTchunksize}{100}%
\usepackage{multicol}
\usepackage{tensor}
\usepackage[usenames,dvipsnames]{xcolor}
% Xy-Pic v3.8 is needed for PDF support
\usepackage[pdf,color,curve,frame,line,poly]{xy}%
  \UseCrayolaColors

\newcommand{\corner}[3]{\ensuremath{\tensor[^{#2}]{#1}{^{#3}}}}
\newcommand{\cornerjoin}{\to}
% alternate:
%\newcommand{\corner}[3]{\ensuremath{\stackrel{#2}{-}#1\stackrel{#3}{-}}}
%\newcommand{\cornerjoin}{\kern-0.25ex\to\kern-0.25ex}

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
  \begin{minipage}{0.5\textwidth}
  {% left column: Xy-Pic diagram
""")
        self.write(LaTeXFile.render_to_xypic(graph))
        self.write(r"""
  }%
  \end{minipage}
  &% right column: Python code
  \begin{minipage}{0.4\textwidth}
""")
        self.write_repr(graph)
        self.write(r"""
  \end{minipage}
\end{tabular}
\vspace{-3em}
""")


    def write_automorphisms(self, automorphisms):
        """
        Output a table of graph automorphisms (if any).
        """
        if len(automorphisms) == 0:
            return
        graph = automorphisms[0].source
        vfmt = 'c' * graph.num_vertices
        efmt = 'c' * graph.num_edges
        bfmt = 'c' * graph.num_boundary_cycles
        self.write(r"""
\begin{center}
  \newcommand\Reven{\rowcolor{MidnightBlue!5}}
  \newcommand\Rodd{\rowcolor{MidnightBlue!25}}
  \begin{longtable}{l|%s|%s|%s}
""" % (vfmt, efmt, bfmt))
        # header
        self.write(r"$A_0$")
        self.write(" & ")
        self.write(str.join("&", [(r"{\slshape %s}" % LaTeXFile._vertex_label(v))
                                  for v in xrange(graph.num_vertices)]))
        self.write(" & ")
        self.write(str.join("&", [(r"{\slshape %s}" % str(graph.edge_numbering[e]))
                                  for e in xrange(graph.num_edges)]))
        self.write(" & ")
        self.write(str.join("&", [(r"$%s$" % LaTeXFile._BOUNDARY_CYCLES_LABELS[b])
                                  for b in xrange(graph.num_boundary_cycles)]))
        self.write(r"\\")
        self.write("\n")
        self.write(r"\hline")
        # one line per autormorphism
        nr = 0
        for a in automorphisms:
            if a.is_identity():
                continue # with next `a`
            nr += 1
            # mark rows with alternating colors
            if nr % 2 == 0:
                self.write(r"\Reven")
            else:
                self.write(r"\Rodd")
            # automorphisms are named A_n
            self.write(r"$A_{%d}$" % nr)
            if a.is_orientation_reversing():
                self.write(r"$\ast$")
            # write vertices permutation
            self.write(" & ")
            self.write(str.join("&", [(r"{%s}" % LaTeXFile._vertex_label(a.pv[v]))
                                      for v in xrange(graph.num_vertices)]))
            # edges permutation
            self.write(" & ")
            self.write(str.join("&", [(r"{%s}" % graph.edge_numbering[a.pe[e]])
                                      for e in xrange(graph.num_edges)]))
            # boundary cycles permutation
            self.write(" & ")
            self.write(str.join("&", [
                (r"$%s$" % LaTeXFile._BOUNDARY_CYCLES_LABELS[
                    graph.boundary_cycles.index(
                        a.transform_boundary_cycle(bcy))])
                 for bcy in graph.boundary_cycles ]))
            self.write(r"\\")
            self.write("\n")
        self.write(r"""
  \end{longtable}
\end{center}
""")

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
    
    def write_boundary_cycles(self, bcys):
        """
        Output a table of the given boundary cycles.
        """
        self.write(r"""
\begin{align*}
""")
        for nr, bcy in enumerate(sorted(bcys, key=LaTeXFile._sort_bcy_marking)):
            self.write(r"""%s &= (%s) \\ """ % (
                LaTeXFile._BOUNDARY_CYCLES_LABELS[nr],
                str.join(r'\cornerjoin',
                         [(r"\corner{%s}{%d}{%d}"
                           % (
                               LaTeXFile._vertex_label(corner[0]),
                               corner[1], # incoming
                               corner[2]  # outgoing
                               ))
                          for corner in bcy]),
                ))
        self.write(r"""
\end{align*}
""")


    def write_differential_start(self, title="Differentials", level=2, omit_null=True):
        """
        Open the "Differentials" section.
        """
        self._d_omit_null = omit_null
        self._d_head = (r"""
\%(sectioncmd)s{%(title)s}
""") % dict(title=title, sectioncmd=LaTeXFile._SECTIONING_COMMANDS[level])
        self._d_parts = [ ]
        self._d_maxterms = 0
        
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

    def write_differential_end(self):
        """
        Close the "Differentials" section.
        """
        if len(self._d_parts) > 0:
            self.write(self._d_head)
            if self._d_maxterms < 4:
                # use two-column layout
                self.write(r"""
\begin{multicols}{2}
""")
                self.write(str.join("\n\n", self._d_parts))
                self.write(r"""
\end{multicols}
""")
            else:
                # use 1-column layout
                self.write(r"""
\begin{longtable}[c]{l}
""")
                self.write(str.join(r"\\", self._d_parts))
                self.write(r"""
\end{longtable}
""")


    def write_markings(self, name, markings, per_row=8):
        """
        Output a table showing how the different markings of the same
        underlying fatgraph do number the boundary components.
        """
        graph = markings.graph
        n = graph.num_boundary_cycles
        N = len(markings.numberings)
        W = 1 + (per_row if N>per_row else N)
        cfmt = '|c' * (W-1)

        # make one table per each group of `per_row` markings
        self.write(r"""
\begin{center}
  \newcommand\Rhead{\rowcolor{Tan!25}}
  \newcommand\Reven{\rowcolor{Tan!10}}
  \newcommand\Rodd{\rowcolor{white}}
  \begin{longtable}{l%(cfmt)s}
    \multicolumn{%(width)d}{l}{} \endfirsthead
    \multicolumn{%(width)d}{l}{\em (continued.)} \endhead
""" % dict(cfmt=cfmt, width=W))
        done = 0
        for ms in iterators.chunks([per_row] * (N / per_row) + [N % per_row],
                                   markings.numberings):
            # if N % per_row == 0, we have an extra cycle; skip it
            if len(ms) == 0:
                break
            # table header lists graph/marking names
            if done > 0:
                self.write(r"""\pagebreak[0]""")
            self.write(r""" \Rhead""")
            for j in xrange(done, done + len(ms)):
                self.write(r""" & $%s^{(%d)}$""" % (name, j))
            self.write(r"""\\""")
            self.write(r"""\nopagebreak""")
            self.write('\n')
            # each row lists the marking of a certain boundary cycle
            for b in range(graph.num_boundary_cycles):
                if b == graph.num_boundary_cycles - 1:
                    self.write(r"""\nopagebreak""")
                if b % 2 == 0:
                    self.write(r"""\Reven""")
                else:
                    self.write(r"""\Rodd""")
                self.write("$%s$" % LaTeXFile._BOUNDARY_CYCLES_LABELS[b])
                self.write(r""" & """)
                self.write(str.join(" & ",
                                    (str(ms[j][b]) for j in range(len(ms)))))
                self.write(r""" \\""")
                self.write('\n')
            done += per_row
        self.write(r"""
  \end{longtable}
\end{center}
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
        vertex_reprs = [ str(V) for V in graph.vertices ]
        m = max(len(vtxt) for vtxt in vertex_reprs)
        self.write(
            str.join("\n", [
                str.join("", [
                    # initial indent
                    "  ",
                    # vertex representation
                    vtxt,
                    # padding to align `#`'s
                    (' ' * (m - len(vtxt) + 1)),
                    # vertex label
                    "# ",
                    LaTeXFile._vertex_label(v),
                    ])
                for v,vtxt in enumerate(vertex_reprs) ]))
        self.write("""\n])\n""")
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
                          % (k,
                             LaTeXFile._ensure_large_odd(len(graph.vertices[k-1])),
                             k,
                             LaTeXFile._rotation_angle(K,k)))

        for l in xrange(graph.num_edges):
            ((v1, i1), (v2, i2)) = graph.edges[l].endpoints
            if v1 != v2:
                result.append((r'"v%d"*+[o]\txt{{%s}}*\frm{o};"v%d"*+[o]\txt{{%s}}*\frm{o}**\crv{"v%dl%d"&"v%dl%d"}?(.75)*\txt{\colorbox{CadetBlue!25}{\scriptsize\bfseries %d}},'
                               % (v1+1, LaTeXFile._vertex_label(v1),
                                  v2+1, LaTeXFile._vertex_label(v2),
                                  v1+1, 1+graph.vertices[v1].index(l),
                                  v2+1, 1+graph.vertices[v2].index(l),
                                  graph.edge_numbering[l])))
            else:
                h = graph.vertices[v1].index(l)
                result.append((r'"v%d"*+[o]\txt{{%s}}*\frm{o};"v%d"*+[o]\txt{{%s}}*\frm{o}**\crv{"v%dl%d"&"v%dl%d"}?(.75)*\txt{\colorbox{CadetBlue!25}{\scriptsize\bfseries %d}},'
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
