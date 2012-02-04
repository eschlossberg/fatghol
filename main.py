#! /usr/bin/env python
"""Command-line front-end for `rg.py`.
"""
__version__ = '5.4'
__docformat__ = 'reStructuredText'


## stdlib imports

import sys
# require Python 2.6 (for the "fractions" module)
if sys.version < '2.6.0':
    sys.stderr.write("Program %s requires Python at least version 2.6,"
                     " but this Python interpreter is version %s. Aborting."
                     % (sys.argv[0], sys.version.split()[0]))
    sys.exit(1)

import cython

from collections import defaultdict
from fractions import Fraction
import gc
import logging
import os
import os.path
import resource
import tempfile


## application-local imports

from const import euler_characteristics, orbifold_euler_characteristics
from combinatorics import minus_one_exp
from graph_homology import FatgraphComplex, NumberedFatgraphPool
from loadsave import load
from rg import (
    Fatgraph,
    MgnGraphsIterator,
    )
from runtime import runtime
from simplematrix import SimpleMatrix
import timing
from utils import concat, positive_int
from valences import vertex_valences_for_given_g_and_n



## actions

def compute_graphs(g,n):
    """
    Compute Fatgraphs occurring in `M_{g,n}`.

    Return a pair `(graphs, D)`, where `graphs` is the list of
    `(g,n)`-graphs, and `D` is a list, the `k`-th element of which is
    the list of differentials of graphs with `k` edges.
    """
    timing.start("compute_graphs(%d,%d)" % (g,n))

    runtime.g = g
    runtime.n = n
    
    logging.info("Stage I:"
                 " Computing fat graphs for g=%d, n=%d ...",
                 g, n)
    G = FatgraphComplex(g,n)
    
    logging.info("Stage II:"
                 " Computing matrix form of boundary operators D[1],...,D[%d] ...",
                 G.length-1)
    D = G.compute_boundary_operators()

    timing.stop("compute_graphs(%d,%d)" % (g,n))
    return (G, D)

    
def compute_homology(g, n):
    """
    Compute homology ranks of the graph complex of `M_{g,n}`.

    Return array of homology ranks.
    """
    timing.start("compute_homology(%d,%d)" % (g,n))

    runtime.g = g
    runtime.n = n

    (G, D) = compute_graphs(g, n)

    logging.info("Stage III: Computing rank of homology modules ...")
    hs = list(reversed(D.compute_homology_ranks()))

    timing.stop("compute_homology(%d,%d)" % (g,n))

    # compare orbifold Euler characteristics
    computed_chi = G.orbifold_euler_characteristics
    logging.info("Computed orbifold Euler characteristics: %s", computed_chi)
    expected_chi = orbifold_euler_characteristics(g,n)
    logging.info("  Expected orbifold Euler characteristics (according to Harer): %s",
                 expected_chi)
    if computed_chi != expected_chi:
        logging.error("Expected and computed orbifold Euler characteristics do not match!"
                      " (computed: %s, expected: %s)" % (computed_chi, expected_chi))

    # compare Euler characteristics
    computed_e = 0
    for i in xrange(len(hs)):
        computed_e += minus_one_exp(i)*hs[i]
    expected_e = euler_characteristics(g,n)
    logging.info("Computed Euler characteristics: %s" % computed_e)
    logging.info("  Expected Euler characteristics: %s" % expected_e)
    if computed_e != expected_e:
        logging.error("Computed and expected Euler characteristics do not match:"
                      " %s vs %s" % (computed_e, expected_e))

    # verify result against other known theorems
    if g>0:
        # from Harer's SLN1337, Theorem 7.1
        if hs[1] != 0:
            logging.error("Harer's Theorem 7.1 requires h_1=0 when g>0")
        ## DISABLED 2009-03-27: Harer's statement seems to be incorrect,
        ## at least for low genus...
        ## # From Harer's SLN1337, Theorem 7.2 
        ## if g==1 or g==2:
        ##     if hs[2] != n:
        ##         logging.error("Harer's Theorem 7.2 requires h_2=%d when g=1 or g=2" % n)
        ## elif g>2:
        ##     if hs[2] != n+1:
        ##         logging.error("Harer's Theorem 7.2 requires h_2=%d when g>2" % n+1)

    return hs


def compute_valences(g,n):
    """
    Compute vertex valences occurring in g,n Fatgraphs.

    Return list of valences.
    """
    runtime.g = g
    runtime.n = n
    logging.info("Computing vertex valences occurring in g=%d,n=%d fatgraphs ...", g, n)
    return vertex_valences_for_given_g_and_n(g,n)



## main

# disable core dumps
resource.setrlimit(resource.RLIMIT_CORE, (0,0))

# parse command-line options
import argparse
parser = argparse.ArgumentParser(
    description="""
Actions:

  graphs G N
    Generate the graphs occurring in M_{g,n}.

  homology G N
    Print homology ranks of M_{g,n}.

  latex G N [-s DIR] [-o FILE]
    Read the listings of M_{g,n} fatgraphs (from directory DIR)
    and output a pretty-print catalogue of the graphs as LaTeX documents.
  
  valences G N
    Print the vertex valences occurring in M_{g,n} graphs.

  shell
    Start an interactive PyDB shell.
  
  selftest
    Run internal code tests and report failures.
    """,
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument('action', metavar='ACTION', default='help',
                    help="Action to perform, see above.")
parser.add_argument('args', metavar='ARG', nargs='*',
                    help="Arguments depend on the actual action, see above.")
# option arguments
if not cython.compiled:
    parser.add_argument("-D", "--debug", nargs='?',
                        dest="debug", default=None, const='debug',
                        help="""Enable debug features:
* pydb -- run Python debugger if an error occurs
* profile -- dump profiler statistics in a .pf file.
Several features may be enabled by separating them
with a comma, as in '-D pydb,profile'.""")
parser.add_argument("-l", "--logfile",
                    action='store', dest='logfile', default=None,
                    help="""Redirect log messages to the named file
(by default log messages are output to STDERR).""")
parser.add_argument("-o", "--output", dest="outfile", default=None,
                    help="Save results into named file.")
parser.add_argument("-s", "--checkpoint", dest="checkpoint_dir", default=None,
                    help="Directory for saving computation state.")
parser.add_argument("-u", "--afresh", dest="restart", action="store_false", default=True,
                    help="Do NOT restart computation from the saved state in checkpoint directory.")
parser.add_argument("-v", "--verbose",
                    action="count", dest="verbose", default=0,
                    help="Print informational and status messages as the computation goes on.")
parser.add_argument('-V', '--version', action='version', version=__version__)
cmdline = parser.parse_args()

# make options available to loaded modules
runtime.options = cmdline

# print usage message if no args given
if 'help' == cmdline.action:
    parser.print_help()
    sys.exit(0)

# configure logging
if cmdline.logfile is None:
    log_output = sys.stderr
else:
    log_output = file(cmdline.logfile, 'a')

if cmdline.verbose == 0:
    log_level = logging.ERROR
elif cmdline.verbose == 1:
    log_level = logging.INFO
else:
    log_level = logging.DEBUG

logging.basicConfig(level=log_level,
                    stream=log_output,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S")
    
# ensure the proper optimization level is selected
if __debug__ and cmdline.debug is None:
    try:
        import os
        os.execl(sys.executable, *([sys.executable, '-O'] + sys.argv))
    finally:
        logging.warning("Could not execute '%s', ignoring '-O' option." 
                        % str.join(" ", [sys.executable, '-O'] + sys.argv))

# enable optional features
if not cython.compiled and cmdline.debug is not None:
    debug = cmdline.debug.split(",")
    if 'pydb' in debug:
        try:
            import pydb
            sys.excepthook = pydb.exception_hook
            logging.debug("PyDB enabled: exceptions will start a debugging session.")
        except ImportError:
            logging.warning("Could not import 'pydb' module - PyDB not enabled.")
    if 'profile' in debug:
        try:
            import hotshot
            pf = hotshot.Profile(__name__ + '.pf')
            pf.start()
            logging.debug("Started call profiling with 'hotshot' module.")
        except ImportError:
            logging.warning("Could not import 'hotshot' - call profiling *not* enabled.")


# hack to allow 'N1,N2,...' or 'N1 N2 ...' syntaxes
for (i, arg) in enumerate(cmdline.args):
    if arg.find(","):
        if arg.startswith('M'):
            cmdline.args[i:i+1] = arg[1:].split(",")
        else:
            cmdline.args[i:i+1] = arg.split(",")

# open output file
if cmdline.outfile is None:
    outfile = sys.stdout
else:
    outfile = open(cmdline.outfile, 'w')


# shell -- start interactive debugging shell
if 'shell' == cmdline.action:
    if cython.compiled:
        logging.error("The 'shell' command is not available when compiled.")
        sys.exit(1)
    else:
        try:
            import pydb
        except ImportError:
            sys.stderr.write("ERROR: Could not import 'pydb' module - Aborting.\n")
            logging.warning("Could not import 'pydb' module - Aborting.")
            sys.exit(1)

        print ("""Starting interactive session (with PyDB %s).

        Any Python expression may be evaluated at the prompt.
        All symbols from modules `rg`, `homology`, `graph_homology`
        have already been imported into the main namespace.

        """ % pydb.__version__)
        pydb.debugger([
            "from homology import *",
            "from graph_homology import *",
            "from rg import *",
            ])
        
# selftest -- run doctests and acceptance tests on simple cases
elif 'selftest' == cmdline.action:
    failures = 0
    
    import doctest
    import imp
    for module in [ 'rg',
                    'homology',
                    'graph_homology',
                    'combinatorics',
                    'iterators',
                    'cyclicseq',
                    ]:
        try:
            file, pathname, description = imp.find_module(module)
            m = imp.load_module(module, file, pathname, description)
            logging.debug("Running Python doctest on '%s' module..." % module)
            # temporarily turn off logging to avoid cluttering the output
            logging.getLogger().setLevel(logging.ERROR)
            # run doctests
            (failed, tested) = doctest.testmod(m, name=module, optionflags=doctest.NORMALIZE_WHITESPACE)
            # restore normal logging level
            logging.getLogger().setLevel(log_level)
            if failed>0:
                failures += 1
                logging.error("  module '%s' FAILED %d tests out of %d." % (module, failed, tested))
                sys.stdout.write("Module '%s' FAILED %d tests out of %d.\n" % (module, failed, tested))
            else:
                if tested>0:
                    logging.debug("  OK - module '%s' passed all doctests." % module)
                    sys.stdout.write("Module '%s' OK, passed all doctests.\n" % module)
                else:
                    logging.warning("  module '%s' had no doctests." % module)
        finally:
            file.close()

    def run_homology_selftest(output=sys.stdout):
        ok = True
        # second, try known cases and inspect results
        for (g, n, ok) in [ (0,3, [1,0,0]),
                            (0,4, [1,2,0,0,0,0]),
                            (0,5, [1,5,6,0,0,0,0,0,0]),
                            (1,1, [1,0,0]),
                            (1,2, [1,0,0,0,0,0]),
                            (2,1, [1,0,1,0,0,0,0,0,0]),
                            ]:
            output.write("  Computation of M_{%d,%d} homology: " % (g,n))
            # compute homology of M_{g,n}
            timing.start("homology M%d,%d" % (g,n))
            hs = compute_homology(g,n)
            timing.stop("homology M%d,%d" % (g,n))
            # check result
            if hs == ok:
                output.write("OK (elapsed: %0.3fs)\n"
                             % timing.get("homology M%d,%d" % (g,n)))
            else:
                logging.error("Computation of M_{%d,%d} homology: FAILED, got %s expected %s"
                              % (g,n,hs,ok))
                output.write("FAILED, got %s expected %s\n" % (hs,ok))
                ok = False
        return ok

    ## run the self-test suite 3 times:
    ##

    ## 1. without any persistence stuff enabled, so we can test the
    ## real results of the algorithm and the performance
    sys.stdout.write("Checking homology algorithm (no checkpointing)\n")
    
    try:
        del runtime.options.checkpoint_dir
    except AttributeError:
        pass

    ok = run_homology_selftest(sys.stdout)
    if not ok:
        failures += 1

    ## 2. run with a temporary checkpoint directory, to test saving of state
    sys.stdout.write("Checking homology algorithm (checkpointing)\n")
        
    runtime.options.checkpoint_dir = tempfile.mkdtemp(prefix="mgn.selftest.")
    
    ok = run_homology_selftest(sys.stdout)
    if not ok:
        failures += 1

    ## 3. run with the same temporary directory, to test that
    ## persistence picks up results of rpevious runs correctly
    sys.stdout.write("Checking homology algorithm (restoring from checkpointed state)\n")
            
    ok = run_homology_selftest(sys.stdout)
    if not ok:
        failures += 1

    # remove anything in the temporary directory
    if failures == 0:
        for entry in os.listdir(runtime.options.checkpoint_dir):
            os.remove(os.path.join(runtime.options.checkpoint_dir, entry))
        os.rmdir(runtime.options.checkpoint_dir)
    else:
        sys.stdout.write("Persisted files left in directory '%s'"
                         % runtime.options.checkpoint_dir)

    # exit code >0 is number of failures
    sys.exit(failures)
        
        
# common code for invocations of `graphs`, `homology` and `valences`
if len(cmdline.args) < 2:
    parser.print_help()
    sys.exit(1)
try:
    g = int(cmdline.args[0])
    if g < 0:
        raise ValueError
except ValueError:
    sys.stderr.write("Invalid value '%s' for argument G: " \
                     "should be positive integer.\n" \
                     % (cmdline.args[0],))
    sys.exit(1)
try:
    n = positive_int(cmdline.args[1])
except ValueError, msg:
    sys.stderr.write("Invalid value '%s' for argument N: " \
                     "should be non-negative integer.\n" \
                     % (cmdline.args[1],))
    sys.exit(1)

# make g,n available to loaded modules
runtime.g = g
runtime.n = n

# ensure checkpoint path is defined and valid
if runtime.options.checkpoint_dir is None:
    runtime.options.checkpoint_dir = os.path.join(os.getcwd(), "M%d,%d.data" % (g,n))
    if not os.path.isdir(runtime.options.checkpoint_dir):
        if os.path.exists(runtime.options.checkpoint_dir):
            logging.error("Checkpoint path '%s' exists but is not a directory. Aborting.",
                          runtime.options.checkpoint_dir)
            sys.exit(1)
        else:
            os.mkdir(runtime.options.checkpoint_dir)
logging.info("Saving computation state to directory '%s'", runtime.options.checkpoint_dir)
if not runtime.options.restart:
    logging.warning("NOT restarting: will ignore any saved state in checkpoint directory '%s'")


# valences -- show vertex valences for given g,n
elif 'valences' == cmdline.action:
    logging.debug("Computing vertex valences occurring in g=%d,n=%d fatgraphs ...", g, n)
    vvs = compute_valences(g,n)
    for vv in vvs:
        outfile.write("%s\n" % str(vv))


# graphs -- create graphs from given g,n but do not compute homology
elif "graphs" == cmdline.action:
    logging.info("Will save graph list files into directory '%s'.",
                 runtime.options.checkpoint_dir)
    graphs, D = compute_graphs(g,n)
    logging.info("Graph family computation took %.3fs.",
                 timing.get("compute_graphs(%d,%d)" % (g,n)))
    

# homology -- compute homology ranks
elif 'homology' == cmdline.action:
    # compute graph complex and its homology ranks
    hs = compute_homology(g, n)
    logging.info("Homology computation took %.3fs.",
                 timing.get("compute_homology(%d,%d)" % (g,n)))

    # print results
    for (i, h) in enumerate(hs):
        outfile.write("h_%d(M_{%d,%d}) = %d\n" % (i, g, n, h))
    if cmdline.outfile is not None:
        logging.info("Results written to file '%s'" % cmdline.outfile)


# latex -- pretty-print graph lists
elif "latex" == cmdline.action:
    if runtime.options.checkpoint_dir is not None:
        dir = runtime.options.checkpoint_dir
    else:
        msg = ("Missing path to directory where graph list files are stored."
               " Set it using the `-s` option.")
        logging.error(msg)
        sys.stderr.write(msg + '\n')
        sys.exit(1)

    max_num_edges = 6*g + 3*n - 6
    min_num_edges = 2*g + n - 1
    max_num_vertices = 4*g + 2*n - 4

    ## read list files from checkpoint directory
    all_graphs = defaultdict(list)
    for graph in MgnGraphsIterator(g,n):
        all_graphs[graph.num_edges].append(graph)

    import output
    if runtime.options.outfile is not None:
        outfile = open(runtime.options.outfile, 'w')
    else:
        outfile = sys.stdout
    outfile = output.LaTeXOutput(outfile,
                                 g=g, n=n, version=__version__,
                                 checkpoint_dir=runtime.options.checkpoint_dir)

    pools = None
    for num_edges, graphs in all_graphs.iteritems():
        outfile.start_section(num_edges, max_num_vertices-(max_num_edges-num_edges))
        pools_prev = pools
        pools = [ NumberedFatgraphPool(G) for G in graphs ]
        # make list of partial sums for mapping matrix indices to graph numbers
        thresholds = [ 0 ]
        for p in (pools_prev or []):
            thresholds.append(len(p) + thresholds[-1])
        def matrix_index_to_G(i):
            if i == 0:
                return (0, 0)
            i0 = 0
            while thresholds[i0] < i:
                i0 += 1
            i0 -= 1
            return (i0, i-i0)
        def labelfn(D, i):
            i0, i1 = matrix_index_to_G(i)
            return ("G_{%d,%d}^{(%d)}" % (num_edges-1, i0, i))

        p = len(all_graphs[num_edges-1]) if (num_edges-1 in all_graphs) else 0
        q = len(all_graphs[num_edges]) if (num_edges in all_graphs) else 0
        r = num_edges - min_num_edges + 1
        matrix_file = os.path.join(dir, ("M%d,%d-D%d.sms" % (g,n,r)))
        if os.path.exists(matrix_file):
            d = SimpleMatrix(p, q)
            d.load(matrix_file)

        k0 = 0
        for j, G in enumerate(graphs):
            pool = pools[j]
            name = ("G_{%d,%d}" % (num_edges, j))
            outfile.start_graph(G, name)
            # print automorphisms
            Aut = list(G.automorphisms())
            if len(Aut) > 1:
                outfile.add_automorphisms(Aut)
            if G.is_oriented():
                # print markings
                if n > 1:
                    outfile.add_markings(name, pool)
                # print differential
                if d is not None and d.num_rows > 0 and d.num_columns > 0:
                    outfile.add_differential_start()
                    for k in xrange(len(pool)):
                        name_k = ("%s^{(%d)}" % (name, k))
                        outfile.add_differential(d, k0 + k, name_k, labelfn=labelfn)
                    outfile.add_differential_end()
                k0 += len(pool)
            outfile.end_graph()
        outfile.end_section()
    outfile.close()


else:
    sys.stderr.write("Unknown action `%s`, aborting.\n" % cmdline.action)
    sys.exit(1)


# try to print profiling information, but ignore failures
# (profiling might not have been requested, in the first place...)
try:
    pf.stop()
    logging.debug("Stopped call profiling, now dumping stats.")
    pf.close()
    import hotshot.stats
    stats = hotshot.stats.load(__name__ + '.pf')
    stats.strip_dirs()
    stats.sort_stats('cumulative', 'calls')
    stats.print_stats(50)
except:
    pass

# print CPU time usage
cputime_s = resource.getrusage(resource.RUSAGE_SELF)[0]
seconds = cputime_s % 60
minutes = int(cputime_s / 60) % 60
hours = int(cputime_s / 3600) % 24
days = int(cputime_s / 86400)
if days > 0:
    elapsed = "%d days, %d hours, %d minutes and %2.3f seconds" % (days, hours, minutes, seconds)
elif hours > 0:
    elapsed = "%d hours, %d minutes and %2.3f seconds" % (hours, minutes, seconds)
elif minutes > 0:
    elapsed = "%d minutes and %2.3f seconds" % (minutes, seconds)
else:
    elapsed = "%2.3f seconds" % seconds
logging.info("Total CPU time used: " + elapsed)

# That's all folks!
logging.debug("Done: %s" % str.join(" ", sys.argv))
