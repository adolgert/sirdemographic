import subprocess
import logging
import numpy as np
import itertools as it
from default_parser import DefaultArgumentParser

logger=logging.getLogger(__file__)


def matrix_run(res, run_cnt):
    seed=0
    total_individual=100000
    for iri in it.product(range(1,res), range(1,res)):
        ir=np.array(total_individual*np.array(iri)/res, dtype=np.int)
        sir=[total_individual, ir[0], ir[1]]
        fname="arr-{0}-{1}-{2}.h5".format(sir[0], sir[1], sir[2])
        to_run=['./sirexp', '-j','4', '--runcnt', str(run_cnt),
            '-s', str(sum(sir)), '-i', str(sir[1]), '-r', str(sir[2]),
            '--seed', str(seed), '--endtime', str(0.2), '--loglevel', 'warning',
            '--beta1', str(0), '--datafile', fname]
        logger.debug(to_run)
        ret=subprocess.Popen(to_run)
        result=ret.communicate()[0] # join the process
        seed+=1


if __name__=='__main__':
    parser=DefaultArgumentParser(description="Run sirexp many times")
    parser.add_function("exp", "Explore initial conditions")
    parser.add_argument("--res", dest="resolution", type=int, action="store",
        default=10, help="How finely to subdivide the SIR space")
    parser.add_argument("--runcnt", dest="runcnt", type=int, action="store",
        default=10, help="How many times to run each simulation")

    args=parser.parse_args()
    if args.exp:
        matrix_run(args.resolution, args.runcnt)
