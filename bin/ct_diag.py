#!/usr/bin/env python3
# Diagnostic: run the REAL infer_hierarchy entry point with flushed, non-redirected stdout so we can
# watch live how many rounds complete and how potential_trees grows per added mutation. run_clonetracer.py
# redirects sys.stdout to an unflushed file, hiding all of this. Not vendored - debug only.
import sys, os, time, torch
from datetime import datetime

def log(m):
    print(f"{datetime.now().strftime('%H:%M:%S')} >>> {m}", flush=True)

json_in  = sys.argv[1] if len(sys.argv) > 1 else 'validate_ct/HD_BM_3_sub8.json'
num_iter = int(sys.argv[2]) if len(sys.argv) > 2 else 60
out_dir  = sys.argv[3] if len(sys.argv) > 3 else 'validate_ct'

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
torch.set_default_tensor_type(torch.DoubleTensor)
log(f"torch {torch.__version__} threads={torch.get_num_threads()} OMP={os.environ.get('OMP_NUM_THREADS')}")

from helper_functions import create_tree_class
log("imported helper_functions")

t = create_tree_class(json_in, 'diag_ih', False, False, False)
log(f"tree class: nmut={len(t.names)} M.shape={tuple(t.M.shape)}")

# Wrap select_tree so each round prints its candidate-tree count + timing live.
orig_select = t.select_tree
round_state = {'n': 0}
def traced_select(num_iter_, init_, *a, **k):
    round_state['n'] += 1
    s = time.time()
    log(f"ROUND {round_state['n']}: muts={list(t.muts)} candidates={len(t.potential_trees)} (num_iter={num_iter_})")
    r = orig_select(num_iter_, init_, *a, **k)
    log(f"ROUND {round_state['n']} done in {time.time()-s:.1f}s -> selected {getattr(t,'tree_indices',None)}")
    return r
t.select_tree = traced_select

log(f"calling infer_hierarchy(num_iter={num_iter}, init={num_iter-100 if num_iter>100 else 20}) ...")
init = num_iter - 100 if num_iter > 100 else 20
s = time.time()
t.infer_hierarchy(num_iter, init, out_dir)
log(f"infer_hierarchy DONE in {time.time()-s:.1f}s")
log(f"outputs: {[f for f in os.listdir(out_dir) if 'diag_ih' in f]}")
