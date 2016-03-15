import sys
import os
import os.path
import math
from subprocess import call, Popen, PIPE
import random
import csv

VERBOSE = False

fname = "fakef"
code_dir = "coded_files"
mfname = "%s/%s_metadata" % (code_dir, fname)
sizeof_long = 8 #sizeof(long) in the C compiler...


def write_metadata(n, k, d, w, a, B, msg_size, tech, rbt):
    mf = open(mfname, 'w')
    rbt_val = 0
    if rbt == 'SYS':
        rbt_val = 1
    elif rbt == 'CYC':
        rbt_val = 2

    for line in [fname, msg_size, msg_size, n, k, d, a, B, w,
            ('1' if tech == 'MBR' else '0'),
            str(rbt_val)
            ]:
        mf.write(str(line) + "\n")
    mf.close()

def gen_fake_files(msg_size, a, B, helper_ids):
    print "generating fake coded files..."
    chunk_size = msg_size / B

    for h in helper_ids:
        for s in range(a):
            chunk_fname = "%s/%s_node%d_sym%d" % (code_dir, fname, h, s)

            # if file exists and has size >= chunk_size, skip creation
            if os.path.isfile(chunk_fname) and os.path.getsize(chunk_fname) >= chunk_size:
                continue

            cmd = "dd if=/dev/urandom of=%s bs=1024 count=%d 2> /dev/null" % (chunk_fname, int(math.ceil(chunk_size / 1024.0)))
            print cmd
            call(cmd, shell=True, stdout=None, stderr=None)

def run_repair(n, k, d, w, helper_ids, failed_id):
    htime = ftime = 0.0

    print "running repair [HELP]..."
    cmd1 = "./repair_pmc help %d %s %s" % (failed_id, mfname, " ".join(map(str, helper_ids)))
    print cmd1
    stdout = Popen(cmd1, shell=True, stdout=PIPE, stderr=sys.stderr).stdout

    for line in stdout:
        if VERBOSE:
            sys.stdout.write(line) # echo output
        if 'Coding Time' in line:
            htime += float(line.split(':')[1].strip())
    stdout.close()

    print "running repair [FIX]..."
    cmd2 = "./repair_pmc fix %d %s %s" % (failed_id, mfname, " ".join(map(str, helper_ids)))
    print cmd2
    stdout = Popen(cmd2, shell=True, stdout=PIPE, stderr=sys.stderr).stdout

    for line in stdout:
        if VERBOSE:
            sys.stdout.write(line) # echo output
        if 'Coding Time' in line:
            ftime += float(line.split(':')[1].strip())

    return htime, ftime


def calc_proper_fsize(targ_size):
    #Determine size of file : needs to be a multiple for B*sizeof(long).
    proper_fsize = targ_size
    residue = targ_size % (B*sizeof_long);
    if residue != 0:
        proper_fsize += B*sizeof_long - residue;
    return proper_fsize

def calc_params(n, k, d, tech):
    if tech == "MSR":
        a   = d-k+1;
        B   = k*a;
    else:
        a   = d;
        B   = k*d - k*(k-1)/2;
    return a, B

def choose_fail(k):
    if rbt == None:
        return random.randrange(n)
    if rbt == 'SYS':
        return random.randrange(k) # choose a random sys node to fail

def choose_helper_list(node_id, rbt):
    helper_list = []
    candidates = list(xrange(n))
    candidates.remove(node_id)
    if rbt == None or rbt == 'SYS':
        # choose d random nodes
        return random.sample(candidates, d)


if __name__ == "__main__":
    if len(sys.argv) < 7:
        print "usage:  python time_repair.py [-v] n k d w {MSR/MBR} node_blocksize [rbt]"
        print ""
        print "\tnode_blocksize : the amount of data each node stores (MB)"
        exit()

    if (sys.argv[1] == '-v'):
        VERBOSE = True
        sys.argv = sys.argv[1:]

    n, k, d, w, bs = map(int, sys.argv[1:5] + [sys.argv[6]])
    tech = sys.argv[5]
    rbt = None
    if len(sys.argv) == 8:
        rbt = sys.argv[7]

    bs *= 1024 * 1024 # MB --> bytes

    a, B = calc_params(n, k, d, tech)

    target_fsize = bs / a * B # file-size s.t. each node stores 'bs' bytes of data

    msg_size = calc_proper_fsize(target_fsize)
    write_metadata(n, k, d, w, a, B, msg_size, tech, rbt)


    failed_id = choose_fail(k)

    helper_ids = choose_helper_list(failed_id, rbt)

    gen_fake_files(msg_size, a, B, helper_ids)
    htime, ftime = run_repair(n, k, d, w, helper_ids, failed_id)
    total_time = htime + ftime
    print ""
    print "TOTAL CODING TIME [s] : ", (htime + ftime)

    # Record the coding time
    csv_name = 'n%d_k%d_d%d_size%d_%s_%s_repair.csv' % (n, k, d, bs / 1024 / 1024, tech, rbt)
    with open(csv_name, 'a') as csvfile:
        csvfile.write(str(total_time) + ',')


