import sys
import os
import os.path
import math
from subprocess import call, Popen, PIPE
import random

VERBOSE = False

fname = "fakef.txt"
fname_content, fname_suffix = fname.split('.')
sizeof_long = 8 #sizeof(long) in the C compiler...
digit = 1


def gen_fake_files(n, k, w, msg_size, helper_ids, j_tech='reed_sol_van'):
    print "generating fake coded files..."
    file_size = msg_size

    # if file exists and has size >= file_size, skip creation
    if not (os.path.isfile(fname) and os.path.getsize(fname) >= file_size):
        cmd = "dd if=/dev/urandom of=%s bs=1024 count=%d 2> /dev/null" % (fname, int(math.ceil(file_size / 1024.0)))
        print cmd
        call(cmd, shell=True, stdout=None, stderr=None)

    # Encode
    m = n - k
    cmd_encode = "./encoder %s %d %d %s %d 0 0" % (fname, k, m, j_tech, w)
    print 'Encoding...'
    call(cmd_encode, shell=True, stdout=PIPE, stderr=sys.stderr)

    # Remove the non-helping files
    for node_id in xrange(n):
        if node_id not in helper_ids:
            node_idx = get_node_name_j(node_id, k)
            node_file_name = '%s_%s.%s' % (fname_content, node_idx, fname_suffix)
            cmd = "rm Coding/%s" % node_file_name
            print cmd
            call(cmd, shell=True, stdout=None, stderr=None)

def get_node_name_j(node_id, k):
    if digit == 2:
        if node_id < k:
            return 'k%02d' % (node_id + 1)
        else:
            return 'm%02d' % (node_id - k + 1)
    if digit == 1:
        if node_id < k:
            return 'k%d' % (node_id + 1)
        else:
            return 'm%d' % (node_id - k + 1)

def run_repair(n, k, w, helper_ids, failed_id, j_tech='reed_sol_van'):
    time = 0.0
    m = n - k
    print "running repair..."
    cmd = './decoder_norecode %s' % fname
    print cmd
    stdout = Popen(cmd, shell=True, stdout=PIPE, stderr=sys.stderr).stdout

    for line in stdout:
        if VERBOSE:
            sys.stdout.write(line) # echo output
        if 'Coding Time' in line:
            time += float(line.split(':')[1].strip())
    stdout.close()

    '''
    print "encoding the falied node..."
    cmd_encode = "./encoder %s %d %d %s %d 0 0" % (fname, k, m, j_tech, w)
    stdout = Popen(cmd_encode, shell=True, stdout=PIPE, stderr=sys.stderr).stdout

    for line in stdout:
        if VERBOSE:
            sys.stdout.write(line) # echo output
        if 'Coding Time' in line:
            time += float(line.split(':')[1].strip())
    stdout.close()
    '''
    return time


def choose_fail(k):
    return random.randrange(k) # choose a random sys node to fail
    #if rbt == None:
        #return random.randrange(n)
    #if rbt == 'SYS':
        #return random.randrange(k) # choose a random sys node to fail

def choose_dc_list(node_id, rbt):
    helper_list = []
    candidates = list(xrange(n))
    candidates.remove(node_id)
    if rbt == None or rbt == 'RAND':
        # choose k random nodes
        return random.sample(candidates, k)
    elif rbt == 'SYS':
        # first k availible nodes
        return candidates[:k]
    elif rbt == 'PARITY':
        # last k availible nodes (parities)
        return candidates[-k:]


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "usage:  python time_repair_j.py [-v] n k w node_blocksize [SYS/RAND/PAIRITY]"
        print ""
        print "\tnode_blocksize : the amount of data each node stores (MB)"
        print "\t[SYS/RAND/PAIRITY] : whom to data-collect from"
        print "(always repair a systematic node)"
        exit()

    if (sys.argv[1] == '-v'):
        VERBOSE = True
        sys.argv = sys.argv[1:]

    global rbt
    rbt = 'RAND'
    if len(sys.argv) == 6:
        rbt = sys.argv[5]

    n, k, w, bs = map(int, sys.argv[1:5])
    global digit
    if n >= 10:
        digit = 2

    bs *= 1024 * 1024 # MB --> bytes

    msg_size = bs * k # file-size s.t. each node stores 'bs' bytes of data

    failed_id = choose_fail(k)

    helper_ids = choose_dc_list(failed_id, rbt)

    gen_fake_files(n, k, w, msg_size, helper_ids)

    time = run_repair(n, k, w, helper_ids, failed_id)

    print ""
    print "TOTAL CODING TIME [s] : ", (time)

    csv_name = 'n%d_k%d_%dMB_repair_j.csv' % (n, k, bs / 1024 / 1024)
    with open(csv_name, 'a') as csvfile:
        csvfile.write(str(time) + ',')
