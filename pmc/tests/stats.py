import os, os.path as osp
import sys
import time
import random


n = 50
k = 20
d = 40
w = 8
MXR     = 'MSR'
fdat    = 'f1M.dat'

Ntrials = 10

jdir    = '../../jerasure/Examples'
pmcdir  = '../filecoder_example'


def time_encoding():
    """
    Returns the time in seconds it takes to do encoding
    """
    t_tot = 0
    curdir = os.getcwd()
    os.chdir(pmcdir)
    for i in xrange(Ntrials):
        os.system('rm -r coded_files')
        tstart = time.time()
        os.system('./encoder_pmc %d %d %d %d %s %s'% (n, k, d, w, MXR, fdat) )
        t_tot += time.time() - tstart
    os.chdir(curdir)
    return t_tot/Ntrials


def time_dc(use_systematic=False):
    """
    Returns the time in seconds it takes to do Data collection.
    If use_systematic == True : uses just the first k-nodes for data-collection.
    """
    # do encoding:
    curdir = os.getcwd()
    os.chdir(pmcdir)
    

    os.system('rm -r coded_files')
    os.system('./encoder_pmc %d %d %d %d %s %s'% (n, k, d, w, MXR, fdat) )
    os.system('mv %s %s_orig'%(fdat, fdat))
    
    t_tot = 0
    for i in xrange(Ntrials):

        klist = xrange(k) if use_systematic else random.sample(xrange(n), k)

        tstart = time.time()
        os.system("./collector_pmc coded_files/%s_metadata %s" % (fdat, ' '.join([str(x) for x in klist]))  )
        t_tot += time.time() - tstart

        os.system('rm %s'%fdat)

    os.system('mv %s %s'%(fdat+'_orig', fdat))

    os.chdir(curdir)
    return t_tot/Ntrials

def time_repair(systematic_helpers=False):
    """
    Returns the time in seconds it takes to do Data collection.
    If use_systematic == True : uses just the first k-nodes for data-collection.
    """
    # do encoding:
    curdir = os.getcwd()
    os.chdir(pmcdir)
    

    os.system('rm -r coded_files')
    os.system('rm -r helper_files')
    os.system('./encoder_pmc %d %d %d %d %s %s'% (n, k, d, w, MXR, fdat) )
    
    t_tot = 0
    for i in xrange(Ntrials):

        fid = random.randint(0, n-1)
        helpers = xrange(d) if systematic_helpers else random.sample([i for i in range(n) if i != fid], d)
        argtuple = (fid, fdat, ' '.join([str(x) for x in helpers]))

        tstart = time.time()
        os.system("./repair_pmc help %d coded_files/%s_metadata %s && ./repair_pmc fix %d coded_files/%s_metadata %s" % (argtuple + argtuple))
        t_tot += time.time() - tstart

        os.system('rm -r helper_files')

    os.chdir(curdir)
    return t_tot/Ntrials



def time_jencoding():
    """
    Returns the time in seconds it takes to do encoding
    """
    curdir = os.getcwd()
    os.chdir(jdir)
    t_tot = 0
    for i in xrange(Ntrials):
        os.system('rm -r Coding')
        tstart = time.time()
        os.system('./encoder %s %d %d reed_sol_van %d 0 0'% (fdat, k, n-k, w) )
        t_tot += time.time() - tstart

    os.chdir(curdir)
    return t_tot/Ntrials


def get_pad_str(d,n):
    return '0'*(n-len(str(d))) + str(d)

def time_jdc(recode=False, use_systematic=False):
    """
    Returns the time in seconds it takes to do Data collection.
    If use_systematic == True : uses just the first k-nodes for data-collection.
    """
    curdir = os.getcwd()
    os.chdir(jdir)


    fdatp, fdatext = os.path.splitext(fdat)
    print fdatp, fdatext

    # do encoding:
    t_tot = 0
    for i in xrange(Ntrials):

        os.system('rm -r Coding')
        os.system('./encoder %s %d %d reed_sol_van %d 0 0'% (fdat, k, n-k, w) )
        
        filenames = os.listdir('Coding')
        
        tn = len(str(n))
        klist = xrange(k) if use_systematic else random.sample(xrange(n), k)
        for ki in xrange(0,n):
            if ki not in klist: ## rm all the non-DC node files
                ki_fname = [fname for fname in filenames if fname.startswith('%s_%s'%(fdatp, 'k' if ki < k else 'm'))  and fname.endswith(get_pad_str( ki+1 if ki < k else ki - k+1,tn)+fdatext)][0]
                os.system('rm Coding/%s'%ki_fname)

        tstart = time.time()
        if recode:
            os.system("./decoder %s" % fdat)
        else:
            os.system("./decoder_norecode %s" % fdat)

        t_tot += time.time() - tstart

    os.chdir(curdir)
    return t_tot/Ntrials



if __name__ == '__main__':
    os.system('cp %s %s/'%(fdat, jdir))
    os.system('cp %s %s/'%(fdat, pmcdir))

    t_dc  = time_dc()
    t_enc  =  time_encoding()
    t_repair = time_repair(systematic_helpers = False)

    t_jdc = time_jdc()
    t_jenc = time_jencoding()
    t_jrepair = time_jdc(recode=True)

    print "Data-collection time (jerasure) [s]: ", t_jdc


    print "Encoding time [s]: ", t_enc
    print "Data-collection time: ", t_dc
    print "Repair time [s]: ", t_repair 



    print "J Encoding time [s]: ", t_jenc
    print "J Data-collection time: ", t_jdc
    print "J Repair time [s]: ", t_jrepair 


    os.system('rm %s/%s'%(jdir, fdat))
    os.system('rm %s/%s'%(pmcdir, fdat))
