# Echo client program
import socket
import sys
import os
import math
import struct
from common import *


coded_files_dir = 'coded_files/'
helper_files_dir = 'helper_files/'

repair_prgm = './repair_pmc' if os.name != 'nt' else 'repair_pmc.exe'
collector_prgm = './collector_pmc' if os.name != 'nt' else 'collector_pmc.exe'

for dir_name in [coded_files_dir, helper_files_dir]:
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)


def sendDCRequest(host, port, filePrefix):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((host, port))
    print 'Connected to host: ', host, ', port: ', port
    sock = BlockingSocket(s)

    #print 'DC Request filePrefix: ', filePrefix
    sock.send("DC")
    filePrefixSize = len(filePrefix)
    #print 'Sending Prefix Size: ', filePrefixSize
    sendInt(sock, filePrefixSize)
    #print 'Sending filePrefix: ', filePrefix
    sock.send(filePrefix)
    numFiles = readInt(sock)
    print 'Number of files to receive: ', numFiles
    for i in xrange(numFiles):
        recvFile(sock, coded_files_dir)



def sendNRRequest(host, port, fid, filePrefix):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((host, port))
    print 'Connected to host: ', host, ', port: ', port
    sock = BlockingSocket(s)
    sock.send("NR")

    filePrefixSize = len(filePrefix)
    #print 'Sending Prefix Size: ', filePrefixSize
    sendInt(sock, filePrefixSize)
    #print 'Sending filePrefix: ', filePrefix
    sock.send(filePrefix)
    # Send the failed node id
    sendInt(sock, fid)

    # Recive actual helper-chunk
    recvFile(sock, helper_files_dir)

    # Recieve metadata (TODO: split this request)
    recvFile(sock, coded_files_dir)
    


    

def dataCollect(all_hosts, dc_node_ids, filename):
    i = 0
    for nid in dc_node_ids:
        i += 1
        host, port = all_hosts[nid]
        print "                       > Collecting from node #%d / %d <" % (i , len(dc_node_ids))
        sendDCRequest(host, port, filePrefix=filename)

    enc_metafile = coded_files_dir + filename + "_metadata"
    os.system("%s %s %s" % (collector_prgm, enc_metafile, ' '.join([str(ni) for ni in dc_node_ids]))) 

def nodeRepair(fid, all_hosts, helper_node_ids, filename):
    i = 0
    for nid in helper_node_ids:
        i += 1
        host, port = all_hosts[nid]
        print "                       > Recieving from helper-node #%d / %d <" % (i , len(helper_node_ids))
        sendNRRequest(host, port, fid, filePrefix=filename)

    enc_metafile = coded_files_dir + filename + "_metadata"
    os.system("%s fix %d %s %s" % (repair_prgm, fid, enc_metafile, ' '.join([str(ni) for ni in helper_node_ids]))) 


def sendAll(all_hosts, filename):
    i = 0
    for hid in all_hosts:
        i += 1
        try:
            host, port = all_hosts[hid]
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.connect((host, port))
            print 'Connected to host: ', host, ', port: ', port
            sock = BlockingSocket(s)
            sock.send("FS")

            print "                       > Sending to node #%d / %d <" % (i , len(all_hosts))
        
            #

            nodeName = "node" + str(hid)
            chunk_fnames = [fname for fname in os.listdir(coded_files_dir) if nodeName in fname]

            # num files being sent
            sendInt(sock, len(chunk_fnames) + 1)

            # first send encoder metadata file
            sendFile(sock, coded_files_dir, filename + "_metadata")

            # then all relevant chunks
            for fname in chunk_fnames:
                sendFile(sock, coded_files_dir, fname)
        except:
            print "*********************   ERROR sending to host #%d @ %s:%d   *********************" % (hid, host, port)


if __name__=='__main__':
    if len(sys.argv) < 4:
        print "usage (data collect) : DC filename node_info_file [participating nodes]"
        print "usage (node repair)  : NR failed_id filename node_info_file [helper nodes]"
        print "usage (file store)   : FS filename node_info_file"
        print "\t nodes : a list of node ids"
        print "\t node_info : the metadata file containing list of hosts"
        print ""
        print "DC/NR : performs data-collection/node-repair."
        print "FS    : sends encoded files in the coded_files dir to the apprioriate node for storage"
        sys.exit(-1)

    reqType = sys.argv[1]
    if reqType not in ["DC", "NR", "FS"]:
        print 'invalid request type ', reqType
        sys.exit(-1)

    if reqType == "DC":
        filename = sys.argv[2]
        meta_fname = sys.argv[3]
        nodes = [int(node_id) for node_id in sys.argv[4:]]
    elif reqType == "NR":
        fid = int(sys.argv[2])
        filename = sys.argv[3]
        meta_fname = sys.argv[4]
        nodes = [int(node_id) for node_id in sys.argv[5:]]
    elif reqType == "FS":
        filename = sys.argv[2]
        meta_fname = sys.argv[3]

    f = open(meta_fname)
    all_hosts = {int(i) : (host, int(port)) for i, host, port in [ln.split(" ") for ln in f.read().splitlines() if len(ln) > 0]}

    f.close()

    if reqType == "DC":
        dataCollect(all_hosts, nodes, filename)
    if reqType == "NR":
        nodeRepair(fid, all_hosts, nodes, filename)
    if reqType == "FS":
        sendAll(all_hosts, filename)
