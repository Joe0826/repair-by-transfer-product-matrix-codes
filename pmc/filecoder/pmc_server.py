# Echo server program
import socket
import os
import sys
import struct
from common import *

coded_files_dir = 'coded_files/'
helper_files_dir = 'helper_files/'

repair_prgm = './repair_pmc' if os.name != 'nt' else 'repair_pmc.exe'

for dir_name in [coded_files_dir, helper_files_dir]:
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

if len(sys.argv) < 3:
    print "usage: port sid [host ip]"
    print ""
    print "    sid : the node id of this server"
    print "    ip  : the local ip of this server (will attempt to resolve if absent)"
    exit(-1)
PORT = int(sys.argv[1])
sid = int(sys.argv[2])    # the server node id

HOST = socket.gethostbyname(socket.gethostname()) if len(sys.argv) < 4 else sys.argv[3]
print 'Host: ', HOST
                                         
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((HOST, PORT))
s.listen(1)

#Helper functions
def sendInt(conn, num):
    conn.send(struct.pack('I', num))

def readInt(conn):
    return int(struct.unpack('I', conn.recv(4))[0])

def handleDCRequest(conn):
    prefixLen = readInt(conn)
    #print 'prefixLen: ', prefixLen
    prefix = conn.recv(prefixLen)
    #print 'prefix: ', prefix
    filenames = [fname for fname in os.listdir(coded_files_dir) if fname.startswith(prefix)]
    print 'Number of files to send: ', len(filenames)
    sendInt(conn, len(filenames))
    for fname in filenames:
        sendFile(conn, coded_files_dir, fname)
    conn.close()


def handleNRRequest(conn):
    prefixLen = readInt(conn)
    print 'prefixLen: ', prefixLen
    prefix = conn.recv(prefixLen)
    print 'prefix: ', prefix
    fid = readInt(conn)

    enc_metafile = coded_files_dir + prefix + "_metadata"
    os.system("%s help %d %s %d" % (repair_prgm, fid, enc_metafile, sid))

    # Send the helper file
    fname = prefix + '_helper_chunk_from' + str(sid) + '_to' + str(fid)
    sendFile(conn, helper_files_dir, fname)

    # Send the original encoding metadata file (TODO: separate this request)
    sendFile(conn, coded_files_dir, prefix + "_metadata")

    conn.close()


def handleFSRequest(conn):
    numFiles = readInt(conn)
    for _ in xrange(numFiles):
        recvFile(conn, coded_files_dir)

#Handle requests
while True:
    connection, addr = s.accept()

    conn = BlockingSocket(connection)
    print 'Connected by: ', addr
    request = conn.recv(2)
    print 'Received request: ', request
    if request == "DC":
        handleDCRequest(conn)
    elif request == "NR":
        handleNRRequest(conn)
    elif request == "FS":
        handleFSRequest(conn)
    else:
        conn.send('0')
    conn.close()
