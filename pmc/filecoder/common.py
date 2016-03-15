import socket
import sys
import os
import struct

class BlockingSocket:
    def __init__(self, baseSocket, buffsize = 8192):
        self.buffsize = buffsize
        self.s = baseSocket
        
    # receives and returns SIZE bytes (blocks)
    def recv(self, size):
        currBytes = ''
        while len(currBytes) < size:
                toread = min(self.buffsize, size - len(currBytes))
                currBytes += self.s.recv(toread)
        return currBytes

    def send(self, data):
        self.s.send(data)

    def close(self):
        self.s.close()

#Helper functions
def sendInt(socket, num):
    socket.send(struct.pack('I', num))

def readInt(socket):
    return int(struct.unpack('I', socket.recv(4))[0])

#Sends filename and then file
def sendFile(conn, folder, fname,blocksize=8192):
    #print 'Sending filename size: ', len(fname)
    sendInt(conn, len(fname))
    conn.send(fname)
    print 'Sending file: ', fname
    #print 'Size of file: ', os.path.getsize(folder + fname)
    sendInt(conn, os.path.getsize(folder + fname))
    f = open(folder + fname, 'rb')
    totalBytes = os.path.getsize(folder + fname)
    currBytes = 0
    while currBytes < totalBytes:
        data = f.read(blocksize)
        conn.send(data)
        currBytes += len(data)
    f.close()

def recvFile(sock, dest_dir):
    filenamesize = readInt(sock)
    #print 'Receiving filename size: ', filenamesize
    filename = sock.recv(filenamesize)
    print 'Receiving: ', filename
    filesize = readInt(sock)
    print 'Size: ', filesize
    data = sock.recv(filesize)
    print 'Saving file in dir: ', dest_dir
    writeToFile(dest_dir + filename, data)

def writeToFile(dest, data):
    fname = open(dest, 'wb')
    fname.write(data)
    fname.close()