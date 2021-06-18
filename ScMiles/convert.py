import struct
import sys

fbin = 'input.vel'
vel = []
binfile = open(fbin,'rb')
context = binfile.read(4)
natom = struct.unpack("i",context)[0]
vel = []
for i in range(natom*3):
    context = binfile.read(8)
    vel.append(-struct.unpack("d",context)[0])
binfile.close()

fout = 'reversed.vel'
fp = open(fout,"wb")
fp.write(struct.pack('i',natom))
for x in vel:
    fp.write(struct.pack('d',x))
