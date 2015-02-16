import copy
import sys
import numpy as np


class QuartetsInfo:
    def __init__(self,c,infile,outfile):
        self.infile = infile
        self.outfile = outfile
        self.quartet_dict_ = None
        self.c = c
        #self.quartet_dict = {}
        #with open(quartetsfile) as f:
            #for line in f:
                #(key,val) = line.split()
                #self.quartet_dict[key] = int(val)

    def quartet_dict(self):
        if self.quartet_dict_:
            return self.quartet_dict_
        d = {}
        with open(self.infile) as f:
            for line in f:
                (key,val) = line.split()
                d[key] = int(val)
        #print d
        self.quartet_dict_ = d
        return d

#s is a quartet string like ((a,b),(c,d));
    def get_freqs(self,Q):
        s = self.quartet_dict()
        #q = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));']
        #u = [0,0,0]
        #for i in range(3):
        #    if q[i] in s:
        #        u[i] = s[q[i]]
        L = copy.copy(Q)
        L = L.replace('(','')
        L = L.replace(')','')
        L = L.replace(';','')
        L = L.split(',')
        u = [0,0,0]
        u[0] = s[Q]
        q = [Q, None, None]
        #q0 = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[1]+'),('+L[3]+','+L[2]+'));','(('+L[1]+','+L[0]+'),('+L[2]+','+L[3]+'));','(('+L[1]+','+L[0]+'),('+L[3]+','+L[2]+'));', '(('+L[2]+','+L[3]+'),('+L[0]+','+L[1]+'));','(('+L[2]+','+L[3]+'),('+L[1]+','+L[0]+'));','(('+L[3]+','+L[2]+'),('+L[0]+','+L[1]+'));','(('+L[3]+','+L[2]+'),('+L[1]+','+L[0]+'));']
        q1 = ['(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[3]+','+L[1]+'));','(('+L[2]+','+L[0]+'),('+L[1]+','+L[3]+'));','(('+L[2]+','+L[0]+'),('+L[3]+','+L[1]+'));','(('+L[1]+','+L[3]+'),('+L[0]+','+L[2]+'));','(('+L[1]+','+L[3]+'),('+L[2]+','+L[0]+'));','(('+L[3]+','+L[1]+'),('+L[0]+','+L[2]+'));','(('+L[3]+','+L[1]+'),('+L[2]+','+L[0]+'));']
        q2 = ['(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));', '(('+L[0]+','+L[3]+'),('+L[2]+','+L[1]+'));','(('+L[3]+','+L[0]+'),('+L[1]+','+L[2]+'));', '(('+L[3]+','+L[0]+'),('+L[2]+','+L[1]+'));','(('+L[1]+','+L[2]+'),('+L[0]+','+L[3]+'));','(('+L[1]+','+L[2]+'),('+L[3]+','+L[0]+'));','(('+L[2]+','+L[1]+'),('+L[0]+','+L[3]+'));','(('+L[2]+','+L[1]+'),('+L[3]+','+L[0]+'));']
        for i in range(8):
            if q1[i] in s:
                u[1] = s[q1[i]]
                q[1] = q1[i]
        for i in range(8):
            if q2[i] in s:
                u[2] = s[q2[i]]
                q[2] = q2[i]
        return u, q
        
    def new_dict(self):
        s = self.quartet_dict()
        d = {}
        for q in s:
            u, q = self.get_freqs(q)
            if u[0] > u[1] and u[0] > u[2]:
                d[q[0]] = ((q[0], u[0]), (q[1], u[1]), (q[2], u[2]))
        m = max(d.values())
        return d
        #print e

    def write_file(self):
        d = self.new_dict()
        H = open(self.outfile,'w')
        for q in np.random.choice(np.array(d.keys()), len(d) * self.c):
            if d[q][0][1]:
                l = str(d[q][0][0]) + ' ' + str(d[q][0][1]) + '\n'
                H.write(l)
            if d[q][1][1]:
                l = str(d[q][1][0]) + ' ' + str(d[q][1][1]) + '\n'
                H.write(l)
            if d[q][2][1]:
                l = str(d[q][2][0]) + ' ' + str(d[q][2][1]) + '\n'
                H.write(l)
        H.close()


