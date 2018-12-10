from ctypes import *
import numpy as np

'''Python interface to gslib2
http://www.fortran90.org/src/best-practices.html#interfacing-with-python

from ctypes import CDLL, POINTER, c_int, c_double
from numpy import empty

fortran = CDLL('./libmyfortranroutines.so')

mesh = empty(N, dtype="double")
fortran.mesh_exp(c_double(r_min), c_double(r_max), c_double(a), c_int(N),
                 mesh.ctypes.data_as(POINTER(c_double)))
                 
https://docs.scipy.org/doc/numpy/reference/routines.ctypeslib.html
                 
numpy.ctypeslib.ndpointer(
'''

FLOAT_ARRAY_ARGTYPE = np.ctypeslib.ndpointer(dtype=np.float32,flags='F')
ROTMATRIX_ARRAY_ARGTYPE = np.ctypeslib.ndpointer(dtype=np.float32,ndim=2,flags='F_CONTIGUOUS')


__gslib2 = CDLL('../fortran/gslib2.so')

#sub array_example(x(:))
__gslib2.array_example.argtypes = (
        FLOAT_ARRAY_ARGTYPE, #x
        c_int,
    )

def array_example(x):
    __gslib2.array_example(x,len(x))

#function setup_superblock(x,y,z,vr,sec,GRID,rotmat,radsqd,MAXSBX,MAXSBY,MAXSBZ)
__gslib2.setup_superblock.argtypes = [
        FLOAT_ARRAY_ARGTYPE, #x
        FLOAT_ARRAY_ARGTYPE, #y
        FLOAT_ARRAY_ARGTYPE, #z
        FLOAT_ARRAY_ARGTYPE, #vr
        FLOAT_ARRAY_ARGTYPE, #sec
        c_int,   #n
        c_int,   #nx
        c_float, #xmn
        c_float, #xsiz
        c_int,   #ny
        c_float, #ymn
        c_float, #ysiz
        c_int,   #nz
        c_float, #zmn
        c_float, #zsiz
        ROTMATRIX_ARRAY_ARGTYPE, #rotmat
        c_float, #radsqd
        c_int, #MAXSBX
        c_int, #MAXSBY
        c_int #MAXSBZ
    ]
__gslib2.setup_superblock.restype = c_void_p

__gslib2.super_block_print.argtypes = [
        c_void_p, #
    ]
__gslib2.super_block_print.restype = None

__gslib2.search_super_block.argtypes = [
        c_void_p, #superblock structure (pointer)
        c_float, #xloc
        c_float, #yloc
        c_float, #zloc
        c_float, #radsqd
        ROTMATRIX_ARRAY_ARGTYPE, #rotmat
        c_int, #ndmax
        c_int, #noct
        FLOAT_ARRAY_ARGTYPE, #x
        FLOAT_ARRAY_ARGTYPE, #y
        FLOAT_ARRAY_ARGTYPE, #z
        c_int,   #nd
        c_int,   #nclose
        FLOAT_ARRAY_ARGTYPE, #close
        FLOAT_ARRAY_ARGTYPE, #infoct
    ]
__gslib2.search_super_block.restype = c_int

def print_superblock(sb):
    __gslib2.super_block_print(sb)

def setup_superblock(
        x,y,z,vr,sec,
        nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,
        rotmat,
        radsqd,
        MAXSBX,
        MAXSBY,
        MAXSBZ):
    
    sb_p = __gslib2.setup_superblock(x,y,z,vr,sec,len(x),nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,
        rotmat,radsqd,MAXSBX,MAXSBY,MAXSBZ)
        
    return sb_p
    


if __name__ == "__main__":
    ndata = 100
    x = np.asfortranarray(np.random.random(ndata),dtype=np.float32)
    print(len(x))
    array_example(x)
    print(x[0:1])
    
    
    y = np.asfortranarray(np.random.random(ndata),dtype=np.float32)
    z = np.asfortranarray(np.random.random(ndata),dtype=np.float32)
    vr = np.asfortranarray(np.random.random(ndata),dtype=np.float32)
    sec = np.asfortranarray(np.random.random(ndata),dtype=np.float32)
    nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz = 10,20.0,40.0,10,20.0,40.0,10,20.0,40.0
    MAXSBX,MAXSBY,MAXSBZ = 10,10,10
    
    rotmat = np.asfortranarray(np.zeros((3,3)),dtype=np.float32)
    rotmat[0,0] = 1.0
    rotmat[1,1] = 1.0
    rotmat[2,2] = 1.0
    
    radsqd = 100.0
    
    sb_p = setup_superblock(x,y,z,vr,sec,
        nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,
        rotmat,
        radsqd,
        MAXSBX,
        MAXSBY,
        MAXSBZ)
        
    print_superblock(sb_p)
