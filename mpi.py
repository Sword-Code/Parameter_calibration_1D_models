try:
    from mpi4py import MPI
except ImportError:
    parallel=False
    print('WARNING: no mpi4py detected. I will proceed in serial mode.')
else:
    parallel=True
    
class Mpi:
        
    def __init__(self, comm=None):
        if parallel:
            if comm is None:
                self.comm=MPI.COMM_WORLD
            else:
                self.comm=comm
            self.rank=self.comm.rank
            self.size=self.comm.size
        else:
            self.rank=0
            self.size=1
            
        
    def gather(self, sendobj, root=0):
        if parallel:
            return self.comm.gather(sendobj, root=root)
        else:
            return [sendobj]
        
    def barrier(self):
        if parallel:
            self.comm.barrier()
            
    def bcast(self, sendobj, root=0):
        if parallel:
            return self.comm.bcast(sendobj, root=root)
        else:
            return sendobj
            
    def print(self, msg, *args, **kwargs):
        if self.rank!=0:
            return
        if 'flush' not in kwargs:
            kwargs['flush']=True
        print(msg, *args, **kwargs)
            
            
mpi=Mpi() 
