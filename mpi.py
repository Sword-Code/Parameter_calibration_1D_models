try:
    from mpi4py import MPI
except ImportError:
    parallel=False
    print('WARNING: no mpi4py detected. I will proceed in serial mode.')
else:
    parallel=True
    
class Mpi:
    
    _root=0
    
    def __init__(self, comm=None, root=0):
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
        self.root=root
    
    @property
    def root(self):
        return _root
    
    @root.setter
    def root(self, new_root):
        if not isinstance(new_root, int):
            raise TypeError(f"mpi.root must be int, but is {type(new_root)} instead")
        if new_root<0:
            raise ValueError(f"mpi.root ({new_root}) must be >= 0")
        if new_root>=self.size:
            raise ValueError(f'mpi.root ({new_root}) must be < MPI size ({mpi.size})')
        self._root=new_root
        
    def gather(self, sendobj, root=None):
        if root is None:
            root=self.root
        if parallel:
            return self.comm.gather(sendobj, root=root)
        else:
            return [sendobj]
        
    def barrier(self):
        if parallel:
            self.comm.barrier()
            
    def bcast(self, sendobj, root=None):
        if root is None:
            root=self.root
        if parallel:
            return self.comm.bcast(sendobj, root=root)
        else:
            return sendobj
            
    def print(self, msg, *args, root=None, **kwargs):
        if root is None:
            root=self.root
        if self.rank!=root:
            return
        if 'flush' not in kwargs:
            kwargs['flush']=True
        print(msg, *args, **kwargs)
            
            
mpi=Mpi() 

def main():
    node_list=mpi.gather(MPI.Get_processor_name())
    if mpi.rank!=0:
        return
    node_dict={}
    s=""
    for i, name in enumerate(node_list):
        if name in node_dict:
            node_dict[name].append(i)
        else:
            node_dict[name]=[i]
        s+=f"rank {i} of {mpi.size} is on node {name}\n"
    print("Summary:")
    for name, ranks in node_dict.items():
        print(f"node {name}: {len(ranks)} ranks ({','.join([str(i) for i in ranks])}) ")
    print()
    print("List of ranks:")
    print(s)
    
if __name__=="__main__":
    main()
