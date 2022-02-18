import numpy as np

def extend_basis(kl_basis, new_length, mode='repeat'):
    kl_basis_length = kl_basis.shape[0]
    new_basis = np.zeros((new_length, kl_basis.shape[1]))
    if new_length < kl_basis_length:
        raise Exception("The new length cannot be shorter than the current basis length!")
    if mode=='repeat':
        # in the repeat mode, we will repeat the last value of the basis until the new_length is matched
        kl_basis_last_row = kl_basis[-1,:]
        new_basis[:kl_basis_length,:] = kl_basis
        new_basis[kl_basis_length:,:] = kl_basis_last_row 
        return new_basis
    else:
        raise Exception("Mode {} is not understood".format(mode))
        